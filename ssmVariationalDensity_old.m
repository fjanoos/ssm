%% variational computation of factor density
function [variational] = ssmVariationalDensity(params, data)
    
    display('computing variational density ...');
    D = data.D; L = data.L; T = data.T;
    Y = data.y; s = data.s; TR = data.TR;    
    K = data.K; tol = data.tol;
    
    H = params.H;
    
    variational.alpha = zeros(T,K);
    % the Sigma_q matrix and its inverse
    variational.invSigma_qz = zeros([D,D,K]);    
    variational.Sigma_qz = zeros([D,D,K]); 
    variational.mu_qz = zeros([D,K,T]);    
    variational.q_x = zeros(T,K);

    sum_H = zeros(D,D);
    for l = 1 : L
        sum_H = sum_H + params.H(:,:,l);
    end

    % -- initialize the variational density parameters --
    
    for k = 1 : K
        % Sigma_qt^{-1} needs no updates       
        variational.invSigma_qz(:,:,k) = params.invSigma_z(:,:,k) ...
                    +sum_H*params.invSigma_eps*sum_H;
   
        variational.Sigma_qz(:,:,k) =CholeskyInverse(...
                                variational.invSigma_qz(:,:,k), eye(D)) ;             
        
        %initialize mu_q(z_t|x_t) as \Sigma_q \Sigma_k\inv\mu_k
        for t = 1 : T
            variational.mu_qz(:,k,t) = variational.Sigma_qz(:,:,k)*...
                        params.invSigma_z(:,:,k)*params.mu_z(:,k);
        end
        
    end

    for t = 1 : T
        % initialize q(x_t=k)
        p = probStateTransitionMatrix(t, params, data); 
        p_stationary = p^100; % one row gives the invariant density
        variational.q_x(t,:) = p_stationary(1,:); %the invariant density
    end % T


    % use the updates to drive the variational estimation and then
    % copy into variational if change is above tolerance
    clear updates;
    updates = variational; 
    %---------------------------------------------------------------------------
    % [fj] -> I observed that when using old q_t to update new q_t+1        
    % "for some reason - the alpha matrix starts to ping pong between
    % two configurations - hence test the change between iterations n and
    % n-2"
    % [fj] this above does not happen when you use dynamic updates (ie. use
    % the current estimate of q_t to update q_t+1) 
    %---------------------------------------------------------------------------
    
    %% iterate to convergence
    it_count = 0;    
    while(1)        
        it_count = it_count + 1;
               
        % -- alpha_qt updates --
        for t = 1 : T
             if t > 1
                % sum_k q(x_{t-1}) ln p(x_t|x_{t-1}=k)
                p_t =  probStateTransitionMatrix(t, params, data);  
                updates.alpha(t,:) = (updates.q_x(t-1,:)*log(p_t+eps)) ;   
             end             
             if t < T
                % sum_k q(x_{t+1}) ln p(x_{t+1}=k|x_t)
                p_tp1 =  probStateTransitionMatrix(t+1, params, data);  
                updates.alpha(t,:) = updates.alpha(t,:) + updates.q_x(t+1,:)*log(p_tp1'+eps);   
             end          
             
             % compute updates.q_x(t)
             updates.q_x(t,:) =  exp(updates.alpha(t,:));
             updates.q_x(t,:) = updates.q_x(t,:)./repmat(sum(updates.q_x(t,:)), 1, K);
        end
         rel_change_alpha = norm(variational.alpha - updates.alpha, 'fro')/norm(variational.alpha, 'fro');
         
         
         % -- \Sigma_qt has no updates --
         
         % --mu_q(z_t|x_t) updates--
         for t = 1 : T
             for k=1:K     
                 updates.mu_qz(:,k,t) = params.invSigma_z(:,:,k)*params.mu_z(:,k);
                 
                 %\sum_l (H_l\Sigma_eps\inv[y_{t+l} - \sum_m H_m E[z_{t+l-m}]])
                 for l_idx = 1:L   % indexing variable
                     l = l_idx-1;  % true value
                     
                     % \sum_m H_m E[z_{t+l-m}]
                     H_m_Ez_m =  zeros(D,1); 
                     for m_idx = 1:L   %indexing variable
                        m = m_idx - 1; %0-indexed value
                        if m == l || t+l-m > T || t+l-m < 1
                            continue;
                        end  
                        H_m_Ez_m = H_m_Ez_m + ...
                            H(:,:,m_idx)*variationalExpectationOfZ(t+l-m, data, updates);
                     end %m_                     
                     
                     % y_{t+l} - \sum_m H_m E[z_{t+l-m}]
                     if t+l <= T
                         H_m_Ez_m = Y(t+l,:)' - H_m_Ez_m;
                     end
                     
                     %H_l\Sigma_eps\inv [ y_{t+l} - \sum_m H_m E[z_{t+l-m}] ]
                     updates.mu_qz(:,k,t) = updates.mu_qz(:,k,t) + ...
                            H(:,:,l_idx)*params.invSigma_eps*H_m_Ez_m;
                 end %l
                 
                 %\Sigma_q \sum_l (H_l\Sigma_eps\inv[y_{t+l} - \sum_m H_m E[z_{t+l-m}]])
                 updates.mu_qz(:,k,t) = updates.Sigma_qz(:,:,k)*updates.mu_qz(:,k,t);
                 
             end% for k
         end% for t
         
         rel_change_mu = norm(variational.mu_qz(:) - updates.mu_qz(:))/norm(variational.mu_qz(:));
                  
         % --do the update or stop if befow tolerance--
         display (sprintf('%d: tolerances alpha %g mu %d', ...
                    [it_count, rel_change_alpha, rel_change_mu] ));
        
         variational = updates;
         
         % end if tolerance reached
         if  rel_change_alpha < tol && rel_change_alpha < tol
             break;
         end     
         
    end% while

    display('computing variational density ... done');

end %function

%% compute E_q_t[z_t] using hte updated values of q
function E_q_z = variationalExpectationOfZ( t, data, updates)
    E_q_z = zeros(data.D,1);
    for k = 1 : data.K
        E_q_z = E_q_z + updates.mu_qz(:,k,t)*updates.q_x(t,k);        
    end
end       