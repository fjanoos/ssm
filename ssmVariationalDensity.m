function [variational] = ssmVariationalDensity(params, data)
%_________________________
% ssmVariationalDensity.m
%
%% variational computation of factor density
    
% $Id: ssmVariationalDensity.m v0.01 2012-06-23 14:25:58 fj $
% $Id: ssmVariationalDensity.m v0.02 2012-06-27 21:50:32 fj $
    
    %%%% propaganda
    myLogo					= cafe_logo( mfilename, 'messg', [ 'subject : ' data.subject_id ' * ' datestr( now, 31)]) ;
    variational.logo				= myLogo.tmp ;

% $$$     display('computing variational density ...');
    fprintf('\ncomputing variational density ...');
    
    D						= data.D ;
    L						= data.L ;
    T						= data.T ;
    Y						= data.y ;
    s						= data.s ;
    TR						= data.TR ;
    K						= data.K ;
    tol						= data.tol ;
    
    H						= params.H ;
    invSigma_eps				= params.invSigma_eps ;
    
    variational.alpha				= zeros(T,K) ;
    % the Sigma_q matrix and its inverse
    variational.invSigma_qz			= zeros([D,D,K]) ;
    variational.Sigma_qz			= zeros([D,D,K]) ;
    variational.mu_qz				= zeros([D,K,T]) ;
    variational.q_x				= zeros(T,K) ;

    sum_H					= zeros(D,D) ;
    Delta					= zeros(D,D) ;
    for l = 1 : L
        sum_H					= sum_H + H(:,:,l) ;
        Delta					= Delta + H(:,:,l)'*invSigma_eps*H(:,:,l) ;
    end
    invDelta					= CholeskyInverse( Delta) ;
    
    % -- initialize the variational density parameters --
    invSigma_k_mu_k				= zeros(D, K) ;
    clear mu_k_invSigma_k_mu_k ln_det_Sigma_k
    for k = 1 : K
        % Sigma_qt^{-1} needs no updates       
        variational.invSigma_qz(:,:,k)		= Delta + params.invSigma_z(:,:,k) ;
        variational.Sigma_qz(:,:,k)		= CholeskyInverse( variational.invSigma_qz(:,:,k)) ;
        
        %\Sigma_k\inv\mu_k and \mu_k\Sigma_k\inv\mu_k needed during updates         
        invSigma_k_mu_k(:,k)			= params.invSigma_z(:,:,k) * params.mu_z(:,k) ;
        mu_k_invSigma_k_mu_k(k)			= params.mu_z(:,k)' * invSigma_k_mu_k(:,k) ;
        ln_det_Sigma_k(k)			= log( det( 2*pi*params.Sigma_z(:,:,k))) ;
        
        %initialize mu_q(z_t|x_t) as \Sigma_q \Sigma_k\inv\mu_k
        mu_qz_k					= variational.Sigma_qz(:,:,k) * invSigma_k_mu_k(:,k) ;
        for t = 1 : T
            variational.mu_qz(:,k,t)		= mu_qz_k ;
        end    
    end
                 
    for t = 1 : T
        % initialize q(x_t=k)
        p					= probStateTransitionMatrix(t, params, data) ; 
        p_stationary				= p^100 ; % one row gives the invariant density
        variational.q_x(t,:)			= p_stationary(1,:) ; %the invariant density
    end %T


    % use the updates to drive the variational estimation and then
    % copy into variational if change is above tolerance
    clear updates ;
    updates					= variational ; 
    %---------------------------------------------------------------------------
    % [fj] -> I observed that when using old q_t to update new q_t+1        
    % "for some reason - the alpha matrix starts to ping pong between
    % two configurations - hence test the change between iterations n and
    % n-2"
    % [fj] this above does not happen when you use dynamic updates (ie. use
    % the current estimate of q_t to update q_t+1) 
    %---------------------------------------------------------------------------
    
    %fj 20120713 progress_line( 'init') ;
	 
    %% iterate to convergence
    it_count					= 0 ;    
    while(1)        
        it_count				= it_count + 1 ;
               
        % -- alpha_qt updates --
        for t = 1 : T
             % sum_k q(x_{t-1}) ln p(x_t|x_{t-1}=k)
             p_t				=  probStateTransitionMatrix( t, params, data) ;  
             if t > 1
                updates.alpha(t,:)		= updates.q_x( t-1,:) * log( p_t+eps) ;   
             else % set q_x(0,:) as equiprobable
                updates.alpha(t,:)		=  repmat( 1.0 / single(K), [1,K]) * log(p_t+eps) ;
             end
             if t < T
                % sum_k q(x_{t+1}) ln p(x_{t+1}=k|x_t)
                p_tp1				=  probStateTransitionMatrix( t+1, params, data) ;
                updates.alpha(t,:)		= updates.alpha( t,:) + updates.q_x(t+1,:) * log(p_tp1'+eps) ;   
             end          
             if exist( 'log_partition', 'var')
                 % add the second set of terms, that act like log-partition terms
                 for k = 1 : K
                    log_partition(k)		= - 0.5*( mu_k_invSigma_k_mu_k(k) + ln_det_Sigma_k(k) ...
                            - updates.mu_qz(:,k,t)' * updates.invSigma_qz(:,:,k) * updates.mu_qz(:,k,t)) ;
                 end
                 % reduce the dynamic range of z-spread to better condition the
                 % exp operation next
                 log_partition			= log_partition - min( log_partition) ;
                 updates.alpha(t,:)		= updates.alpha( t,:) + log_partition ;
             end
             
             % compute updates.q_x(t)
             updates.q_x(t,:)			= exp( updates.alpha( t,:)) ;
             updates.q_x(t,:)			= updates.q_x(t,:)./repmat( sum( updates.q_x(t,:)), 1, K) ;
	     
	    %% progress_line( 1, T, t) ;	     
        end
        rel_change_alpha			= norm( variational.alpha - updates.alpha, 'fro') / norm( variational.alpha, 'fro') ;
        
         % -- \Sigma_qt has no updates --
         fprintf('\n') ;
	 
	 % --mu_q(z_t|x_t) updates--
         for t = 1 : T
             % compute \chi_t
             chi_t				= zeros( D,1) ;
             %\sum_m (H_m\Sigma_eps\inv(y_{t+m} - \sum_{l\neq m} H_l \nu_{t+m-l}])
             for m_idx = 1:L   %indexing variable
                 m				= m_idx - 1 ; %0-indexed value                 
                 if t+m > T
                     continue
                 end
                 
                 % keeps track of  (y_{t+m} - \sum_{l\neq m} H_l \nu_{t+m-l}])
                 acc_t_m			= Y( t+m,:)'  ; 
                 for l_idx = 1:L   % indexing variable
                      l				= l_idx-1 ;  % true value
                      if m == l 
                          continue ;
                      end
                      if t+m-l < 1
                          % replicate the first time-point backwards
                          nu_qt_ml		= zeros( size(acc_t_m)) ; % variationalMomentsOfZ(updates, 1, K ) ;
                      else
                          nu_qt_ml		= variationalMomentsOfZ( updates, t+m-l, K ) ;
                      end
                      acc_t_m			= acc_t_m - H(:,:,l_idx) * nu_qt_ml ;
                 end %l
                 chi_t				= chi_t + H(:,:,m_idx)*invSigma_eps*acc_t_m ;
             end %m
             %chi_t =invDelta*chi_t ; % this invDelta is redundanat as it cancels out later
                     
             for k=1:K     
                 % did not include Delta with chi_t as i did not premultiply it
                 % with invDelta earlier                 
                 updates.mu_qz(:,k,t)		= variational.Sigma_qz(:,:,k) * (chi_t + invSigma_k_mu_k(:,k)) ; 
             end% for k
	     %% fj20120713 progress_line( 1, T, t) ;	     
         end% for t
         
         rel_change_mu				= norm( variational.mu_qz(:) - updates.mu_qz(:)) / norm( variational.mu_qz(:)) ;
                  
         % --do the update or stop if befow tolerance--
% $$$          display (sprintf('%d: tolerances alpha %g mu %d', ...
% $$$                     [it_count, rel_change_alpha, rel_change_mu] )) ;
         fprintf(['\n\n%d: tolerances alpha %g mu %d\n'], it_count, rel_change_alpha, rel_change_mu )  ;
	 %fj 20120713 progress_line( 'line') ;
        
         variational				= updates ;
         
         % end if tolerance reached
         if  rel_change_alpha < tol && rel_change_mu < tol
             break ;
         end     
         
    end% while
    
% $$$     display('computing variational density ... done') ;
    fprintf('\ncomputing variational density ... DONE\n') ;

end %function


% %% compute E_q_t[z_t] using updated values of q
% % -- use variationalMomentsOfZ instead, with one output parameter --
% % function E_q_z = variationalExpectationOfZ( t, data, updates)
% %     E_q_z = zeros(data.D,1) ;
% %     for k = 1 : data.K
% %         E_q_z = E_q_z + updates.mu_qz(:,k,t)*updates.q_x(t,k) ;        
% %     end
% % end     
