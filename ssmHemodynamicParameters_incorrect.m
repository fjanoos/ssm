% function to estimate hemodynamic parameters h
% deal with each dimension of Y independently
function [H_estimate] = ...
            ssmHemodynamicParameters(data, params, variational)
        
    display('computing hemodynamic parameters ...');
    D = data.D; L = data.L; T = data.T;
    Y = data.y; s = data.s; TR = data.TR;    
    K = data.K; tol = data.tol;
    
    invSigma_eps    = params.invSigma_eps;
    Sigma_h         = data.Sigma_h ; 
    mu_h            = data.mu_h ;
        
    % initialize
    H = params.H;

    %%%% PISTI : 2012-04-03
% $$$     wb_h = waitbar(0,'computing hemodynamic parameters ...');
    for d = 1 : D             
        % we'll solve Ax = b using conjugate gradients
        A = sparse( L, L) ; %LxL matrices for the current d
        b = sparse( L, 1) ; %L vector

        for t = 1 : T
            y_t = Y(t,:); 
	    %%%% PISTI : 2012-04-03
% $$$             waitbar(((d-1)*T+t)/D/T);
                
            % construct the tilde vectors element wise (in l)
            for l_idx = 1 : L   %indexing variable
                l = l_idx - 1;  %0-indexed value 
                % corresponds to the 'l'-th entry for the dxd' dimension, at time t
                clear Lambda_t_d_dp  Lambda_t_d_d     
            
                if t-l < 1
                    continue;
                end
                
                % compute nu_q_{t-l} and Lambda_q_{t-l} as per eqn 12.2.12.
                [nu_qt_l, Lambda_qt_l] = variationalMomentsOfZ(variational, t-l, K);
                               
                for dp = 1 : D   
                    % the l-th line of the
                    %\Lambda_t[d,d'] = E[\tilde{z}_t(d) \tilde{z}_t(d')\tr] matrix
                    Lambda_t_d_dp_l = zeros(1,L);
                    Lambda_t_d_dp_l(l_idx) = Lambda_qt_l(d,dp) ;
                    for m_idx = 1 : L   %indexing variable
                        m = m_idx - 1;  %0-indexed value
                        if t-m < 1 || l == m
                            continue;
                        end 
                        % the m-th element of nu_t[d']
                        nu_qt_m = variationalMomentsOfZ(variational, t-m, K);
                        % the l-th element of nu_t[d] into the m-th element of
                        % nu_t[d']
                        Lambda_t_d_dp_l(m_idx) = nu_qt_l(d)*nu_qt_m(dp);                         
                    end
                    
                    if d == dp
                        % update A
                        A(l_idx,:) = A(l_idx,:) ...
                                    + invSigma_eps(d,d)*Lambda_t_d_dp_l;
                        b(l_idx) = b(l_idx) ...
                                        + invSigma_eps(d,d)*y_t(d)*nu_qt_l(d);
                    else
                        % b = \sum_t \{ \sum_d' \Sigma_eps\inv y_t[d'] \nu_t[d] 
                        %       - \sum_{d'\neq d}\Sigma_eps\inv \Lambda_t[d,d']h[d']     
                        h_dp = squeeze(H(dp,dp,:));
                        b(l_idx) = b(l_idx) + invSigma_eps(dp',d)*( ...
                                y_t(dp)*nu_qt_l(d) -  Lambda_t_d_dp_l*h_dp ...
                                                                  );
                    end
                end %for dp
                
            end % l
        end %t
        
        % now add in the prior contributions
        A = Sigma_h*A + speye(size(A));
        b = Sigma_h*b + mu_h;
        % note - deviating from the equations in the paper, I am
        % multiplying both sides by Sigma_h(d,d,) to avoid adding a very
        % poorly conditioned inv(Sigma_h) later.
        
        [h,flag,relres,iter,resvec]  = pcg(A,b, 1e-3, 50) ;
        display ( [ 'd = ',                             num2str(d), ...
                    ' cdg flag = ',                   num2str(flag),...
                    ' relative resdiual norm = ',     num2str(relres), ...
                    ' iterations = ',                 num2str(iter)]);
 
        H(d,d,:) = h;
        
    end %d
    
    H_estimate.H = H;
    H_estimate.rel_change = norm( H_estimate.H(:) - params.H(:))/norm( params.H(:));

    %%%% PISTI : 2012-04-03
% $$$     close(wb_h) ;
    display('computing hemodynamic parameters ... DONE');
        
end

