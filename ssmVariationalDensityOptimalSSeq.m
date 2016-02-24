function [variational] = ssmVariationalDensityOptimalSSeq(params, data, sseq)
%____________________________________
% ssmVariationalDensityOptimalSSeq.m
%
%% variational computation of density for optimal state sequence estimation
%% identical to ssmVariationalDensity in z parameter

% $Id: ssmVariationalDensityOptimalSSeq.m v0.01 2012-06-23 15:15:21 fj $

%%%% propaganda
myLogo						= cafe_logo( mfilename, 'messg', [ 'subject : ' data.subject_id ' * ' datestr( now, 31)]) ;
variational.logo				= myLogo.tmp ;

% $$$     display('computing sseq variational density ...') ;
    fprintf('\ncomputing sseq variational density ...') ;
    D = data.D ; L = data.L ; T = data.T ;
    Y = data.y ; s = data.s ; TR = data.TR ;    
    K = params.K ; tol = data.tol ;
    x = sseq.x ;
    
    H						= params.H ;
    invSigma_eps				= params.invSigma_eps ;
        
    % the Sigma_q matrix and its inverse
    variational.invSigma_qz			= zeros([D,D,K]) ;    
    variational.Sigma_qz			= zeros([D,D,K]) ; 
    variational.mu_qz				= zeros([D,T]) ;    
   
    sum_H					= zeros(D,D) ;
    Delta					= zeros(D,D) ;
    for l = 1 : L
        sum_H					= sum_H + H(:,:,l) ;
        Delta					= Delta + H(:,:,l)' * invSigma_eps * H(:,:,l) ;
    end
    invDelta					= CholeskyInverse( Delta) ;
    
    % -- initialize the variational density parameters --
    invSigma_k_mu_k				= zeros( D, K) ;
    mu_qz_initial				= zeros( D,K) ;
    clear mu_k_invSigma_k_mu_k ln_det_Sigma_k
    for k = 1 : K
        % Sigma_qt^{-1} needs no updates       
        variational.invSigma_qz(:,:,k)		= Delta + params.invSigma_z(:,:,k) ;
        variational.Sigma_qz(:,:,k)		= CholeskyInverse( variational.invSigma_qz(:,:,k)) ;
        
        %\Sigma_k\inv\mu_k 
        invSigma_k_mu_k(:,k)			= params.invSigma_z(:,:,k) * params.mu_z(:,k) ;
        
        % used to initialize mu_q(z_t) as \Sigma_q \Sigma_k\inv\mu_k
        mu_qz_initial(:,k)			= variational.Sigma_qz(:,:,k)*invSigma_k_mu_k(:,k) ;
    end
    
    %initialize mu_q(z_t) as \Sigma_q \Sigma_k\inv\mu_k
    for t = 1 : T
        variational.mu_qz(:,t)			= mu_qz_initial( :,x(t) ) ;       
    end    
    
    % use the updates to drive the variational estimation and then
    % copy into variational if change is above tolerance
    clear updates ;
    updates					= variational ; 
    
    progress_line( 'init') ;
    
    %% iterate to convergence
    it_count					= 0 ;    
    while(1)        
        it_count				= it_count + 1 ;  
        
        % -- \Sigma_qt has no updates --
         
         % --mu_q(z_t|x_t) updates--
         for t = 1 : T
	     
	     progress_line( 1, T, t) ;
	     
             % compute \chi_t
             chi_t				= zeros(D,1) ;
             %\sum_m (H_m\Sigma_eps\inv(y_{t+m} - \sum_{l\neq m} H_l \nu_{t+m-l}])
             for m_idx = 1:L   %indexing variable
                 m				= m_idx - 1 ; %0-indexed value                 
                 if t+m > T
                     continue
                 end
                 
                 % keeps track of  (y_{t+m} - \sum_{l\neq m} H_l \nu_{t+m-l}])
                 acc_t_m			= Y(t+m,:)'  ; 
                 for l_idx = 1:L   % indexing variable
                      l				= l_idx-1 ;  % true value
                      if m == l 
                          continue ;
                      end
                      if t+m-l < 1
                          % replicate the first time-point backwards
                          % nu_qt_ml		= variationalMomentsOfZ(updates, 1, K ) ;
                          nu_qt_ml		= zeros(size(acc_t_m)) ; 
                      else
                          %nu_qt_ml		= variationalMomentsOfZ(updates, t+m-l, K ) ;
                          nu_qt_ml		= updates.mu_qz(:, t+m-l) ;
                      end
                      acc_t_m			= acc_t_m - H(:,:,l_idx) * nu_qt_ml ;
                 end %l
                 chi_t				= chi_t + H(:,:,m_idx) * invSigma_eps * acc_t_m ;
             end %m
             %chi_t				=invDelta*chi_t ; % this invDelta is redundanat as it cancels out later
                     
             % did not include Delta with chi_t as i did not premultiply it with invDelta earlier
             updates.mu_qz(:,t)			= variational.Sigma_qz(:,:,x(t)) * ( chi_t + invSigma_k_mu_k(:, x(t)) ) ;
         end% for t
         
         rel_change_mu				= norm(variational.mu_qz(:) - updates.mu_qz(:)) / norm( variational.mu_qz(:)) ;
                  
         % --do the update or stop if befow tolerance--
% $$$          display( sprintf('%d: tolerances  mu %d', [ it_count, rel_change_mu] )) ;
         fprintf([ '\n\n%d: tolerances  mu %d\n'], it_count, rel_change_mu) ;
	 progress_line( 'line') ;
       
         variational				= updates ;
         %the sseq from the previous optimization run is needed to determine the 
         %correct Sigma_qz at each time point
         variational.x_prev			= x ;

         % end if tolerance reached
         if  rel_change_mu < tol
             break ;
         end     
         
    end% while
    
 
% $$$     display('computing variational density ... done') ;
    fprintf('\ncomputing variational density ... DONE\n') ;
    