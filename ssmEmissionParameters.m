function [vartheta] = ssmEmissionParameters(data, params, variational)
%_________________________
% ssmEmissionParameters.m
%
% function to estimate emmission parameters mu_k and Sigma_k
    
% $Id: ssmEmissionParameters.m v0.01 2012-06-23 14:32:13 fj $
    
    %%%% propaganda
    myLogo					= cafe_logo( mfilename, 'messg', [ 'subject : ' data.subject_id ' * ' datestr( now, 31)]) ;
    vartheta.logo				= myLogo.tmp ;
    
    D = data.D ; L = data.L ; T = data.T ;
    Y = data.y ; s = data.s ; TR = data.TR ;    
    K = data.K ; tol = data.tol ;      
    
    vartheta.mu_z				= zeros(size(params.mu_z)) ;
    vartheta.Sigma_z				= zeros(size(params.Sigma_z)) ;
    vartheta.invSigma_z				= zeros(size(params.invSigma_z)) ;
    q_x						= variational.q_x ;
    mu_qz					= variational.mu_qz ;
    Sigma_qz					= variational.Sigma_qz ;
    
    for k = 1 : K
        den					= 0 ;    % the effective number of time-points in state k
        for t = 1 : T
            vartheta.mu_z(:,k)			= vartheta.mu_z(:,k) + q_x(t,k)*mu_qz(:,k,t) ;            
            den					= den +  q_x(t,k) ;
        end %t
        vartheta.mu_z(:,k)			= vartheta.mu_z(:,k)/den ;
        for t = 1 : T
            vartheta.Sigma_z(:,:,k)		= vartheta.Sigma_z(:,:,k) ...
                + q_x(t,k)*( ...
		    Sigma_qz(:,:,k) + mu_qz(:,k,t)*mu_qz(:,k,t)' ...
                    - vartheta.mu_z(:,k)*mu_qz(:,k,t)' ...
                    - mu_qz(:,k,t)*vartheta.mu_z(:,k)' ...
                    + vartheta.mu_z(:,k)*vartheta.mu_z(:,k)' ...
                    ) ;
        end %t
        vartheta.Sigma_z(:,:,k)			= vartheta.Sigma_z(:,:,k)/den ;
       
        vartheta.rel_change_mu(k)		= norm(vartheta.mu_z(:) - params.mu_z(:)) / norm(params.mu_z(:)) ;
        vartheta.rel_change_Sigma(k)		= norm(vartheta.Sigma_z(:) - params.Sigma_z(:))/norm(params.Sigma_z(:)) ;
                
        if vartheta.Sigma_z(:,:,k) ~= vartheta.Sigma_z(:,:,k)'
            error( 'ssmEmissionParameters : vartheta.Sigma_z not symmetric') ;
        end
    end %k    
    vartheta.invSigma_z				= CholeskyInverse( vartheta.Sigma_z) ;
    vartheta.invSigma_z				= Sparsify( vartheta.invSigma_z) ;
    vartheta.Sigma_z				= CholeskyInverse( vartheta.invSigma_z) ;

%end
