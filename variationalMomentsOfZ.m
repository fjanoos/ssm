function [ nu_qt, Lambda_qt] = variationalMomentsOfZ( variational, t, K)
%_________________________
% variationalMomentsOfZ.m
%
% get the nu_qt and Lambda_qt matrices of eqn. E.2.12
% one output argument gives only nu_qt
    
% $Id: variationalMomentsOfZ.m v0.01 2012-06-23 15:19:13 fj $
    
    q_xt					= variational.q_x( t,:) ;
    mu_qzt					= variational.mu_qz( :,:,t) ;
    Sigma_qzt					= variational.Sigma_qz ;
    
    nu_qt					= zeros( size( mu_qzt( :,1))) ;
    Lambda_qt					= zeros( size( Sigma_qzt( :,:,1))) ;
    for k = 1 : K
        nu_qt					= nu_qt + q_xt( k) * mu_qzt( :,k) ;
        if nargout > 1
            Lambda_qt				= Lambda_qt + q_xt( k) * ( Sigma_qzt( :,:,k) + mu_qzt( :,k) * mu_qzt( :,k)') ;
        end
    end        
end
