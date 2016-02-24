function [mu_h, Sigma_h] = buildHRFPrior(TR)
%_________________
% buildHRFPrior.m
%
%
%% build the HRF prior
    
% $Id: buildHRFPrior.m v0.01 2012-06-23 14:50:44 fj $

%%%% propaganda
% $$$ myLogo						= cafe_logo( mfilename) ;
    
    
    [hrf,p] = spm_hrf(TR);
    
    % 10% variation in all terms
    Sigma_p =  diag( [ 3, 6, 2, 2, 1.5, 2, 0]);
    %Sigma_p(6,6) = 1; % to adjust for the fact that p(6) = 0
    mu_h = zeros(size(hrf));
    Sigma_h = zeros([length(hrf),length(hrf)]);
    
    for it = 1:1000
        p_ = p + randn(size(p))*Sigma_p ;
        p_(1:5) = abs(p_(1:5)); % avoid negatives in p(1), p(2), p(3), p(4), p(5) 
        p_(7) = p(7);
        hrf_ = spm_hrf(TR,p_);
        mu_h = mu_h + hrf_ ;
        Sigma_h = Sigma_h + hrf_ * hrf_' ;
    end
    
    mu_h = mu_h / 1000;
    Sigma_h = Sigma_h/999 - mu_h * mu_h';
    
    %Sigma_h = Sigma_h + diag( max( abs(mu_h), 0.05));
    
end
    
