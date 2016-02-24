function pm = probStateTransitionMatrix(t, params, data)
%_____________________________
% probStateTransitionMatrix.m
%
% transition from state p[t-1,t] = p(x_t|x_{t-1})
% returns the transition matrix
    
% $Id: probStateTransitionMatrix.m v0.01 2012-06-23 15:01:51 fj $

    K = data.K; T = data.T;
    if isempty( find( data.u_idx == t ) ) 
        s_t = data.s(t,:);
    else
        if isfield( params, 'u') 
            s_t = params.u(t,:); % params.u(t,:);
        else
            s_t = data.s(t,:);
        end
    end
    
    pm = zeros(K,K);
    log_pm = zeros(K,K);
    for j = 1 : K
        omega_j = params.omega(:,j);
        for i = 1: K
            w_ij = params.W(:,i,j);
            pm(i,j) = exp(s_t*(omega_j + w_ij));
            log_pm(i,j) = s_t*(omega_j + w_ij);
        end
    end
    log_pm = log_pm - repmat(max(log_pm,[],2), 1, K);
    pm = exp(log_pm);
    %normalize
    pm = pm ./ repmat( sum(pm,2), 1, K);
end

