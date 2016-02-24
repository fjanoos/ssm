function [sseq] = ssmVariationalViterbi( data, params, variational)
%_________________________
% ssmVariationalViterbi.m
%
    
% $Id: ssmVariationalViterbi.m v0.01 2012-06-23 15:17:35 fj $
    
%%%% propaganda
myLogo						= cafe_logo( mfilename) ;
sseq.logo					= myLogo.tmp ;

% $$$     display( 'Viterbi ...');
    fprintf( '\nViterbi ...');
    
    progress_line( 'init') ;
    
    D = data.D;   T = data.T;  K = data.K; 

    log_p_opt = zeros(T,K);     % matrix of optimal sequence probabilities
    psi_opt   = zeros(T,K);     % matrix of optimal state (backtrack)

    % initialize t = 1 with the the invariant density
    pm = probStateTransitionMatrix(1, params, data); 
    ps = pm^100; % one row gives the invariant density
    log_p_opt(1,:) = log(ps(1,:)); 
    [d_,i_] = max(log_p_opt(1,:));
    psi_opt(1,:) = repmat( i_, [1,K] );

    %% The forward pass of the viterbi algorithm return the matrices
    for t = 2 : T
	
	progress_line( 2, T, t) ;
	
        % compute [ln p(x_t|x_t-1) + int q(z_t)\ln p(z_t|x_t)dz_t   
        log_p_t = computeLogProb(t,data, params, variational);
        % for state k at time t
        for k = 1: K
            % the row-index of the log_pt(:,k) gives the state at time point t-1
            [log_p_opt(t,k), psi_opt(t,k)] = max( log_p_t(:,k) + log_p_opt(t-1,:)' );
        end    
    end

    %% Backtracking
    sseq.x = zeros([T,1]);
    [sseq.log_p_opt, sseq.x(T)] = max( log_p_opt(T,:) ) ;

    for t = T-1 :-1: 1
        % read out the psi_opt matrix in backward ordering
        sseq.x(t) =  psi_opt(t, sseq.x(t+1));
    end
    
    progress_line( 'term') ;
    
% $$$     display( 'Viterbi ... done');
    fprintf( '\nViterbi ... DONE\n');
end % function

%% 
function lp = computeLogProb(t, data, params, variational)
    K = data.K; 
    Sigma_z = params.Sigma_z;
    invSigma_z = params.invSigma_z;
    mu_z = params.mu_z;
    
    mu_qzt = variational.mu_qz(:,t);
    Sigma_qzt = variational.Sigma_qz(:,:,variational.x_prev(t)); 
    Sigma_mu_qzt = Sigma_qzt + mu_qzt*mu_qzt'; %temp store
     
    % compute \ln p(x_t|x_{t-1}) + \int q(z_t)\ln p(z_t|x_t)dz_t 
    log_pm = log(probStateTransitionMatrix(t, params, data));

    log_p_emit = zeros(1,K); % keep track of the second term

    for k = 1 : K
        log_p_emit(k) =  -0.5*log(det(2*pi*Sigma_z(:,:,k))) ...
                        - 0.5*mu_z(:,k)'*invSigma_z(:,:,k)*mu_z(:,k) ...
                        + mu_z(:,k)'*invSigma_z(:,:,k)*mu_qzt ...
                        - 0.5*trace( invSigma_z(:,:,k)*Sigma_mu_qzt );
    end
    
    % add in the emmision probability to the state-transition prob matrix
    lp = log_pm + repmat( log_p_emit, [K,1]) ; 
end
