function [data, params_true, hidden] = ssmBuildSimulation(sim_params)
%______________________
% ssmBuildSimulation.m
%
% build a simulated data-set based on simulation parameters
% output structure has the same format at 'data' used in ssmEstimateParameters
% params_true structure has the same format as 'params_opt' returned by ""
% hidden struct contains the hidden variables x and z

% $Id: ssmBuildSimulation.m v0.01 2012-06-23 15:05:39 fj $

%%%% propaganda
myLogo						= cafe_logo( mfilename) ;

display('Building simulation run ...');
if nargin<1
    sim_params.T = 100;
    sim_params.D = 1;
    sim_params.ch = 2;
    sim_params.K = 5;
    sim_params.lambda_w = 10^-0.5;
    sim_params.TR = 2;
    sim_params.FWHM = 4;                              % stimulus smoothing in TR units
    sim_params.beta = 10;                             % inverse noise variance
    
    [sim_params.mu_h, sim_params.Sigma_h] = buildHRFPrior(sim_params.TR);    
    sim_params.L = length(sim_params.mu_h);
    
    display('using default parameters ...');
    display(sim_params);
end

T = sim_params.T;               data.T = T;
D = sim_params.D;               data.D = D;
L = sim_params.D;               data.D = D;
ch = sim_params.ch;             data.ch = ch;
K = sim_params.K;               data.K = K;
TR =sim_params.TR;              data.TR = TR;
lambda_w =sim_params.lambda_w;  data.lambda_w = lambda_w;

mu_h = sim_params.mu_h;         data.mu_h = mu_h;
Sigma_h = sim_params.Sigma_h;   data.Sigma_h = Sigma_h; 
L = sim_params.L;               data.L = L; 

block_length = 10;
num_blocks = T/block_length;
data.y = zeros(T+L,D);
data.tol = 5e-2;
data.u_idx = find(kron(randi(5, [num_blocks,1]),ones(block_length,1))==1);
if find( data.u_idx == 1 )
    % avoid t=1 in U
    data.u_idx(find( data.u_idx == 1 )) = block_length+1;
end

%% -- set up model parameters--

% state-transition params
params_true.omega = randn([ch,K]);
params_true.W = sqrt(lambda_w)*randn([ch,K,K]);         % indexed as ch, i, j

% HRF matrices
params_true.H = zeros([D,D,L]);
for d = 1 : D 
    h = mu_h + Sigma_h^0.5*randn(size(mu_h));
    params_true.H(d,d,:) = h;
end

% emmission model
params_true.mu_z = zeros( [D,K] ) ;
params_true.Sigma_z = repmat(zeros(D),[1,1,K]);
for k = 1 : K
    params_true.Sigma_z(:,:,k) = wishrnd(eye(D),T)/T;    
    % add a 40dB SNR mulitplier
    params_true.mu_z(:,k) = 100*params_true.Sigma_z(:,:,k)^0.5*randn([D,1]);
end
params_true.invSigma_z = CholeskyInverse(params_true.Sigma_z);
params_true.invSigma_z = Sparsify(params_true.invSigma_z);
params_true.Sigma_z = CholeskyInverse(params_true.invSigma_z);


params_true.Sigma_eps =  wishrnd(eye(D)/sim_params.beta,T)/T;
params_true.invSigma_eps = CholeskyInverse(params_true.Sigma_eps);
params_true.invSigma_eps = Sparsify(params_true.invSigma_eps);
params_true.Sigma_eps = CholeskyInverse(params_true.invSigma_eps);

%% -- build the data-stream --

% The stimulus sequence - using Gaussian smoothing in time
gf_ = fspecial( 'gaussian', [10,1], sim_params.FWHM);
data.s = randn( [sim_params.T, ch] );
data.s = imfilter( data.s ,gf_ ,'replicate' );
params_true.u = data.s;                                 % hidden observations

for t = 1 : T
    pm = probStateTransitionMatrix(t, params_true, data);
    if t == 1
        invar_ = pm^20;            
        md = invar_(1,:);                               % the invariant density
    else
        md = pm(hidden.x(t-1),:);
    end
    % the state sequence
    k = sampleFromMultinomial( md );
    hidden.x(t) = k;
    
    % the emission probability
    hidden.z(:,t) = params_true.mu_z(:,k) ...
                        + params_true.Sigma_z(:,:,k)^0.5*randn([D,1]);
    % hrf convolution
    for l = 0 : L-1
       data.y(t+l,:) = data.y(t+l,:) + (params_true.H(:,:,l+1)*hidden.z(:,t))';      
    end
    % noise
    data.y(t,:) = data.y(t,:) + randn([1,D])*params_true.Sigma_eps^0.5;
end

% put these guys into the global workspace
global params_true hidden data

save( 'simulation_data.mat', 'data', 'params_true', 'hidden');
end %function

% sample from a multinomial distribution
function k = sampleFromMultinomial( md )

    md_cdf = cumsum(md);
    k = find( md_cdf > rand, 1, 'first' );     %quantile function inv_cdf(y) = inf_x ( f(x) > y )

end %function
% % % test code for the mulinomial sampler
% % md_emp = zeros(5,1);
% % md_cdf = cumsum(md);
% % for n = 1 : 100000
% %     k = find( md_cdf >= rand, 1, 'first' );     %quantile function inv_cdf(y) = inf_x ( f(x) > y )
% %     md_emp(k) = md_emp(k)+1;
% % end
% % md_emp = md_emp / 100000