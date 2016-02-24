function [ params_opt, data, results ] =  ssmEstimateFull( data, env, num_cvs, block_length, tol)
% arguments:
% -  fs_coords - the feature-space embedding of the fmri session
% -  T - the number of time points
% -  spm_file - the name of SPM mat file
% -  num_cvs - (optional, default 5) number of cross validation steps for hyper-
% parameter selection. typically set between 5-10.
% -  block_length - (optional, default 5) the length of a block of missing stim-
% uli (in tr units) for the cross validation method of hyper-parameter selection.
% typically set between 4 - 8.
% -  tol - (optional) the relative tolerance of the various parameter updates spec-
% ifying the termination criteria for iterative optimization.  typically set at 0.001
% (indicating termination for a 0.1% change in parameter value).
% returns:
% -  error_missing - an array giving the distribution of the missing stimulus
% prediction error rate (over multiple cross validations) for the optimal model.
% -  log_likelihood - an array giving the distribution of the data log likelihood
% (over multiple cross validations) for the optimal model.
% -  params_opt - a structure of the ssm parameters averaged multiple cvs for
% the optimal model. it contains the ?elds:
% -  k_opt - optimal model size (number of states).
% -  lambda_opt - optimal value of the hyper-parameter controlling the es-
% timation of the state-transition parameters.
% -  lambda_opt - optimal value of the hyper-parameter controlling the es-
% timation of the state-transition parameters.
% -  w - the 3-d array giving the state transition parameters w, where w(i,j,:)
% = w i,j .
% -  omega - the 2-d array giving the state transition parameters ?, where
% omega(i,:) = ? i .
% -  mu - An 2-d matrix where each column gives the mean value of the activa-
% tion patterns corresponding to each state.
% -  Sigma - An 3-d matrix where each 2-d sub-matrix gives the variance for
% each state.
% -  h - An 2-d matrix where each column gives the estimated hemodynamic
% response ?lter for each element of the feature (data embedding) space.
% -  Sigma_eps - A 2-d matrix giving the noise variance in the activation

% $Id: ssmEstimateFull.m v0.02 2011-03-27 13:35:00 fj $
% $Id: ssmEstimateFull.m v0.03 2012-05-29 10:01:14 fj $
% $Id: ssmEstimateFull.m v0.04 2012-06-07 10:38:19 fj $ adjusted path for temp*.mat

%%%% propaganda
myLogo						= cafe_logo( mfilename, 'messg', [ 'subject : ' data.subject_id ' * ' datestr( now, 31)]) ;
data.logo.( mfilename)				= myLogo.tmp ;	%%%% direct LOGO insert

data.tictoc.( mfilename).begin			= datestr( now, 31) ;
    
%%%% these list shall encompass the LOGOs of all relevant mfiles used for this specific SSM run
data.logo.list					= { ...
    'CholeskyInverse' , ...
    'cafe_central' , ...
    'filter_string', ...
    'progress_line' , ...
    'spmapSetup' , ...
    'ssmActivationMapFS' , ...
    'ssmBuildActivationVolume' , ...
    'ssmBuildStimulusMaps' , ...
    'ssmComputeMI' , ...
    'ssmComputeSpatialMaps' , ...
    'ssmComputeStateSpaceBalance' , ...
    'ssmEmissionParameters' , ...
    'ssmEstimateFull' , ...
    'ssmEstimateHiddenStimulus' , ...
    'ssmEstimateNoiseParameter' , ...
    'ssmEstimateOptimalSSeq' , ...
    'ssmEstimateParameters' , ...
    'ssmEstimateStateTransitionParameters' , ...
    'ssmHemodynamicParameters' , ...
    'ssmSetupEnvironment' , ...
    'ssmVariationalDensity' , ...
    'ssmVariationalDensityOptimalSSeq' , ...
    'ssmVariationalViterbi' , ...
    'xyz2ijk' } ;
myLogo						= cafe_logo( mfilename, 'messg', [ 'LOGOs for subject : ' data.subject_id ' * ' datestr( now, 31)]) ;
for i = 1 : size( data.logo.list, 2)
    myLogo					= cafe_logo( data.logo.list{i}, 'NaN') ;
    if isempty( myLogo) == 0
	fprintf(['\tok  %s\n'], data.logo.list{i}) ;
	data.logo.( data.logo.list{i})		= myLogo.tmp ;
    else
	fprintf(['\t-?- %s\n'], data.logo.list{i}) ;
	data.logo.( data.logo.list{i})		= myLogo ;
    end
end

T						= data.T ;
TR						= data.TR ;
    
do_estimate					= [] ;
    
if nargin == 2
    num_cvs					= 2 ;	%%%% PISTI <?> : does this CVS have any role ?
    block_length				= 4 ;	% should be an exact factor of T
    tol						= 5e-2 ;
end

num_blocks					= T / block_length ;
if ( num_blocks - int16( num_blocks))
    error( 'block_length must be an exact factor of T') ;
end
    
% --set up model hemodynamic hyper-parameters
[mu_h, Sigma_h]					= buildHRFPrior( TR) ;
L						= length( mu_h) ;

data.L						= L ;
data.mu_h					= mu_h ;
data.Sigma_h					= Sigma_h ;
% add eps to condition the matrix
data.invSigma_h					= CholeskyInverse( Sigma_h, eye( size( Sigma_h))) ;

% hyper-parameter tuning
pred_err					= Inf ; 
hp_opt.K					= 0 ;
hp_opt.log_lambda_w				= 0 ;
params_opt					= [] ;
results						= [] ;

%cafe_logo( mfilename, 'messg', sprintf([ ' ...computing... "%s" * %s'], data.subject_id , datestr( now, 31))) ;

if exist( 'do_estimate', 'var')
    bal_factor					= ssmComputeStateSpaceBalance( data.set) ;
% $$$     data.logo.( bal_factor.logo.name)		= bal_factor.logo ;    
    
%    for log_lambda_w = bal_factor.ll_range
     for log_lambda_w = 10:25   %% fj 20120713 -play with this range for each new data-set.   
        for K = bal_factor.K_range 
            data.K				= K ;
            data.lambda_w			= 10^log_lambda_w ;             
% $$$             display( '--------------------------------------------------------') ;
% $$$             display(sprintf('Parameter estimation for K = %d lambda = %d', data.K, log_lambda_w)) ;
% $$$             fprintf([ '\n--------------------------------------------------------\nParameter estimation for K = %d lambda = %d\n'], data.K, log_lambda_w) ;
            fprintf([ 'Parameter estimation for K = %d lambda = %d\n'], data.K, log_lambda_w) ;
	    
            % keep only 1/5th of the time points as missing. use block-bootstrap
            data.u_idx				= find( kron( randi( 5, [ num_blocks,1]), ones( block_length,1))==1) ;
            if find( data.u_idx == 1 )
                % avoid t=1 in U
                data.u_idx( find( data.u_idx == 1 )) = block_length + 1 ;
            end

            if ~exist('params_true', 'var')
                results{end+1}.params		= ssmEstimateParameters( data, env) ;
% $$$ 		data.logo.( results{end}.params.logo.name)	= results{end}.params.logo ;	%%%% END is good enough !
		
            else
                results{ end+1}.params		= ssmEstimateParameters( data, env, params_true) ;
% $$$ 		data.logo.( results{end}.params.logo.name)	= results{end}.params.logo ;	%%%% END is good enough !
            end
	    
	    %%%% returning and updating tictoc structure
	    data.tictoc				= results{ end}.params.tictoc ;
	    data.tictoc.( mfilename).log	= datestr( now, 31) ;
            
            % -- save the work so far --
% $$$	    save (fullfile( env.work_path, 'temp2.mat')) ;
            save( fullfile( env.tmp_path, 'temp2.mat')) ;

	    %cafe_logo( mfilename, 'messg', sprintf([ ' ...computing... "%s" * %s'], data.subject_id , datestr( now, 31))) ;
	    
            results{end}.K			= K ; 
            results{end}.lambda_w		= 10^log_lambda_w ;
            
            if  results{end}.params.pred_error.all_tot < pred_err
                params_opt			= results{ end}.params ;
                pred_err			= results{ end}.params.pred_error.all_tot ;
            end
% $$$             display(sprintf('DONE Parameter estimation for K = %d lambda = %d', K, log_lambda_w)) ;
            fprintf([ 'DONE Parameter estimation for K = %d lambda = %d'], K, log_lambda_w) ;
            
            break
        end
    end
end %if
    
cafe_logo( mfilename, 'messg', sprintf([ ' ...computing... "%s" * %s'], data.subject_id , datestr( now, 31))) ;

% compute optimal state sequence    
if ~exist('params_true', 'var')
    params_opt.sseq				= ssmEstimateOptimalSSeq( params_opt, data) ;
% $$$     data.logo.( params_opt.sseq.logo.name)	= params_opt.sseq.logo ;
else
    params_opt					= params_true ;
    params_opt.opt_seq				= ssmEstimateOptimalSSeq( params_opt, data, hidden) ;
% $$$     data.logo.( params_opt.opt_seq.logo.name)	= params_opt.opt_seq.logo ;
end

params_opt.Z					= data.y ;

save( fullfile( env.tmp_path, 'results.mat'), 'params_opt') ;

data.tictoc.( mfilename).end			= datestr( now, 31) ;

cafe_logo( mfilename, 'messg', sprintf([ ' ...done with :  "%s" * %s'], data.subject_id , datestr( now, 31))) ;
