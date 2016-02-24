close all; clear all; clc;
% set up model estimation code
tic;

% $$$ env.subject_id = 'LH_20070608' ;	% complete
% $$$ env.subject_id = 'BG_20061205' ;	% sibelius

% $$$ env.subject_id = 'AJ_20061121' ;	% rachmaninoff ; enumerate
% $$$ env.subject_id = 'AR_20070227' ;	% enumerate
% $$$ env.subject_id = 'BC_20070518' ;	% enumerate
% $$$ env.subject_id = 'BK_20070613' ;	% enumerate
% $$$ env.subject_id = 'BZE_20070102' ;	% rachmaninoff ; enumerate
% $$$ env.subject_id = 'CM_20070511' ; % rachmaninoff
% $$$ env.subject_id = 'CZ_20070524' ; % enumerate
% $$$  env.subject_id = 'DJ_20061219' ; % enumerate
% $$$  env.subject_id = 'DR_20061128' ; % rachmaninoff
% $$$  env.subject_id = 'EM_20061124' ; % enumerate
% $$$  env.subject_id = 'FC_20070223' ; % rachmaninoff
% $$$  env.subject_id = 'MP_20070601' ; % enumerate
% $$$  env.subject_id = 'NP_20061212' ; % rachmaninoff
% $$$  env.subject_id = 'SD_20070411' ; % enumerate ?
% $$$  env.subject_id = 'SK_20061215' ; % rachmaninoff
% $$$  env.subject_id = 'SL_20070614' ; % enumerate
% $$$  env.subject_id = 'SS_20070504' ; % enumerate
 env.subject_id = 'SW_20070530' ; % rachmaninoff

%The zip file greated by GIFT - no zip extension !
env.ica_zip_fname = '_mean_component_ica_s_all_';  

% -- these have to be set for each subject --
% p1 controls the balance between the previous state versus the current stimulus
% on the outcome of the current state. Increasing p1 increases the effect of
% the current stimulus on the current state - which generally results in 
% greater variance 
% p2 controls the balance between prediction error and data (fMRI signal) fitting. 
% A lower p2 makes the model optimize for prediction error at the cost of
% explaining (fitting) the fMRI data. Increasing p2 generally increases the
% desire of the model to explain the fMRI data as good as possible - which
% also increases the variance in the number of state per subject (because
% more states generally results in better fit.
set.p1 = 2;     
set.p2 = 15; 

% set up the env structure needed to do estimation
env = ssmSetupEnvironment(env)

% Run preprocessing if necessary
if ~exist(env.preproc_file, 'file') || ~exist(env.preproc_file_2, 'file')
    env = ssmMSPreprocessing(env); 
end


% the lenght of a block for predcition. typically 25% of all TRs treated
% as hidden (to be predicted)
block_length = 4; % should be an exact factor of T 

% decrease tolerance to refine results
tol = 5e-2;

% load the design information from the FIR design file.    
basis_struct = load(env.preproc_file_2); 
s = basis_struct.SPM.xX.X(:,1:36);

% number of time-points
T = sum(basis_struct.SPM.nscan); % the sum is in case of multiple runs
env.TR = str2num(basis_struct.SPM.xsDes.Interscan_interval(1:end-3));
if exist(env.fs_file, 'file')        
    fs_coord_struct = spm_vol(env.fs_file);   
    fs_coords = fs_coord_struct.private.dat(:,:);
else
    b_st = load(env.preproc_file); 
    fs_coords = b_st.SPM.xX.X;;
end

D = size(fs_coords,2);
if size(fs_coords,1)~=T
    error('Time-series dimensions not correct');
end

% transform the stimulus into [-1,+1] range though the
%inverse probit transform    
clear s_norm;
s_norm = zeros(size(s));
quantile_values  = linspace(0.01, 0.99,10 );
try
    inv_probit = icdf( 'norm',quantile_values, 0, 1);    
    s_q = quantile( s, quantile_values);
     % iterate thru all stimulus channels
    for ch = 1:size(s,2)
        for elem = 1 : size(s,1)
            qntl = find(s_q(:,ch)<=s(elem,ch),1,'last');
           try
            s_norm(elem,ch) =  inv_probit(qntl);
           catch e
              s_norm(elem,ch) = s(elem,ch);
           end
        end
    end   
catch e   
    s_norm = s;
    %inv_probit = ltqnorm(quantile_values);    
end


data.D = D;
data.T = T;
data.TR = env.TR;
data.s = s_norm;
data.ch =size(s,2); %number of channels
data.y = fs_coords;
data.tol = tol;
data.basis = b_st.SPM;
data.set = set;
[params_opt, data, results]=  ssmEstimateFull(data, env, 2, block_length, tol);
toc;

save (fullfile( env.work_path, 'ssm_output.mat'), ...
                            'params_opt', 'data', 'results', 'env');

display('----------------------------------------------------------------------');                        
display( 'computing spatial maps');                        
display('----------------------------------------------------------------------');   

% this function creates an estimate of Z at each time point.
% writes into a folder /ssm/K=%02d
ssmComputeSpatialMaps(env, params_opt, data); 
% these functions create maps corresponding to particular "conditions" (i.e. spm
% contrasts)
ssmBuildStimulusMaps(env, params_opt, 'encode-all', 2.5);
ssmBuildStimulusMaps(env, params_opt, 'target-all', 3);
ssmBuildStimulusMaps(env, params_opt, 'foil-all', 3);
ssmBuildStimulusMaps(env, params_opt, 'encode-load-1', 2.5);
ssmBuildStimulusMaps(env, params_opt, 'encode-load-3', 2.5);
ssmBuildStimulusMaps(env, params_opt, 'encode-load-5', 2.5);
ssmBuildStimulusMaps(env, params_opt, 'target-load-1',4);
ssmBuildStimulusMaps(env, params_opt, 'target-load-3',4);
ssmBuildStimulusMaps(env, params_opt, 'target-load-5',4);
ssmBuildStimulusMaps(env, params_opt, 'foil-load-1',4);
ssmBuildStimulusMaps(env, params_opt, 'foil-load-3',4);
ssmBuildStimulusMaps(env, params_opt, 'foil-load-5',4);

return;

% not included here is ssmComputeMI. Run this file after you have estimated
% the models for multiple subjects.

% also - you have write specific small scripts to query / analyze results.


