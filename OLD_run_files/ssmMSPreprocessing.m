function env = ssmMSPreprocessing(env)
% run pre-processing after env has been setup
% This function essentially performs two functions:
% a) build an SPM HRF design (with regressors similar to those that you will
%    use in ssm analysis) - does estimation - and computes beta-maps. These
%    beta-maps are used for initializing ssm estimation
%    Also add contrasts of interest here ... these will be extracted to create
%    spatial maps in ssmBuildStimulusMaps
% b) build an SPM FIR design (with same leading regressors as above) -- 
%    this is read and used to initialize the ssm stimulus values.

if ~strcmp( env.CURR_PLATFORM, 'numbra')
    return;
end

addpath('/home/pisti/munka/USA/Boston/BWH/ms/analysis/');

VECT.ExtrSubj={env.subject_id};

% HRF analysis
display( 'running first (HRF) stage ... ');
env.preproc_1 = preprocMSDataSet( VECT,1);

% FIR analysis (no estimate)
display( 'running second (FIR) stage ... ');
env.preproc_2 = preprocMSDataSet( VECT,2);

%%%% PISTI : 2012-04-02 - preparing motion_parameters in SPM for GIFT component sorting
display( 'setting up (ICA) stage ... ');
env.preproc_3 = preprocMSDataSet( VECT,3);

end
