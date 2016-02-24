function env = ssmSetupEnvironment( env) ;
% set up the paths etc for the SSM computation
% also change ATLAS_TPL_VI and ATLAS_LBL_VI entries if needed in order
% to use for another reference brain space (or private brain space).
% regarding ICA - change the path and file name of the timecourses .hdr/img file
%   (after unzipping manually)
% 

% $Id: ssmSetupEnvironment.m v0.01 2012-04-19 00:34:00 fj $
% $Id: ssmSetupEnvironment.m v0.02 2012-05-29 09:43:22 fj $
% $Id: ssmSetupEnvironment.m v0.03 2012-06-07 20:00:01 pisti $ corrected problem with ICA filename
% $Id: ssmSetupEnvironment.m v0.04 2012-06-17 18:29:27 pisti $ switching order of lines for env.subject_base and env.data_base
    
tic
env.tictoc.start				= datestr( now, 31) ;

%%%% platforms
AVAILABLE_PLATFORMS = {'bhramand-win64', 'glenn', 'numbra'};
switch computer 
    case 'GLNXA64'
	env.CURR_PLATFORM			= AVAILABLE_PLATFORMS{3};
    case 'PCWIN64'
	env.CURR_PLATFORM			= AVAILABLE_PLATFORMS{1};
    % in case you have another platform, add its selection criterion here
end

%%%% paths
switch env.CURR_PLATFORM
    case  'bhramand-win64'
	env.fMRI_base				= 'D:\research\fMRI';
	env.data_base				= 'E:\data\ms';
    case 'numbra'
	
	%%%% set up these two paths here for your computing environment
% $$$ 	env.fMRI_base				= '/home/fjanoos/fMRI/';
% $$$	env.data_base				= '/data/kodaly/r/ms/'
	env.fMRI_base				= '/home/pisti/programs/matlab/';
	env.data_base				= '/data/kodaly/r/ms/' ;		%%%% is re-defined later !
	
end

%%%% modify these paths only if you know what you're doing...
% $$$ env.atlas_path				= fullfile( env.fMRI_base, fullfile( 'tools', 'Atlases')) ;	%%%% Atlases
% $$$ env.code_base				= fullfile( env.fMRI_base, fullfile('scripts', 'statespace')) ;	%%%% SSM
% $$$ env.spm_path				=  fullfile( env.fMRI_base, fullfile('tools', 'spm8')) ;	%%%% SPM8
env.atlas_path					= fullfile( env.fMRI_base, 'Atlases') ;				%%%% Atlases
env.code_base					= fullfile( env.fMRI_base, 'ssm') ;				%%%% SSM
env.spm_path					= fullfile( env.fMRI_base, 'spm') ;				%%%% SPM = SPM8

%%%% add the stuff
% $$$ addpath( '/home/fjanoos/programs/matlab')
% $$$ addpath('/home/pisti/programs/matlab') ;
addpath( env.fMRI_base) ;
addpath( genpath( env.code_base)) ;
addpath( env.spm_path) ;

%%%% private stuff
addpath('/home/pisti/programs/matlab/Cafe') ;

%%%% propaganda
myLogo						= cafe_logo( mfilename) ;
env.logo.( mfilename)				= myLogo.tmp ;	%%%% direct LOGO insert

%%%% loading subject specific information
%%%% PISTI 2012-04-18
% $$$ env.all_subject_info = load( fullfile( env.data_base, 'ms_subjects.mat')) ;       
% $$$ env.subject_info = env.all_subject_info.myV.( env.subject_id) ;
VECT.ExtrSubj					= { env.subject_id} ;
myV						= cafe_central( VECT) ;
env.subject_info				= myV.( env.subject_id) ;
env.data_base					= env.subject_info.db.path.data ;	%%%% adjust path !
env.subject_base				= fullfile( env.data_base, ['fMRI_', env.subject_id]) ;
env.study_id					= 'v2fj' ; 

%%%% path : work
env.work_path					= fullfile( env.subject_base, 'ssm') ;
if ~exist( env.work_path, 'file')	% 'dir')
    mkdir( env.subject_base, 'ssm') ;
end

%%%% path : tmp
env.tmp_path					= fullfile( env.work_path , 'tmp') ;
if ~exist( env.tmp_path, 'file')		% 'dir')
    mkdir( env.work_path, 'tmp') ;
end

%%%% SPM : preprocessing directories
env.preproc_base				= fullfile( env.work_path, 'preproc_hrf') ;	%%%% HRF
if ~exist( env.preproc_base, 'file')	% 'dir')
    mkdir( env.work_path, 'preproc_hrf') ;
end
env.preproc_base_2				= fullfile( env.work_path, 'preproc_fir') ;	%%%% FIR
if ~exist( env.preproc_base_2, 'file')	% 'dir')
    mkdir( env.work_path, 'preproc_fir') ;
end

%%%% SPM : preprocessing SPM.mat files
env.preproc_file				= fullfile( env.preproc_base, 'SPM.mat') ;
% $$$ env.preproc_file_2			= fullfile( env.preproc_base, 'SPM.mat') ;
%%%% PISTI : 2012-04-02
env.preproc_file_2				= fullfile( env.preproc_base_2, 'SPM.mat') ;

%%%% path : ICA/GIFT
env.ica_path					= fullfile( env.work_path, 'ica') ;
env.ica_path_zip				= fullfile( env.ica_path, env.ica_zip_fname) ;	%%%% ica_zip_fname imported

%%%% data : UNZIP here GIFT/ICA zip files
if ~exist( env.ica_path_zip, 'file')	% 'dir')
    mkdir( env.ica_path, env.ica_zip_fname) ;
end
try
    unzip( [ env.ica_path_zip,'.zip'], env.ica_path_zip ) ;
catch e
% $$$     display( 'Error in unzip - figure this out ... ') ;
    fprintf( '\nError in unzip - figure this out ... \n') ;
    display( e) ;
end

%%%% data : looking up GIFT/ICA files
k_1						= strfind( env.ica_zip_fname, '_component') ;
%%%% v0.03 2012-06-07 : ICA fname was incorrect...
% $$$ k_2					= strfind( env.ica_zip_fname, '_component') ;
k_2						= length( '_component') ;
if isempty( k_1) || isempty( k_2)
% $$$     display( 'ssmSetupEnvironment - what the hell have you specified for the ica fname') ;
    fprintf( '\nssmSetupEnvironment - what the hell have you specified for the ica fname ???\n') ;
end
%%%% v0.03 2012-06-07 : ICA fname was incorrect...
%%%% example 1 : _mean_component_ica_s_all_	_mean_timecourses_ica_s_all_.hdr
%%%% example 2 : _sub01_component_ica_s1_	_sub01_timecourses_ica_s1_.hdr
% $$$ time_course_fname				= [ env.ica_zip_fname(1:k_1), 'timecourses' , env.ica_zip_fname(k_2:end),'.hdr'] ;	% KAPUTT
time_course_fname				= [ env.ica_zip_fname(1:k_1), 'timecourses' , env.ica_zip_fname(k_1+k_2:end),'.hdr'] ;	% GOOD
env.fs_file					= fullfile( env.ica_path_zip, time_course_fname) ;

%%%% PISTI : 2012-04-02 - SPM.mat for ICA
env.preproc_file_3				= fullfile( env.ica_path, 'SPM.mat') ;

%%%% SSM : random number generator
ssmResetRNG( env) ;

% -- for normalization to atlas space -- typcially ICBM templates
global MAT_ICBM_IJK2XYZ DIMS LBL_VOL BRD_VOL MAT_ICBM_IJK_2_WFU_IJK MAT_WFU_IJK2_ICBM_IJK MAT_ICBM_XYZ2IJK ATLAS_TPL_VI

% Change these two names below if you want to use another "reference" structural
% image
% The reference brain (normalized space) 
ATLAS_TPL_VI					= spm_vol( fullfile( env.atlas_path, 'ICBM_Template.2mm.hdr')) ;
% The labels (grey / white / csf) for the above brain
ATLAS_LBL_VI					= spm_vol( fullfile( env.atlas_path, 'ICBM_Template_Labels.2mm.hdr')) ;
% internal structures -- do not modify
MAT_ICBM_IJK2XYZ				= ATLAS_LBL_VI.mat ;		% MNI to ijk
% some fuckall inversion here
% MAT_ICBM_IJK2XYZ(1,:)				= -MAT_ICBM_IJK2XYZ(1,:) ; 
MAT_ICBM_XYZ2IJK				= inv(MAT_ICBM_IJK2XYZ) ;	% MNI to ijk
% MAT_ICBM_XYZ2IJK(1,1)				= -MAT_ICBM_XYZ2IJK(1,1) ; 
% WFU coordinates
ATLAS_BRD_VI					= spm_vol( fullfile( env.atlas_path, 'dilate_1.nii')) ;
MAT_WFU_IJK2XYZ					= ATLAS_BRD_VI.mat ;		% MNI to ijk
% MAT_WFU_IJK2XYZ(1,:)				= -MAT_WFU_IJK2XYZ(1,:) ; 
MAT_WFU_XYZ2IJK					= inv(MAT_WFU_IJK2XYZ) ;	% MNI to ijk
% MAT_WFU_XYZ2IJK(1,1)				= -MAT_WFU_XYZ2IJK(1,1) ; 

% convert the ICBM space to the WFU space
MAT_WFU_IJK2_ICBM_IJK				= MAT_ICBM_XYZ2IJK*MAT_WFU_IJK2XYZ ;
MAT_ICBM_IJK_2_WFU_IJK				= MAT_WFU_XYZ2IJK*MAT_ICBM_IJK2XYZ ;

load ( fullfile( env.atlas_path, 'labels.mat')) ;				% ATLAS_LBL_NAMES

% labels to keep
KEEP_LIST					= ATLAS_LBL_NAMES{1}([1:48,53,54,56,57] ) ;

% output volume
DIMS						= ATLAS_TPL_VI.dim; 
BRD_VOL						= ATLAS_BRD_VI.private.dat(:,:,:) ; 	% WFU
lbl_vol						= ATLAS_LBL_VI.private.dat(:,:,:) ; 	% ICBM
LBL_VOL						= zeros( DIMS) ;			% ICBM coords
% ATLAS_VOL					= ATLAS_TPL_VI.private.dat(:,:,:) ;	% ICBM
for lbl = 1 : length([KEEP_LIST])
    LBL_VOL( find( lbl_vol == KEEP_LIST( lbl)))	= lbl;
 end

% new_vi = ATLAS_LBL_VI;
% new_vi.private = [];
% new_vi.fname = fullfile(atlas_path, 'mask.2mm.nii') ;
% new_vi = spm_create_vol(new_vi) ;
% new_vi = spm_write_vol(new_vi, ATLAS_VOL) ;
