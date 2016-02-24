close all; clear all; clc;
% set up model estimation code

% reset the rng
stream = RandStream.getGlobalStream;
reset(stream);


subject_id = 'BG_20061205';
study_id   = '1';
study_id = 'D';
data_base = 'E:\data\ms';
code_base = 'D:\research\fMRI\scripts\statespace';


TR = 2;

rng('default'); %reset random number generator

addpath(code_base);
addpath(genpath('D:\research\matlab'));
subject_path = fullfile( data_base, ['fMRI_',subject_id]);
work_path = fullfile(subject_path, 'work');
fs_file = fullfile( fullfile(work_path, ...
            '_sub01_component_ica_s1_'), '_sub01_timecourses_ica_s1_.hdr');
        
        
subject_info = load( fullfile(subject_path,  [subject_id,'.mat']));       
fs_coord_struct = spm_vol(fs_file);   
fs_coords = fs_coord_struct.private.dat(:,:);

spm_file = fullfile(work_path, 'SPM.mat');

block_length = 4; % should be an exact factor of T

tol = 5e-2;


return;

%% garbage code
subject_info.myV.BZE_20070102.SIRPv0.D.probedigit_time(find((subject_info.myV.BZE_20070102.SIRPv0.D.load_level==5)&(subject_info.myV.BZE_20070102.SIRPv0.D.probedigit_correctans==1)'))
inv_probit(find(s_q(:,ch)==s(elem,ch),1,'first'));