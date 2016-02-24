function ssmBuildStimulusMaps( env, params, condition, rf)
%________________________
% ssmBuildStimulusMaps.m
%
% convert the params structure into an activation volume.
%   
% the 4-th parameter is a relaxation factor that increases the amount
% of visible activity ... similar to opening up the p-value in spm

% $Id: ssmBuildStimulusMaps.m v0.01 2012-05-30 17:53:10 fj $

%%%% propaganda
myLogo						= cafe_logo( mfilename, 'messg', [ 'subject : ' params.subject_id ' * ' datestr( now, 31)]) ;
    
global  ATLAS_TPL_VI 

% the activation map information
sptmap						= spmapSetup( env, condition, params) ;

% in feature space coordimates
% rf is a "relaxation factor" - a threshold similar to opening up the p value in SPM
fs_map						= ssmActivationMapFS( sptmap , rf) ;

num_maxima					= size( fs_map.dat,1) ;

for m = 1 :  num_maxima

    fprintf('.')
    
    if ~isempty( fs_map.dat{ m,5})
        sz_cluster( m)				= fs_map.dat{ m,5} ;
    else
        sz_cluster( m-1)			= sz_cluster( m-1)/2 ;
        sz_cluster( m)				= sz_cluster( m-1) ;
    end
    
    val_peak( m)				= fs_map.dat{ m,9} ;

    ijk_coords(:,m)				= xyz2ijk( fs_map.dat{ m,12}) ;

    rd_cluster( m)				= ( sz_cluster( m)/pi)^(1/3) ;
end   

for m = 1 : num_maxima
    if any ( ijk_coords(:,m) <= 0 )
        val_peak( m)				= 0 ;
        ijk_coords(:,m)				= [1;1;1] ;
        rd_cluster( m)				= 0 ;
    end
end   

act_vol						= ssmBuildActivationVolume( env, ijk_coords, val_peak, ceil(rd_cluster), 1) ;

act_vi						= ATLAS_TPL_VI ;
act_vi.dt					= [ 16,0] ;
act_vi.private					= [] ;
act_vi.fname					= fullfile( env.work_path,[ condition,'.nii']) ;
act_vi.descrip					= condition ;
act_vi						= spm_create_vol( act_vi) ;
act_vi						= spm_write_vol( act_vi, act_vol) ;
