function ssmComputeSpatialMaps( env, params, data)
%_________________________
% ssmComputeSpatialMaps.m
%
% this function creates an estimate of Z at each time point - writes
% into a folder /ssm/K=%02d
    
% $Id: ssmComputeSpatialMaps.m v0.01 2012-05-30 17:33:38 fj $

%%%% propaganda
myLogo						= cafe_logo( mfilename, 'messg', [ 'subject : ' data.subject_id '* ' datestr( now, 31)]) ;

display('Bulding ssm spatial activity volumes (Z -layer) ... ') ;

D						= data.D ;
T						= data.T ;
Z						= params.Z ;

basis_set					= data.basis.Vbeta ;
K						= params.K ;

% a template for the other volumes
vol_info_Z_t					= spm_vol( fullfile( env.preproc_base, basis_set(1).fname)) ;

ssm_map_path					= fullfile( env.work_path, sprintf( 'K=%02d',K)) ;
if ~exist( ssm_map_path)
    mkdir ( env.work_path, sprintf('K=%02d',K)) ;
end

fprintf('\n')

% reconstruct the spatial maps from the estimated parameters the basis set
Phi						= {};
for d = 1 : D-4
    tmp						= spm_vol(fullfile( env.preproc_base, basis_set(d).fname)) ;
    Phi{d}					= tmp.private.dat(:,:,:) ;
    fprintf(':')
end

cafe_logo( mfilename, 'messg', sprintf([ ' ...computing... "%s" * %s'], env.subject_id , datestr( now, 31))) ;

for t = 1 : T
    vol_info_Z_t.fname				= fullfile( ssm_map_path, sprintf('t=%03d.nii',t)) ;
    vol_info_Z_t.descrip			= sprintf('ssm Z_t for %s with K=%02d at t = %03d', env.subject_id, K, t) ;
    rmfield( vol_info_Z_t, 'private') ;
    vol_info_Z_t				= spm_create_vol( vol_info_Z_t) ;
    vol_Z_t					= zeros( vol_info_Z_t.dim) ;
        
    for d = 1 : D-4
	vol_Z_t					= vol_Z_t+Z(t,d)*Phi{d} ;
    end
    spm_write_vol( vol_info_Z_t, vol_Z_t) ;
    fprintf('.')
end

fprintf('\n')
display( 'Bulding ssm spatial activity volumes (Z -layer) ... DONE ') ;

% $$$ end
