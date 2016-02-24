function spmap = spmapSetup( env, condition, params)
%______________
% spmapSetup.m
%
% the ref condition is obtained by saving spmGetResults for a
% reference subject from spm_results_ui - along with the contrast
% number for the condition

% $Id: spmapSetup.m v0.01 2012-05-31 15:08:12 fj $
% $Id: spmapSetup.m v0.02 2012-05-31 15:40:15 pisti $ changed path to ref-*.mat files
    
%%%% propaganda
myLogo						= cafe_logo( mfilename) ;

% $$$ v0.02 2012-05-31 : changing path for ref-*.mat files - keep them now in the SPM HRF directory
% $$$ spmap					= load ( fullfile( env.data_base, [ 'ref-',condition,'.mat'] ) ) ;
spmap						= load ( fullfile( env.preproc_base, [ 'ref-',condition,'.mat'] ) ) ;

spmap.fname					= fullfile( env.preproc_base, sprintf( 'spmT_%04d.img', spmap.cond_number) ) ;
spmap.xSPM.swd					= env.preproc_base ;
spmap.xSPM.title				= condition ;
spmap.xSPM.Vspm					= spm_vol( spmap.fname ) ;

% $$$ end
