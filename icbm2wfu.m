function wfu_coords = icbm2wfu(icbm_coords)
%____________
% icbm2wfu.m
%
% transform icbm ijk (image volume) to wfu ijk coordinates

% $Id: icbm2wfu.m v0.01 2012-06-23 14:54:25 fj $

%%%% propaganda
% $$$ myLogo						= cafe_logo( mfilename) ;

global  MAT_ICBM_IJK2XYZ DIMS LBL_VOL BRD_VOL MAT_ICBM_IJK_2_WFU_IJK ...
        MAT_WFU_IJK2_ICBM_IJK

[m,n] = size(icbm_coords);
if m == 3
    % pad a row of ones
    icbm_coords = [ icbm_coords; ones(1,n) ];
else
    error('xyz_coords wrong size');
end
[p,q] = size(MAT_ICBM_IJK_2_WFU_IJK);
if p ~= 4 || q ~= 4 
    error('MAT_XYZ2IJK wrong size');
end

% flip this sign
wfu_coords = MAT_ICBM_IJK_2_WFU_IJK*icbm_coords;
wfu_coords = round(wfu_coords(1:3,:));

