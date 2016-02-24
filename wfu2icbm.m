function icbm_coords = wft2icbm(wfu_coords)
%____________
% wft2icbm.m
%
% transform wfu ijk (image volume) to icbm ijk coordinates

% $Id: wft2icbm.m v0.01 2012-06-23 15:21:20 fj $
    
%%%% propaganda
% $$$ myLogo						= cafe_logo( mfilename) ;



global  MAT_ICBM_IJK2XYZ DIMS LBL_VOL BRD_VOL MAT_ICBM_IJK_2_WFU_IJK ...
        MAT_WFU_IJK2_ICBM_IJK

[m,n] = size(wfu_coords);
if m == 3
    % pad a row of ones
    wfu_coords = [ wfu_coords; ones(1,n) ];
else
    error('xyz_coords wrong size');
end
[p,q] = size(MAT_WFU_IJK2_ICBM_IJK);
if p ~= 4 || q ~= 4 
    error('MAT_XYZ2IJK wrong size');
end

% flip this sign
icbm_coords = MAT_WFU_IJK2_ICBM_IJK*wfu_coords;
icbm_coords = round(icbm_coords(1:3,:));

