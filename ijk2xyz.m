function xyz_coords = ijk2xyz(ijk_coords)
% transform icbm ijk (image volume) coordinates to xyz (real world/physical)%
% xyz = 3xn  vector in physical space
% ijk = 3xn  vector in ICBM coordinate image space

% NOTE: the affine matrix seems to have a flip/reflection about
% the x-axis - the coordinates from the affine transform
% do not match (or logically make sense with) those of MRIcro

                                    k -> superior
global  ATLAS_TPL_VI ATLAS_LBL_VI ATLAS_LBL_NAMES BRD_VOL ...
        LBL_VOL ...
        AVAILABLE_PLATFORMS CURR_PLATFORM ...
        numbra_path code_path spm_path atlas_path results_path ...
        MAT_WFU_IJK2_ICBM_IJK MAT_WFU_IJK2XYZ MAT_WFU_XYZ2IJK ...
        MAT_ICBM_IJK2XYZ MAT_ICBM_XYZ2IJK DIMS ...
        DEF_RAND_STREAM rand_stream_state

[m,n] = size(ijk_coords);
if m == 3
    % pad a row of ones
    ijk_coords = [ ijk_coords; ones(1,n) ];
else
    error('xyz_coords wrong size');
end
[p,q] = size(MAT_ICBM_IJK2XYZ);
if p ~= 4 || q ~= 4 
    error('MAT_IJK2XYZ wrong size');
end

xyz_coords = MAT_ICBM_IJK2XYZ*ijk_coords;
xyz_coords = round(xyz_coords(1:3,:));