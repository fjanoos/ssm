function ijk_coords = xyz2ijk( xyz_coords)
% transform xyz (real world/physical) to ijk (image volume) coordinates
% include a 4x4 tranformation matrix 
% mat_ijk2xyz: transform from ijk to xyz (opposite direction)
% xyz = 3xn input vector in physical space
% ijk = 4xn output vector in WFU image space

% NOTE: the affine matrix seems to have a flip/reflection about
% the x-axis - the coordinates from the affine transform
% do not match (or logically make sense with) those of MRIcro

% In MRICro:
% The ijk coord system is  increasing i -> right
%                                     j -> anterior
%                                     k -> superior

    
% $Id: xyz2ijk.m v0.01 2012-05-31 15:10:21 fj $

% $$$ %%%% propaganda
% $$$ eval([ 'data.logo.' mfilename ' = cafe_logo( mfilename) ;' ])
    
    
    
global   MAT_ICBM_XYZ2IJK

[ m,n]						= size( xyz_coords) ;

if m == 3
    % pad a row of ones
    xyz_coords					= [ xyz_coords ; ones( 1,n) ] ;
else
    error('xyz_coords wrong size') ;
end

[ p,q]						= size( MAT_ICBM_XYZ2IJK) ;

if p ~= 4 || q ~= 4 
    error('MAT_XYZ2IJK wrong size') ;
end

% flip this sign
ijk_coords					= MAT_ICBM_XYZ2IJK * xyz_coords ;
ijk_coords					= round( ijk_coords( 1:3,:)) ;
