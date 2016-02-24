function act_vol = ssmBuildActivationVolume( env, ijk_coords_icbm, values, radius_in_vox, bck_noise)
% $Id: ssmBuildActivationVolume.m v0.01 2012-05-30 18:36:02 fj $

%%%% propaganda
eval([ 'env.logo.' mfilename ' = cafe_logo( mfilename) ;' ])
    
% ijk_coords : arranged column wise (in icbm space)
% size: size of focus in voxels (gaussian)
% values: average value 

global  MAT_ICBM_IJK2XYZ DIMS LBL_VOL BRD_VOL

global ri rand_nums ang_tol rot_mat

% % DEF_RAND_STREAM.State = rand_stream_state ;
ssmResetRNG(env) ;

ang_tol						= 10 / 180 * pi ;

rand_nums					= rand( 10000, 1) ;
ri						= 1 ; 

% first add background noise
noise_vol					= gaussianSmooth( bck_noise * rand(DIMS) ,3)+ 0.01*bck_noise * rand(DIMS) ;
noise_vol					= noise_vol .* ( LBL_VOL > 0 ) ;

n_foci						= size( ijk_coords_icbm,2) ;

%also get icbm coordinates
ijk_coords_wfu					= icbm2wfu( ijk_coords_icbm) ;

act_vol_smooth_foci				= zeros( DIMS) ;
act_vol_noisy_foci 				= zeros( DIMS) ;
act_vol_noisy_smooth_foci 			= zeros( DIMS) ;

%%%% PISTI : 2012-04-03
% $$$ h = waitbar(0, 'Building activation map ') ;
for fc = 1 : n_foci   
    %%%% PISTI : 2012-04-03
% $$$     waitbar( fc /n_foci, h) ;
    fprintf('.')
    
    % focus coordinates
    fc_ijk_wfu					= ijk_coords_wfu( :,fc) ;
    fc_ijk_icbm					= ijk_coords_icbm( :,fc) ; %icbm coords
    
    for d = 1 : 3
        if fc_ijk_icbm(d) < 1 ;		fc_ijk_icbm(d) = 1 ; end
        if fc_ijk_icbm(d) > DIMS(d) ;	fc_ijk_icbm(d) = DIMS(d) ; end        
    end
        
    % focus radius
    fc_rad_vox					= radius_in_vox( fc) ;
    if fc_rad_vox <= 0 ; fc_rad_vox = 1 ; end
    % focus value
    fc_value					= values( fc) ;
    %the focus label and brodmann number
    fc_lbl					= LBL_VOL( fc_ijk_icbm(1), fc_ijk_icbm(2), fc_ijk_icbm(3)) ;
    fc_brd					= BRD_VOL( fc_ijk_wfu(1), fc_ijk_wfu(2), fc_ijk_wfu(3)) ;
    
      
    %wfu coordinates
    act_vol_noisy_smooth_foci(fc_ijk_icbm(1),fc_ijk_icbm(2),fc_ijk_icbm(3)) = ...
        act_vol_noisy_smooth_foci(fc_ijk_icbm(1),fc_ijk_icbm(2),fc_ijk_icbm(3)) + fc_value  ;
    act_vol_noisy_foci(fc_ijk_icbm(1),fc_ijk_icbm(2),fc_ijk_icbm(3)) = ...
        act_vol_noisy_foci(fc_ijk_icbm(1),fc_ijk_icbm(2),fc_ijk_icbm(3)) + fc_value  ;
    
                        
    %preferential angles to keep     
    ri						= mod( ri + 11, 9900) + 1 ; 
    fc_angles_azim				= 2 * pi * rand_nums( ri:ri+4) - pi ;
    ri						= mod( ri + 11, 9900) + 1 ; 
    fc_angles_elev				= pi * rand_nums( ri:ri+4) - pi/2 ;
    ri						= mod( ri + 11, 9995) + 1 ; 
    
    cos_azim					= cos( fc_angles_azim) ;
    sin_azim					= sin( fc_angles_azim) ;
    cos_elev					= cos( fc_angles_elev) ;
    sin_elev 					= sin( fc_angles_elev) ;
    
    
    % rotation matrices
    for rm = 1 : 4 
        rot_mat(:,:,rm)				= ...
	    [ cos_azim(rm) sin_azim(rm) 0 ; -sin_azim(rm) cos_azim(rm) 0 ; 0 0 1 ] * ...
	    [ cos_elev(rm) 0 sin_elev(rm) ; -sin_elev(rm) 0  cos_elev(rm) ; 0 1 0 ] ;
    end
    
                        
% $$$     for i = floor(-3*(fc_rad_vox+4*rand_nums(ri))):ceil(+3*(fc_rad_vox+4*rand_nums(ri+1)))
    for i = floor(-4*(fc_rad_vox)):ceil(+4*(fc_rad_vox))
% $$$ 	ri = mod( ri + 2, 9995) + 1 ; 
	if i+fc_ijk_wfu(1) < 1 || i+fc_ijk_wfu(1) > DIMS(1) ...
		|| i+fc_ijk_icbm(1) < 1 || i+fc_ijk_icbm(1) > DIMS(1) 
	    continue  
	end
        
% $$$ 	for j = floor(-3*(fc_rad_vox+4*rand_nums(ri))):ceil(+3*(fc_rad_vox+4*rand_nums(ri+1)))
	for j = floor(-4*(fc_rad_vox)):ceil(+4*(fc_rad_vox))
% $$$ 	    ri = mod( ri + 2, 9995) + 1 ; 
            if j+fc_ijk_wfu(2) < 1 || j+fc_ijk_wfu(2) > DIMS(2) ...
		    ||  j+fc_ijk_icbm(2) < 1 || j+fc_ijk_icbm(2) > DIMS(2)
		continue ;  
            end
	    
% $$$ 	    for k = floor(-3*(fc_rad_vox+4*rand_nums(ri))):ceil(+3*(fc_rad_vox+4*rand_nums(ri+1)))
            for k = floor(-4*(fc_rad_vox)):ceil(+4*(fc_rad_vox))
% $$$ 		ri = mod( ri + 2, 9995) + 1 ; 
                if k+fc_ijk_wfu(3) < 1 || k+fc_ijk_wfu(3) > DIMS(3)...
			|| k+fc_ijk_icbm(3) < 1 || k+fc_ijk_icbm(3) > DIMS(3)
                    continue ; 
                end               
                
                test_flg			= ...
		    ( LBL_VOL( fc_ijk_icbm(1)+i, ...
			       fc_ijk_icbm(2)+j, ...
			       fc_ijk_icbm(3)+k ) == fc_lbl ) & ...
		    ( BRD_VOL( fc_ijk_wfu(1)+i,...
			       fc_ijk_wfu(2)+j,...
			       fc_ijk_wfu(3)+k ) == fc_brd ) ;
		
                pv				= buildPattern( i,j,k, fc_rad_vox, fc_angles_azim, fc_angles_elev, test_flg ) ;
                
		act_vol_noisy_smooth_foci(fc_ijk_icbm(1)+i,fc_ijk_icbm(2)+j,fc_ijk_icbm(3)+k) = ...
		    act_vol_noisy_smooth_foci(fc_ijk_icbm(1)+i,fc_ijk_icbm(2)+ j,fc_ijk_icbm(3)+k) + fc_value*pv ;
		
		act_vol_noisy_foci(fc_ijk_icbm(1)+i,fc_ijk_icbm(2)+j,fc_ijk_icbm(3)+k) = ...
		    act_vol_noisy_foci(fc_ijk_icbm(1)+i,fc_ijk_icbm(2)+ j,fc_ijk_icbm(3)+k) + fc_value*pv ;
            end % for k
        end %for j
    end %for i
end %for fc

%%%% PISTI : 2012-04-03
% $$$ close(h) % waitbar
fprintf('\n\n')

% smooth out act_vol_smooth_foci
act_vol_noisy_smooth_foci			= 3 * gaussianSmooth( act_vol_noisy_smooth_foci ,0.5 ) ;
act_vol						= noise_vol + ( act_vol_noisy_smooth_foci ).* ( LBL_VOL > 0 ) ;

end % function


% --------------------------------------------------------------------------------
% helper functions
function vol_out = gaussianSmooth( vol_in, RAD )
% filter size
%RAD = 3 ; %in voxels
    x						= [ -3*RAD:3*RAD] ;
    g						= exp(-x.^2 /2/ RAD^2 ) / sqrt(2*pi) / RAD ;      
    
    vol_out					= shiftdim( conv1d( g, vol_in),1) ;
    vol_out					= shiftdim( conv1d( g, vol_out),1) ;
    vol_out					= shiftdim( conv1d( g, vol_out),1) ;
end %function   

% --------------------------------------------------------------------------------
function vol_out = conv1d( h_row, vol_in )
    for m = 1 : size( vol_in, 3) 
        vol_out(:,:,m)				=  conv2( h_row, 1, vol_in(:,:,m), 'same') ;
    end
end


% --------------------------------------------------------------------------------
% see if the current voxel is close to the focus and if it is
% within a functional parcellation area
function yes_no = randomlyKeepVoxel(i,j,k,fc_rad_vox, test_flg)
    yes_no					= 0 ;
% $$$     if rand < 0.2
% $$$ 	% throw away 20% of the voxels
% $$$ 	return
% $$$     end
    if (i^2 + j^2 + k^2) < rand_nums(ri) *(2*fc_rad_vox)^2
        mod( ri + 1, 9999) + 1 ; 
        if test_flg == 0 
            % the test failed - keep with 20% probability
            if rand_nums(ri) < 0.35  ; 
		yes_no				= 1 ;  
            end
            mod( ri + 1, 9999) + 1 ; 
        else
            % the test succeeded - always keep
            yes_no				= 1 ;
        end
    end
    mod( ri + 1, 9999) + 1 ; 
end

% --------------------------------------------------------------------------------
function pv = buildPattern( i,j,k, fc_rad_vox, fc_angles_azim, fc_angles_elev, test_flg )
    
    global ri rand_nums ang_tol rot_mat
    
   pv						= 0 ;

   if test_flg == 0  & rand_nums(ri) < 0.5
       % the test failed - keep with 50% probability
       ri					= mod( ri + 1, 9999) + 1 ;
       return
   end
   
   ri						= mod( ri , 9995) + 1 ;
   
   % pv is a mixture of guassians
   pv						= 0 ;
   for rm = 1 : 4
       x					= [i,j,k]*rot_mat(:,:,rm) ;
       ri					= mod( ri + 1, 9999) + 1 ;
       if ((x(1).^2)/((6*fc_rad_vox)^2) ...
	   +(x(2).^2)/((fc_rad_vox)^2) ... 
	   +(x(3).^2)/((fc_rad_vox)^2)) > rand_nums(ri)/2
           % distance is too large - drop the voxel with increasing prob
           continue ;
       end
       
       pv					=  pv+ rand_nums(ri)*exp( -(x(1).^2)/((6*fc_rad_vox)^2) ...
							 -(x(2).^2)/((fc_rad_vox)^2) ... 
							 -(x(3).^2)/((fc_rad_vox)^2) ) ;
       ri					= mod( ri + 1, 9999) + 1 ;
   end
   
% $$$    pv = pv * rand_nums(ri) ;
% $$$    ri = mod( ri + 1, 9999) + 1 ;
% $$$    
% $$$    
% $$$    pv =  exp( -(i.^2)/2/((fc_rad_vox+rand_nums(ri))^2) ...
% $$$ 	      -(j.^2)/2/((fc_rad_vox+rand_nums(ri+1))^2) ... 
% $$$ 	      -(k.^2)/2/((fc_rad_vox+rand_nums(ri+2))^2) ) ;
% $$$    ri = mod( ri + 3, 9999) + 1 ;
% $$$    
% $$$    azim_ang_dist = min([ abs(fc_angles_azim - atan2(j,i)) ;
% $$$ 		       abs(abs(fc_angles_azim - atan2(j,i)) - 2*pi)]) ;
% $$$    elev_ang_dist = min(abs(fc_angles_elev - atan2(k^2,i^2+j^2) ) ) ;
% $$$    
% $$$    ri = mod( ri , 9999) + 1 ;
% $$$    if  (azim_ang_dist > ang_tol  || elev_ang_dist > ang_tol) & ...
% $$$ 	   rand_nums(ri) < 0.5
% $$$        pv = 0 ;
% $$$    end
     
   ri						= mod( ri + 1, 9999) +1 ;

   if isnan(pv)
       'stop here'
   end
   
end

% $$$ %%%%%%%%% GARBAGE %%%%
% $$$     
% $$$ %     [x_grid,y_grid,z_grid] = meshgrid( spacing, spacing, spacing  ) ;
% $$$ %     
% $$$ %     for k = 1 : size( x_grid, 3)
% $$$ %         
% $$$ %         if(fc_ijk(3) + vox_spacing(k)< 0) || (fc_ijk(3)+ vox_spacing(k)>DIMS(3))
% $$$ %             continue ;
% $$$ %         end
% $$$ %         
% $$$ %         act_blk_smooth = fc_value*exp( -(x_grid(:,:,k).^2)/2/(fc_rad_mm^2+rand) ...
% $$$ %                               -(y_grid(:,:,k).^2)/2/(fc_rad_mm^2+rand) ... 
% $$$ %                               -(z_grid(:,:,k).^2)/2/(fc_rad_mm^2+rand) ) ...
% $$$ %                   /( 2*pi*fc_rad_mm^2)^1.5 ;
% $$$ %                
% $$$ %         [a_l, a_m] = size(act_blk) ;
% $$$ %         
% $$$ %         
% $$$ %         
% $$$ %         % insert this block into the volume
% $$$ %         
% $$$ %     end
