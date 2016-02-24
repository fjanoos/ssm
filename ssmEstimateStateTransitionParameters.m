function [stp_estimate] = ssmEstimateStateTransitionParameters( data, params, variational)
%________________________________________
% ssmEstimateStateTransitionParameters.m
%
% function to estimate the state transition parameters using bound maximization
% TODO: 
%  Because of the normalization condition the  weight  vector  for  one  of  
% the  classes  need  not  be estimated. Change the code to avoid estimation
% of the K-the class
    
% $Id: ssmEstimateStateTransitionParameters.m v0.01 2012-06-23 14:27:57 fj $
        
    %%%% propaganda
    myLogo					= cafe_logo( mfilename, 'messg', [ 'subject : ' data.subject_id ' * ' datestr( now, 31)]) ;
    stp_estimate.logo				= myLogo.tmp ;

% $$$     display( 'Estimating State Transition Parameters ...');
    fprintf( '\nEstimating State Transition Parameters ...');

    T						= data.T ;
    s						= data.s ;
    u_idx					= data.u_idx ;
    u						= params.u ;
    CH						= data.ch ;
    K						= data.K ;
    
    % accumulate s_t*s_t'
    s_stranspose				= zeros( CH, CH) ;
    for t = 1 : T
        s_stranspose				= s_stranspose + s( t,:)'*s( t,:) ;     
    end
    
    % this permuting is to allow the (:) operator to list the 3d array with
    % ordering ( row-idx = CH; column-idx = j; slice-idx = i )
    clear temp ;
    temp					= permute( params.W, [ 1,3,2] ) ;
    stp_estimate.W				= params.W ;
    stp_estimate.omega				= params.omega ;
    stp_estimate.w				= [ temp(:) ; params.omega(:)] ;	% in vector form
    %stp_estimate.u				= params.u ;				% use s in this phase for now
    stp_estimate.s				= data.s ;				% use s in this phase for now
    
    %%progress_line( 'init') ;
    
    % use the method http://www.mathworks.com/help/toolbox/optim/ug/brhkghv-7.html
    % to pass additional parameters
    opt_func__ = @(w) optFunc( w, data, variational) ;
    % -- Use optimset/fminuc with Matlab optimization toolbox. --
    %opt_options = optimset('GradObj','on', ...'DerivativeCheck','on',...
    %                         'FunValCheck', 'on', 'Hessian','on',...
    %                          'PlotFcns', @optimplotfunccount,...
    %                         'Diagnostics','on', ...
    %                         'Display', 'notify-detailed') ;
    %[w,fval, exitval, opt_output, grd] = fminunc(opt_func__, stp_estimate.w, opt_options) ;
    % else use this ...
    [ w, opt_output, fval ]			= steep( stp_estimate.w, opt_func__, data.tol, 25);

    %%progress_line( 'term') ;
    
% $$$     display( ['Finished W parameter optimization with details']) ;
    fprintf( '\nFinished W parameter optimization with details\n') ;

    %%%% if verbose...
% $$$     display( opt_output) ;	%%%% column length = 100, see loop below in computeQ
    
    rel_change_w				= norm( w - stp_estimate.w) / norm( stp_estimate.w) ;
    stp_estimate				= updateSTPEstimate( w , K, CH) ;
    stp_estimate.rel_change			= rel_change_w ;
% $$$     display( 'Estimating State Transition Parameters ... done') ;
    fprintf( '\nEstimating State Transition Parameters ... DONE\n') ;
    
    %%%% make sure we exist...
    stp_estimate.logo				= myLogo.tmp ;


end %ssmEstimateStateTransitionParameters

          
%%%% function computeQ -----------------------------------------------------------------------------
function [ grads] = computeQ( data, stp_estimate, variational, hessian_flag)
%@params 
%   hessian_flag : if true, compute hessian
 
    K						= data.K ;
    T						= data.T ;
    s						= data.s ;
    CH						= data.ch ;
    W						= stp_estimate.W ;
    omega					= stp_estimate.omega ;
    q_x						= variational.q_x ;
    lambda_w					= data.lambda_w ;
    
    grads.Q					= 0 ;					% the (+ve) log-likelihood
    grads.Qgrad_W				= zeros( size( W) ) ;			% gradient wrt w
    grads.Qgrad_omega				= zeros( size( omega)) ;		% wrt omega
    hessian					= sparse( K * ( K+1) * CH, K * ( K+1) * CH) ;	% hessian matrix
    
    % the B matrix needed for bound optimization - precompute
    clear initialize_B
    persistent B ;    
    if isempty( B)
        %initialize_B = true ;
        %B =  sparse(K*(K+1),K*(K+1)) ;  % hessian matrix ;
        %display('initializing B matrix') ;
    end
    
% $$$     %h = waitbar(0, 'computing Q', 'Name', 'ssmEstimateStateTransitionParameters') ;

    %%%% PISTI : 100 corresponds with opt_output
    for t = 2 : mod( T, 100): T

% $$$         %waitbar( t/T, h) ;
	%%progress_line( 2, T, t) ;
	
% returns pm(i,j) = p(x_{t}=j|x_{t-1}=i)
        pm					= probStateTransitionMatrix( t, stp_estimate, data) ;        
        for i = 1 : K
            for j = 1 : K           
                grads.Q				= grads.Q  + q_x( t-1,i)*q_x( t,j)*log( pm( i,j) +eps ) ;
                
                grads.Qgrad_W( :,i,j)		= grads.Qgrad_W( :,i,j)...
		    + q_x( t-1,i)*( q_x( t,j)-pm( i,j))*s( t,:)' ;
                grads.Qgrad_omega( :,j)		=  grads.Qgrad_omega( :,j) ...
		    + q_x( t-1,i)*( q_x( t,j)-pm( i,j))*s( t,:)' ;
		
                            
                ij_idx				= [(( i-1)*K+( j-1))*CH+1:(( i-1)*K+( j-1))*CH+CH] ;
                j_idx				= [( K*K+j-1)*CH+1:( K*K+j-1)*CH+CH] ;
                %set up the bound-optimization B matrix once
                if exist( 'initialize_B', 'var')
                    B( ij_idx , ij_idx)		= B( ij_idx , ij_idx)+0.5 ;   
                    B( ij_idx, j_idx)		= B( ij_idx, j_idx) + 0.5 ;  
                    B( j_idx, ij_idx)		= B( ij_idx, j_idx) ;
                    B( j_idx, j_idx)		= B( j_idx, j_idx) + 0.5 ;  
                end
                if hessian_flag
                    % initialize the diagonal \Hess_w_{i,j}_{i,j}
                    hessian( ij_idx , ij_idx)	= hessian( ij_idx , ij_idx) + q_x( t-1,i)*pm( i,j)*s( t,:)'*s( t,:) ;
                    %\Hess_w_{i,j}_omega_{j}
                    hessian( ij_idx, j_idx)	= hessian( ij_idx, j_idx) + q_x( t-1,i)*pm( i,j)*s( t,:)'*s( t,:) ;  
                    hessian( j_idx, ij_idx)	= hessian( ij_idx, j_idx) ;
                    %\Hess_omega_{j}_{j}       
                    hessian( j_idx, j_idx)	= hessian( j_idx,j_idx) + q_x( t-1,i)*pm( i,j)*s( t,:)'*s( t,:) ;                                                                 
                end %if hessian_flag
                
                for jp = 1 : K                            
                    ijp_idx			= [(( i-1)*K+jp-1)*CH+1:(( i-1)*K+jp-1)*CH+CH] ;
                    jp_idx			= [( K*K+jp-1)*CH+1:( K*K+jp-1)*CH+CH] ;
                    if hessian_flag
                        hessian( ij_idx , ijp_idx)	= hessian( ij_idx , ijp_idx) - q_x( t-1,i)*pm( i,j)*pm( i,jp)*s( t,:)'*s( t,:) ;
                        hessian( ij_idx , jp_idx)	= hessian( ij_idx , jp_idx) - q_x( t-1,i)*pm( i,j)*pm( i,jp)*s( t,:)'*s( t,:) ;
                        hessian( jp_idx , ij_idx)	= hessian( ij_idx , jp_idx) ;
                        hessian( j_idx , jp_idx)	= hessian( j_idx , jp_idx) - q_x( t-1,i)*pm( i,j)*pm( i,jp)*s( t,:)'*s( t,:) ;
                    end
                    if exist( 'initialize_B', 'var')
                        B( ij_idx , ijp_idx)	= B( ij_idx , ijp_idx)-0.5/K ;   
                        B( ij_idx, jp_idx)	= B( ij_idx, jp_idx) -0.5/K ;  
                        B( jp_idx, ij_idx)	= B( ij_idx, jp_idx) ;
                        B( j_idx, jp_idx)	= B( j_idx, jp_idx) -0.5/K ; 
                    end
                end %jp     
	    end   %for j
        end  %for i              
    end
    %delete( h) ;
    
    % the prior on W    
    for j = 1 : K
        for i = 1 : K        
            w_ij				= W( :,i,j) ;
            grads.Q				= grads.Q - 0.5*lambda_w*( w_ij'*w_ij) ;
            grads.Qgrad_W( :,i,j)		= grads.Qgrad_W( :,i,j)- lambda_w*W( :,i,j) ;
        end
    end
                           
    % this permuting is to allow the (:) operator to list the 3d array with
    % ordering ( row-idx = CH; column-idx = j; slice-idx = i )
    temp					= permute( grads.Qgrad_W , [ 1,3,2] ) ;
    grads.Qgrad_w				= [ temp( :) ; grads.Qgrad_omega( :) ] ;
    
    if hessian_flag
        grads.Qhessian				= -hessian - lambda_w * [ speye( K^2*CH) , sparse( K^2*CH,K*CH)  ; sparse( K*CH,K*CH), sparse( K*CH,K^2*CH) ] ;
        if exist( 'initialize_B', 'var') 
            B					= -kron( B,s_stranspose)-lambda_w * [ speye( K^2*CH) , sparse( K^2*CH,K*CH) ; sparse( K*CH,K*CH), sparse( K*CH,K^2*CH) ] ;
        end
        %[ev] = eig( grads.Qhessian-B) ;
        grads.B					= B ;
    end
end %function computeGradientQ


%%%% reshapes w into W and omega -------------------------------------------------------------------
function [ stp_estimate] = updateSTPEstimate( w , K, CH)
    
    temp					= reshape( w( 1:CH*K^2), [ CH, K, K] ) ;
    stp_estimate.W				= permute( temp, [ 1,3,2] ) ;
    stp_estimate.omega				= reshape( w( CH*K^2+1:CH*K^2+CH*K), [ CH, K]) ;
    stp_estimate.w				= w ;
    
end %updateSTPEstimate


%%%% function to be fed into the optimizer fminunc -------------------------------------------------
function [ f,g,H]= optFunc( w, data, variational)

    [ stp_estimate]				= updateSTPEstimate( w , data.K, data.ch) ;
    if nargout > 2
         hessian_flag				= true ;
    else
        hessian_flag				= false ;
    end
    [ grads]					= computeQ( data, stp_estimate, variational, hessian_flag) ;
    f						= -grads.Q ;
    
% $$$     display( ['Q = ',num2str(f)] );
% $$$     fprintf( ['Q = %s\n'], num2str(f)) ;
% $$$     fprintf( ['%s '], num2str(f)) ;
    fprintf( ['  Q = %s\n'], num2str( f)) ;
    
    if nargout > 1
        g					= -grads.Qgrad_w ;
    end
    if nargout > 2
        H					= -grads.Qhessian ;
    end
end
