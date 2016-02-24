function [u_estimate] = ssmEstimateHiddenStimulus( data, params, variational)
%_____________________________
% ssmEstimateHiddenStimulus.m
%
% function to estimate the hidden stimulus 

% $Id: ssmEstimateHiddenStimulus.m v0.01 2012-06-23 14:30:27 fj $
% $Id: ssmEstimateHiddenStimulus.m v0.02 2012-07-05 12:38:50 fj $ using progress_line
    
    %%%% propaganda
    myLogo					= cafe_logo( mfilename, 'messg', [ 'subject : ' data.subject_id ' * ' datestr( now, 31)]) ;
    u_estimate.logo				= myLogo.tmp ;

% $$$     display( ['Estimating Hidden stimulus '] ) ;
    fprintf( [ '\nEstimating Hidden stimulus ...'] ) ;
    
    K						= data.K ;
    T						= data.T ;
    s						= data.s ;
    CH						= data.ch ;
    W						= params.W ;
    omega					= params.omega ;
    q_x						= variational.q_x ;
    lambda_w					= data.lambda_w ;
    u_idx					= data.u_idx ;
    u_len					= length( u_idx) ;
    
    u_estimate.W				= params.W ;
    u_estimate.omega				= params.omega ;
    u_estimate.u				= params.u ;
    u_estimate.err.ch				= zeros( T, CH) ;
    u_estimate.err.ch_tot			= zeros( 1, CH) ; ;
    u_estimate.err.all				= zeros( T, 1) ;
    u_estimate.err.all_tot			= 0 ;

    %%%% PISTI : 2012-04-03
% $$$     h = waitbar(0, 'Estimating u', 'Name', 'ssmEstimateHiddenStimulus') ;
    %%progress_line( 'init') ;
    
    for idx = 1 : u_len
	%%%% PISTI : 2012-04-03
% $$$         waitbar( idx/u_len, h) ;
	%%progress_line( 1, u_len, idx) ;
	
        t					= u_idx( idx) ;
        % initialize the unobserved stimulus
        opt_func__				= @(u) optFunc( u, t, data, u_estimate, variational) ;
% Use optimset/fminuc with Matlab optimization toolbox. --
%         opt_options = optimset('GradObj','on', ...% 'DerivativeCheck','on',...
%                         'FunValCheck', 'on', 'Hessian','on',...                       
%                         ...% 'Diagnostics','on', ...
%                         'Display', 'notify-detailed' ) ;
%         [u_t,fval, exitval, opt_output, grd] = ...
%                     fminunc(opt_func__, u_estimate.u(t,:), opt_options) ;
% else use this ...
        [ u_t, opt_output, fval ]		= cgtrust( u_estimate.u( t,:)', opt_func__, data.tol) ;        
        % the hidden stimulus
        u_estimate.u( t,:)			= u_t' ;
        % error per stimulus channel at t
        u_estimate.err.ch( t,:)			= abs( u_estimate.u( t,:) - s( t,:))./ abs( s( t,:) + ( s( t,:)==0) ) ;
        % total error per stimulus channel
        u_estimate.err.ch_tot			= u_estimate.err.ch_tot + u_estimate.err.ch( t,:) ;
        % total error at time point        
        u_estimate.err.all( t)			= sum( u_estimate.err.ch( t,:)) ;
        % total error at all 
        u_estimate.err.all_tot			= u_estimate.err.all_tot + u_estimate.err.all( t) ;
    end
    %%%% PISTI : 2012-04-03
% $$$     close( h)
    %%progress_line( 'term') ;
    
    u_estimate.err.all_tot			= u_estimate.err.all_tot / u_len ;
    u_estimate.err.ch_tot			= u_estimate.err.ch_tot / u_len ;
    
% $$$     display( ['Estimating Hidden stimulus ... DONE'] ) ;
% $$$     display( ['Error per channel = ', num2str(u_estimate.err.ch_tot) ] ) ;
% $$$     display( ['Total Error = ', num2str(u_estimate.err.all_tot)] ) ;            
    fprintf( [ '\nEstimating Hidden stimulus ... DONE'])  ;
    fprintf( [ '\nError per channel = %s'], filter_string( num2str( u_estimate.err.ch_tot), '  ', ' ')) ;
    fprintf( [ '\nTotal Error = %s\n'], num2str( u_estimate.err.all_tot))  ;
end


%%%% compute the log-likelihood of u_t, gradient and hessian ---------------------------------------
function [ Q_t] = computeQ( t, data, u_estimate, variational, hessian_flag)

    K						= data.K ;
    T						= data.T ;
    s						= data.s ;
    CH						= data.ch ;
    W						= u_estimate.W ;
    omega					= u_estimate.omega ;
    u_t						= u_estimate.u( t,:) ;
    q_x						= variational.q_x ;
    
    % output values
    Q_t.fnc_val					= 0 ;			% the free energy
    Q_t.grad					= zeros( CH, 1) ;	% gradient wrt w
    Q_t.hessian					= zeros( CH, CH) ;	% hessian matrix
        
    pm						= probStateTransitionMatrix( t, u_estimate, data)  ;
    for i = 1 : K
	for j = 1 : K
	    W_ij				= W( :,i,j) + omega( :,j)  ;
	    Q_t.fnc_val				= Q_t.fnc_val + q_x( t - 1,i) * q_x( t,j) * log( pm( i,j) + eps ) ;
	    Q_t.grad				= Q_t.grad + q_x( t - 1,i) * q_x( t,j) * W_ij ;
	    
	    for k = 1 : K
		W_ik				= W( :,i,k) + omega( :,k) ;
		Q_t.grad			= Q_t.grad - q_x( t - 1,i) * q_x( t,j) * pm( i,k) * W_ik ;
		
		if hessian_flag   
		    Q_t.hessian			= Q_t.hessian + q_x( t - 1,i) * q_x( t,j) * pm( i,k) * W_ik * W_ik' ;
		    
		    for l = 1 : K
			W_il			= W( :,i,l) + omega( :,l) ;
			Q_t.hessian		= Q_t.hessian - q_x( t - 1,i) * q_x( t,j) * pm( i,k) * W_ik * pm( i,l) * W_il' ;
		    end %l    
		end
	    end %k
	end %j
    end %i
    
    Q_t.hessian					= -Q_t.hessian ;
end


%%%% function to be fed into the optimizer fminunc -------------------------------------------------
function [ f,g,H]= optFunc( u, t, data, u_estimate, variational)

    u_estimate.u( t,:)				= u ;
    %u_estimate.s( t,:)				= u ;
    if nargout > 2
         hessian_flag				= true ;
    else
        hessian_flag				= false ;
    end    
    
    [ Q_t]					= computeQ( t, data, u_estimate, variational, hessian_flag) ;

    f = -Q_t.fnc_val  ; 
    if nargout > 1
        g					= -Q_t.grad ;
    end
    if nargout > 2
        H					= -Q_t.hessian ;
    end
    %display( [ 'Q = ',num2str( f)] ) ;
    
end
