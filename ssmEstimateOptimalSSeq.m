function [sseq] = ssmEstimateOptimalSSeq( params, data, hidden)
%__________________________
% ssmEstimateOptimalSSeq.m
%
% Estimates optimal state sequence given a set of parameters

% $Id: ssmEstimateOptimalSSeq.m v0.01 2012-05-30 19:07:19 fj $

%%%% propaganda
myLogo						= cafe_logo( mfilename, 'messg', [ 'subject : ' data.subject_id ' * ' datestr( now, 31)]) ;
sseq.logo					= myLogo.tmp ;
    
D						= data.D ;
L						= data.L ;
T						= data.T ;
Y						= data.y ;
s						= data.s ;
TR						= data.TR ;
ch						= data.ch ;
K						= data.K ;
lambda_w					= data.lambda_w ;
tol						= data.tol ;

% initialize
sseq						= [] ;
sseq.x						= randi( K, [T,1] ) ;
sseq.log_p_opt					= -10e5 ; %some small number
    
it_cnt						= 0 ;
hf1						= figure ;
if exist( 'hidden', 'var')
    hf2						= figure ;
end
clear variational rel_change rel_err
rel_err						= {} ;

while(1)
    it_cnt 					= it_cnt + 1 ;
    % Variational E-Step
    variational_updated				= ssmVariationalDensityOptimalSSeq(params, data, sseq) ;
% $$$     sseq.logo.( variational_updated.logo.name)	= variational_updated.logo ;	%%%% LOGO huckepack...
    
    if exist( 'variational', 'var')
	rel_change{it_cnt}.Sigma_qz		= norm(variational_updated.Sigma_qz (:) - variational.Sigma_qz (:))...
	    /norm(variational.Sigma_qz (:)) ;
	rel_change{it_cnt}.mu_qz		= norm(variational_updated.mu_qz  (:) - variational.mu_qz  (:))...
	    /norm(variational.mu_qz (:)) ;
	figure( hf1) ;
	subplot( 3,3,1) ; hold all ; plot( it_cnt, rel_change{ it_cnt}.mu_qz) ;        
	title( 'variational rel-change \mu_{q_z}') ;     
	subplot( 3,3,2) ; hold all ; plot( it_cnt, rel_change{ it_cnt}.Sigma_qz) ;        
	title( 'variational rel-change \Sigma_{q_z}') ;     
    end
    variational					= variational_updated ;       
    
    
    % use for ground truth testing
    if exist( 'hidden', 'var')
	rel_err{it_cnt}				= 0 ;
	for t = 1:T
	    rel_err{ it_cnt}			= rel_err{ it_cnt} + (hidden.x( t)~=sseq.x( t)) ;
	end
	rel_err{ it_cnt}			= rel_err{ it_cnt} / T ;
	figure( hf2) ;
	hold all ; plot( it_cnt, rel_err{ it_cnt}, 'k.') ;
	title( 'Rel. error state-sequence estimate') ;
    end
    
    sseq_update					= ssmVariationalViterbi( data, params, variational) ;
% $$$     sseq.logo.( sseq_update.logo.name)		= sseq_update.logo ;	%%%% LOGO huckepack...
    
% $$$     display ([' p_opt = ', num2str(sseq_update.log_p_opt)]) ;
    fprintf([ '\np_opt = %s\n\n'], num2str( sseq_update.log_p_opt)) ;
    
    if sseq_update.log_p_opt >= sseq.log_p_opt 
	warning (' your maximization has failed ') ;
    end
    if abs( sseq_update.log_p_opt - sseq.log_p_opt)/abs( sseq_update.log_p_opt) < tol ;
	break
    end
    
    sseq					= sseq_update ;
    
end %while

% $$$ end
