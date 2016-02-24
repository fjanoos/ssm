function bal_factor = ssmComputeStateSpaceBalance( set)
%_______________________________
% ssmComputeStateSpaceBalance.m
%

% $Id: ssmComputeStateSpaceBalance.m v0.01 2012-03-27 13:24:00 fj $
% $Id: ssmComputeStateSpaceBalance.m v0.02 2012-05-29 10:20:11 fj $
    
%%%% propaganda
myLogo						= cafe_logo( mfilename) ;
bal_factor.logo					= myLogo.tmp ;

K						= double( ceil( set.p2 + set.p1*randn(1) )) ;
bal_factor.K_range				= [K : 1 : K+1] ;
bal_factor.ll_range				= [0] ;

% $$$ end




% 2012-06-22	function bal_factor = ssmComputeStateSpaceBalance( set)
% 2012-06-22	%_______________________________
% 2012-06-22	% ssmComputeStateSpaceBalance.m
% 2012-06-22	%
% 2012-06-22	
% 2012-06-22	% $Id: ssmComputeStateSpaceBalance.m v0.01 2012-03-27 13:24:00 fj $
% 2012-06-22	% $Id: ssmComputeStateSpaceBalance.m v0.02 2012-05-29 10:20:11 fj $
% 2012-06-22	% $Id: ssmComputeStateSpaceBalance.m v0.03 2012-06-17 17:09:23 fj $ temporary fix together with CholeskyInverse.m
% 2012-06-22	    
% 2012-06-22	%%%% propaganda
% 2012-06-22	eval([ 'env.logo.' mfilename ' = cafe_logo( mfilename) ;' ])
% 2012-06-22	    
% 2012-06-22	K						= double( ceil( set.p2 + set.p1*randn(1) )) ;
% 2012-06-22	
% 2012-06-22	%%%% out 2012-06-16
% 2012-06-22	% $$$ bal_factor.K_range				= [K : 1 : K+1] ;
% 2012-06-22	% $$$ bal_factor.ll_range				= [0] ;
% 2012-06-22	
% 2012-06-22	%%%% v0.03 2012-06-17 : temporary bug fix together with CholeskyInverse.m fix ...
% 2012-06-22	bal_factor.K_range				= [K - set.p1 : 1 : K+set.p1];
% 2012-06-22	bal_factor.ll_range				= [-1:1];
% 2012-06-22	
% 2012-06-22	
% 2012-06-22	% $$$ end
