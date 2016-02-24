function MI = ssmComputeMI( data_a, params_a, data_b, params_b)
%________________
% ssmComputeMI.m
%
% compute the symmetric MI between two subjects a and b using the
% given data.
% 
% Run this file after you have estimated the models for multiple
% subjects.  also - you have write specific small scripts to query /
% analyze results.

% $Id: ssmComputeMI.m v0.01 2012-05-29 15:10:59 fj $

%%%% propaganda
myLogo						= cafe_logo( mfilename) ;


K_a						= params_a.K ;
K_b						= params_b.K ;
    
sseq_a						= ssmEstimateOptimalSSeq( params_a, data_a) ;
sseq_b						= ssmEstimateOptimalSSeq( params_b, data_a) ;
x_a						= sseq_a.x ;
x_b						= sseq_b.x ;
if length( x_a) ~= length( x_b)
    error ( 'ssmComputeMI: Not same length sequences') ;
end

MI_1						= mutualinfo( x_a, x_b) ;

sseq_a						= ssmEstimateOptimalSSeq( params_a, data_b) ;
sseq_b						= ssmEstimateOptimalSSeq( params_b, data_b) ;
x_a						= sseq_a.x ;
x_b						= sseq_b.x ;
if length( x_a) ~= length( x_b)
    error ( 'ssmComputeMI: Not same length sequences') ;
end

MI_2						= mutualinfo( x_a, x_b) ;

MI						= 0.5*( MI_1 + MI_2) ;

% $$$ end
