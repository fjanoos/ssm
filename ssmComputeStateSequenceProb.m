function prob = ssmComputeStateSequenceProb(data, params, sseq)
%_______________________________
% ssmComputeStateSequenceProb.m
%
% compute the probability of a state sequence using the generative model

% $Id: ssmComputeStateSequenceProb.m v0.01 2012-06-23 15:07:52 fj $

%%%% propaganda
% $$$ myLogo						= cafe_logo( mfilename) ;

x = sseq.x;

% track p(x)
ln_p_x = 0;

% initialize t = 1 with the the invariant density
pm = probStateTransitionMatrix(1, params, data); 
ps = pm^100; % one row gives the invariant density
ln_p_x = log( ps(1,x(1)) ); 

for t = 2 : T
   ln_pm = log( probStateTransitionMatrix(t, params, data) );  
   ln_p_x =  ln_p_x + ln_pm(x(t-1), x(t) ; % state transition probabilities    
end

prob = ln_p_x;

end