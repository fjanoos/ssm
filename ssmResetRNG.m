function ssmResetRNG(env)
%_______________
% ssmResetRNG.m
%
% reset the rng

% $Id: ssmResetRNG.m v0.01 2012-06-23 15:12:05 fj $
    
switch env.CURR_PLATFORM
    case  'bhramand-win64'
        stream = RandStream.getGlobalStream;
        reset(stream);
        rng('default'); %reset random number generator
    case 'numbra'
        rng('default')
end


   
    
   



end