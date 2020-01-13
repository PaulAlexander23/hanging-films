function t = setupTimePoints(tFinal, tStep, debug)
    if nargin < 3
        debug = false;
    end
    
    if ~debug
        t = 0:tStep:tFinal;
        if rem(tFinal, tStep) ~= 0
            t = [t, tFinal];
        end
    else
        t = [0,tFinal];
    end
    
end