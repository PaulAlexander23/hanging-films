function t = setupTimePoints(arguments, debug)
    if nargin < 3
        debug = false;
    end
    
    if ~debug
        t = 0:arguments.tStep:arguments.tFinal;
        if t(end) ~= arguments.tFinal
            t = [t, arguments.tFinal];
        end
    else
        t = [0, arguments.tFinal];
    end
    
end