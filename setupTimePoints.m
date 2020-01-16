function t = setupTimePoints(arguments)
    
    t = 0:arguments.tStep:arguments.tFinal;
    
    if t(end) ~= arguments.tFinal
        t = [t, arguments.tFinal];
    end
end