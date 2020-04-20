function [y, t, timeTaken] = iterateTimeStepper(odeFunction, t, y0, timeStepper)
    
    tic
    [t, y] = timeStepper(odeFunction, t, y0);
    timeTaken = toc;

    y = y';
end
