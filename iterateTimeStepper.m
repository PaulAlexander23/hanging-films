function [y, t, timeTaken] = iterateTimeStepper(odeFunction, t, y0, timeStepper)
    
    tic
    [y, t] = odeMatrixSolver(odeFunction, t, y0, timeStepper);
    timeTaken = toc;
    
end