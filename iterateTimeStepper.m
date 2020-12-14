function [y, t, timeTaken] = iterateTimeStepper(odeFunction, t, y0, timeStepper)
    
    tic;
    [t, y] = timeStepper(odeFunction, t, y0);
    timeTaken = toc;

    % Transpose back to column vector because of matlab's default behaviour.
    y = y.';
end
