function [y, timeTaken] = nonlinearSolve(odeFunction, y0, nonlinearSolver)
    
    tic;
    y = nonlinearSolver(odeFunction, y0);
    timeTaken = toc;

    % Transpose back to column vector because of matlab's default behaviour.
    y = y.';
end
