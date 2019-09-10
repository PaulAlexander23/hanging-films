function [yMatrix, t] = odeMatrixSolver(odefunMatrix, t, y0Matrix, timestepper)
    shape = size(y0Matrix);
    
    y0Vector = reshapeToVector(y0Matrix, shape);
    odefunVector = @(t, yVector) vectoriseFunction(odefunMatrix, t, yVector, shape);
    
    [t, yVector] = callTimestepper(odefunVector, t, y0Vector, timestepper);
    
    yMatrix = reshapeToShape(yVector, shape);
end

function fVector = vectoriseFunction(odefunMatrix, t, yVector, shape)
    yMatrix = reshapeToShape(yVector, shape);
    
    fMatrix = odefunMatrix(t, yMatrix);
    
    fVector = reshapeToVector(fMatrix, shape);
end

function Y = reshapeToShape(y, shape)
    Y = squeeze(reshape(y, [shape, numel(y) / prod(shape)]));
end

function y = reshapeToVector(Y, shape)
    y = squeeze(reshape(Y, prod(shape), []));
end

function [t, yVector] = callTimestepper(odefunVector, t, y0Vector, timestepper)
    [t, yVector] = timestepper(odefunVector, t, y0Vector);
    yVector = yVector.';
end