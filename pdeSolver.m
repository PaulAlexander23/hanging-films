function [y, t] = pdeSolver(pdefun, t, x, y0, diffMethod, timestepper)
    
    shape = size(y0);
    y0 = reshape(y0, [prod(shape), 1]);
    
    [t, y] = timestepper(@(t, y) callPDEFunctionWithReshape(pdefun, ...
        t, x, y, diffMethod, shape), t, y0);
    y = y';
    y = squeeze(reshape(y, [shape, length(t)]));
end

function F = callPDEFunctionWithReshape(pdefun, t, x, y, diffMethod, shape)
    f = pdefun(t, x, reshapeInput(y, shape), diffMethod);
    F = reshapeOutput(y, f, shape);
end

function out = reshapeInput(y, shape)
    out = reshape(y, [shape, numel(y) / prod(shape)]);
end

function F = reshapeOutput(y, f, shape)
    F = reshape(f, [prod(shape), numel(y) / prod(shape)]);
end