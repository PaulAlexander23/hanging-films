function [y, t] = pdeSolver(pdefun, t, domain, y0, timestepper)
    shape = size(y0);
    y0 = reshape(y0, [prod(shape), 1]);
    
    [t, y] = timestepper(@(t, y) callPDEFunctionWithReshape(pdefun, ...
        t, domain, y, shape), t, y0);
    y = y';
    y = squeeze(reshape(y, [shape, length(t)]));
end

function F = callPDEFunctionWithReshape(pdefun, t, domain, y, shape)
    f = pdefun(t, domain, reshapeInput(y, shape));
    
    F = reshapeOutput(y, f, shape);
end

function out = reshapeInput(y, shape)
    out = reshape(y, [shape, numel(y) / prod(shape)]);
end

function F = reshapeOutput(y, f, shape)
    F = reshape(f, [prod(shape), numel(y) / prod(shape)]);
end