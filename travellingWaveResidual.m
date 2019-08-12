function [res, c] = travellingWaveResidual(domain, h, params, c)
    
    dhdt = fwibl1(domain, h, params);
    
    dydx = domain.diff(h(1:end/2,:), [1, 0]');
    dFdx = domain.diff(h(end/2+1:end,:), [1, 0]');
    
    mean(dhdt./[dydx; dFdx],'all')
    
    if nargin < 4
        c = fminunc(@(c) norm(dhdt./[dydx; dFdx] + c), 2);
    end
    norm(dhdt + c * [dydx; dFdx])
    res = dhdt + c * [dydx; dFdx];
end