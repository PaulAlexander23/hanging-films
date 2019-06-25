function F = fbenney1d(domain,y,params)
    delta = params(1);
    theta = params(2);
    Re = params(3);
    C = params(4);
    
    P = y * cot(theta) - 1/2 * 1/C * domain.diff(y,2);
    
    Q = (2/3 * y.^3 + 8/15 * Re * delta * y .^ 6 .* domain.diff(y, 1))...
        - 2/3 * delta * y.^3 .* domain.diff(P,1);
    
    F = -domain.diff(Q, 1);
end
