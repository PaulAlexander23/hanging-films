function F = fbenney1d(domain,y,params)
    P = y * cot(params.theta) - 1/2 * 1/params.C * domain.diff(y,2);
    
    Q = (2/3 * y.^3 + 8/15 * params.Re * y .^ 6 .* domain.diff(y, 1))...
        - 2/3 * y.^3 .* domain.diff(P,1);
    
    F = -domain.diff(Q, 1);
end
