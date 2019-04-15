function F = fbenney1d(x,y,params,method)
    deg = [1,2];
    delta = params(1);
    theta = params(2);
    Re = params(3);
    We = params(4);
    C = params(5);
    
    dy = feval(method,x,y,deg);
    
    P = y * cot(theta) - 1/2 * 1/C * (dy(:,2));
    
    Q = (2/3 * y.^3 + 8/15 * Re * delta * y .^ 6 .* dy(:,1))...
        - 2/3 * delta * y.^3 .* feval(method,x,P,1);
    
    F = feval(method,x,Q,1);
    
end