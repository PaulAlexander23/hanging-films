function F = fbenney1d(x,y,params,method)
    deg = [1,2];
    delta = params(1);
    theta = params(2);
    Re = params(3);
    C = params(4);
    
    dy = feval(method,x,y,deg);
    
    P = y * cot(theta) - 1/2 * 1/C * (dy{2});
    
    Q = (2/3 * y.^3 + 8/15 * Re * delta * y .^ 6 .* dy{1})...
        - 2/3 * delta * y.^3 .* cell2mat(method(x,P,1));
    
    F = -cell2mat(method(x,Q,1));
    
end