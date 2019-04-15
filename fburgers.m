function F = fburgers(x,y,params,method)
    deg = [0,1;2,0;0,2]';
    D = params(1);
    dy = feval(method,x,y,deg);
    F = 10*dy(:,:,1)... * y.*dy
        - D*(dy(:,:,2) + dy(:,:,3));
end