function F = fdiffusion(x,y,params,method)
    deg = [2,0;0,2]';
    D = params(1);
    dy = feval(method,x,y,deg);
    F = - D*(dy(:,:,1) + dy(:,:,2));
end