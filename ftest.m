function F = ftest(x,y,params,method)
    deg = [1,0;2,0;0,2]';
    D = params(1);
    
    e1 = zeros(1,1,2);
    e1(1,1,1) = 1;
    dy = feval(method,x,y,deg);
    
    Q = y .* e1 - D * grad(x,dy(:,:,2)+dy(:,:,3));
    F = div(x,Q);
    
    function out = div(x,y)
        divdeg = cat(3,[1,0;0,0],[0,0;0,1]);
        d2y = diff_ps_2d(x,y,divdeg);
        out = d2y(:,:,1) + d2y(:,:,2);
    end
    
    function out = grad(x,y)
        graddeg = [1,0;0,1]';
        out = diff_ps_2d(x,y,graddeg);
    end
end