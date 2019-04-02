function p = construct_pressure(x,y,params,method)
    theta = params(2);
    C = params(5);
    
    dy = method(x,y,[2,0;0,2]');
    
    p = 2*y*cot(theta) - 1/C * (dy{1} + dy{2});
    % p = 2*y*cot(theta) - We * R(y) - 1/C * (dy{1} + dy{2});
end