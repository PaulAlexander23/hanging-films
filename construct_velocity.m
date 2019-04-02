function [z,u,v,w] = construct_velocity(x,y,params,method)
    delta = params(1);
    theta = params(2);
    Re = params(3);
    C = params(5);
    
    zN = 100;
    z = zeros(1,1,zN);
    z(1,1,:) = linspace(0,1,zN);
    
    z = y .* z;

    dy = method(x,y,[1,0]');
    
    p = construct_pressure(x,y,params,method);
    
    u0 = 2 * y .* z - z.^2;
    u1 = (1/2*z.^2 - z .* y) .* cell2mat(...
        method(x,p,[1,0]')) + Re * (1/6 * z.^4 .* y - ...
        2/3 * z.^3 .* y.^2 + 4/3 * z .* y.^4).*dy{1};
    u = u0 + delta*u1;
    
    v1 = (1/2*z.^2 - z .* y) .* cell2mat(method(x,p,[0,1]'));
    v = delta * v1;
    
    w0 = -z.^2 .* dy{1};
    w1 = -cumsum(cell2mat(method(x,u1,[1,0]')) + ...
        cell2mat(method(x,v1,[0,1]')),3) * (z(2)-z(1));
    w = w0 + delta * w1;
end