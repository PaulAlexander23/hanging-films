function L = fwibl1(x, Y, params, method)

    delta = params(1);
    theta = params(2);
    Re = params(3);
    C = params(4);

    y = Y(1:end/2,:,:);
    F1 = Y(end/2+1:end,:,:);
    
    deg = [1, 0; 0, 1; 2, 0; 0, 2]';
    dy = method(x, y, deg);
    dF1 = method(x, F1, [1, 0]');
    P = y * cot(theta) - 1/2 * 1 / C * (dy{3} + dy{4});

    dP = method(x, P, [1, 0; 0, 1]');
    
    F2 = - 2/3 * delta * y.^3 .* dP{2};
    
    Ly = -cell2mat(method(x, F1, [1, 0]')) - ...
        cell2mat(method(x, F2, [0, 1]'));
    LF1 = (9 * dy{1} .* F1.^2 - 17 * F1 .* dF1{1}) ./ (7 * y) + ...
        (10 * y - 15 * F1 ./ y.^2 - y .* dP{1}) / (6 * Re);
   
    L = cat(1, Ly, LF1);
end
