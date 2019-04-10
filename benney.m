function [F, J] = benney(x, y, params, method, getD)

    delta = params(1);
    theta = params(2);
    Re = params(3);
    %We = params(4);
    C = params(5);

    %L = cellfun(@(x) x(end),x);
    e1 = zeros(1, 1, 2);
    e1(1, 1, 1) = 1;

    deg = [1, 0; 0, 1; 2, 0; 0, 2]';
    dy = method(x, y, deg);

    P = y * cot(theta) - 1/2 * 1 / C * (dy{3} + dy{4});
    %P = H * cot(theta) - We * R(H, L) - 1/2 * 1/C * lap(H, L);

    gradP = method(x, P, [1, 0; 0, 1]');
    gradP = cat(3, gradP{1}, gradP{2});

    Q = (2/3 * y.^3 + 8/15 * Re * delta * y.^6 .* dy{1}) .* e1 ...
        -2/3 * delta * y.^3 .* gradP;

    F = -cell2mat(method(x, Q(:, :, 1), [1, 0]')) - ...
        cell2mat(method(x, Q(:, :, 2), [0, 1]'));

end
