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

    if nargout > 1
        deg2 = [1, 0; 0, 1; 2, 0; 0, 2]';
        D = getD(deg2);

        Lap = D{3} + D{4};
        Lap_Dx = D{1} * Lap;
        Lap_Dy = D{2} * Lap;

        n = numel(y);

        v = reshape(y, n, 1);
        dv = cellfun(@(z) reshape(z, n, 1), dy, ...
            'UniformOutput', false);

        q1_h = spdiags(2 * v.^2 .* (1 - dv{1} * cot(theta) + ...
            (Lap_Dx * v) / (2 * C)), 0, n, n) + ...
            2 * v.^3/3 .* (-cot(theta) * D{1} + Lap_Dx / (2 * C));

        q1_h = q1_h + ...
            Re * 8 * v.^6/15 .* D{1} + ...
            spdiags(Re * 16 * v.^5 .* dv{1} / 5, 0, n, n);

        q2_h = spdiags(2 * v.^2 .* (-dv{2} * cot(theta) + ...
            (Lap_Dy * v) / (2 * C)), 0, n, n) + ...
            2 * v.^3/3 .* (-cot(theta) * D{2} + Lap_Dy / (2 * C));

        D = getD([1, 0; 0, 1]');

        J = -(D{1} * q1_h + D{2} * q2_h);
    end

end
