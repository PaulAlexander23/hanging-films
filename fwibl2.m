function dYdt = fwibl2(domain, Y, params)
    Y = domain.reshapeToDomain(Y);

    h = Y(1:end/3, :, :);
    F1 = Y(end/3+1:2*end/3, :, :);
    F2 = Y(2*end/3+1:end, :, :);
    
    beta = - beta * pi/180;

    theta = pi/2 + beta;
    g = 9.81;

    ReWave = h_N^3 * g * sin(theta) / (3 * nu^2);
    KaWave = sigma / (rho * (g * sin(theta))^(1/3) * nu^(4/3));
    CtWave = cot(theta);

    delta = (3 * ReWave)^(11/9) / (KaWave^(1/3));
    eta = (3 * ReWave)^(4/9) / (KaWave^(2/3));
    zeta = CtWave * (3 * ReWave)^(2/9) / (KaWave^(1/3));

    dhdx = domain.diff(h, [1; 0]);
    d2hdx2 = domain.diff(h, [2; 0]);
    d3hdx3 = domain.diff(h, [3; 0]);
    dhdy = domain.diff(h, [0; 1]);
    d2hdy2 = domain.diff(h, [0; 2]);
    d3hdy3 = domain.diff(h, [0; 3]);
    d2hdxy = domain.diff(h, [1; 1]);
    d3hdx2y = domain.diff(h, [2; 1]);
    d3hdxy2 = domain.diff(h, [1; 2]);

    dF1dx = domain.diff(F1, [1; 0]);
    d2F1dx2 = domain.diff(F1, [2; 0]);
    dF1dy = domain.diff(F1, [0; 1]);
    d2F1dy2 = domain.diff(F1, [0; 2]);
    d2F1dxy = domain.diff(F1, [1; 1]);
    
    dF2dx = domain.diff(F2, [1; 0]);
    d2F2dx2 = domain.diff(F2, [2; 0]);
    dF2dy = domain.diff(F2, [0; 1]);
    d2F2dy2 = domain.diff(F2, [0; 2]);
    d2F2dxy = domain.diff(F2, [1; 1]);
    



    dhdt = - dF1dx - dF2dy;

    dF1dt = delta * (9/7 * F1.^2 ./ h * dhdx - 17/7 * F1 ./ h * dF1dx) ...
        + (5/6 * h - 5/2 * F1 ./ h.^2 + delta * ( ...
        - 8/7 * F1 ./ h .* dF2dy - 9/7 * F2 ./ h .* dF1dy ...
        + 9/7 * F1 .* F2 ./ h.^2 .* dhdy) ...
        + eta * (4 * F1 ./ h.^2 * dhdx.^2 - 9/2 * dF1dx .* dhdx ./ h ...
        - 6 * F1 ./ h .* d2hdx2 + 9/2 * d2F1dx2 ...
        + 13/4 * F2 ./ h.^2 .* dhdx .* dhdy - dF1dy .* dhdy ./ h ...
        - 43/16 * dF2dx ./ h .* dhdy - 13/16 * dF2dy ./ h .* dhdx ...
        + 3/4 * F1 ./ h.^2 .* dhdy.^2 - 23/16 * F1 .* d2hdy2 ./ h ...
        - 73/16 * F2 .* d2hdxy ./ h + d2F1dy2 + 7/2 * d2F2dxy) ...
        - 5/6 * zeta * h .* dhdx + 5/6 * h .* (d3hdx3 + d3hdxy2)) ...
        ./ (1 - delta/70 * F1 .* dhdx);

    dF2dt = delta * (9/7 * F2.^2 ./ h * dhdy - 17/7 * F2 ./ h * dF2dy) ...
        - 5/2 * F2 ./ h.^2 ...
        + delta * (- 8/7 * F2 ./ h .* dF1dx - 9/7 * F1 ./ h .* dF2dx ...
        + 9/7 * F1 .* F2 ./ h.^2 .* dhdx) ...
        + eta * (4 * F2 ./ h.^2 * dhdy.^2 - 9/2 * dF2dy .* dhdy ./ h ...
        - 6 * F2 ./ h .* d2hdy2 + 9/2 * d2F2dy2 ...
        + 13/4 * F1 ./ h.^2 .* dhdx .* dhdy - dF2dx .* dhdx ./ h ...
        - 43/16 * dF1dx ./ h .* dhdx - 13/16 * dF1dx ./ h .* dhdy ...
        + 3/4 * F2 ./ h.^2 .* dhdx.^2 - 23/16 * F2 .* d2hdx2 ./ h ...
        - 73/16 * F1 .* d2hdxy ./ h + d2F2dx2 + 7/2 * d2F1dxy) ...
        - 5/6 * zeta * h .* dhdy + 5/6 * h .* (d3hdx2y + d3hdy3);
      

    dYdt = cat(1, dhdt, dF1dt, dF2dt);

    dYdt = domain.reshapeToVector(dYdt);
end
