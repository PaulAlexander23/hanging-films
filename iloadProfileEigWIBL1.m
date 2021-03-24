function y0 = iloadProfileEigWIBL1(x, filename, perturabtionAmplitude, params)
    data = load(filename,'x','y');

    y2 = reshape(data.y, length(data.y)/3, 3).';
    h = kron(y2, ones(length(x{1}),1));

    h0 = periodicInterp2(x{1}, data.x{2} * x{2}(end)/data.x{2}(end), h(1:end/3,:), x{1}, x{2}, 'spline');
    f10 = periodicInterp2(x{1}, data.x{2} * x{2}(end)/data.x{2}(end), h(1+end/3:2/3*end,:), x{1}, x{2}, 'spline');
    f20 = periodicInterp2(x{1}, data.x{2} * x{2}(end)/data.x{2}(end), h(1+2*end/3:end,:), x{1}, x{2}, 'spline');

    y0 = [h0; f10; f20];

    domain = FDDomain({x{2}}, [1,2,3,4], 4, "central");
    eigF = real(eigenfunctionWIBL1(domain, h0(1,:).', f10(1,:).', f20(1,:).', params, 2*pi/x{1}(end)));
    eigF = eigF ./ max(abs(eigF(1:end/3)));
    eigF2 = reshape(eigF, length(eigF)/3, 3).';
    perturbation = - perturabtionAmplitude * kron(eigF2, cos(2*pi/x{1}(end) * x{1}));

    y0 = y0 + perturbation;
end
