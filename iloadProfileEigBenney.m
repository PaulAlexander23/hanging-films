function y0 = iloadProfileEigBenney(x, filename, perturabtionAmplitude, params)
    data = load(filename,'x','y');

    y2 = data.y.';
    h = kron(y2, ones(length(x{1}),1));

    y0 = periodicInterp2(x{1}, data.x{2} * x{2}(end)/data.x{2}(end), h, x{1}, x{2}, 'spline');

    domain = FDDomain({x{2}}, [1,2,3,4], 4, "central");
    eigF = real(eigenfunctionBenney(domain, y0(1,:).', params, 2*pi/x{1}(end)));
    eigF = eigF.' ./ max(abs(eigF));
    perturbation = - perturabtionAmplitude * kron(eigF, cos(2*pi/x{1}(end) * x{1}));

    y0 = y0 + perturbation;
end
