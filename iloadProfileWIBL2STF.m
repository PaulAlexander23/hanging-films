function y0 = iloadProfileWIBL2STF(x, filename, perturabtionAmplitude)
    load(filename,'y');

    y2 = reshape(y, length(y)/3, 3).';

    y0 = kron(y2, ones(length(x{1}),1));

    perturbation = - perturabtionAmplitude * cos(2*pi/x{1}(end) * x{1});

    y0(1:end/3,:) = y0(1:end/3,:) + perturbation;
end
