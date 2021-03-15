function y0 = iloadProfileWIBL1STF(x, filename, perturabtionAmplitude)
    load(filename,'y');

    y2 = reshape(y, length(y)/2, 2).';

    y0 = kron(y2, ones(length(x{1}),1));

    perturbation = - perturabtionAmplitude * cos(2*pi/x{1}(end) * x{1});

    y0(1:end/2,:) = y0(1:end/2,:) + perturbation;
end
