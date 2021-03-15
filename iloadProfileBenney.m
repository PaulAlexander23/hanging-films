function y0 = iloadProfileBenney(x, filename, perturabtionAmplitude)
    load(filename,'y');

    perturbation = - perturabtionAmplitude * cos(2*pi/x{1}(end) * x{1});

    y0 = y.' + perturbation;
end
