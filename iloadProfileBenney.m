function y0 = iloadProfileBenney(x, filename)
    load(filename,'y');

    y0 = repmat(y.', [length(x{1}),1]);
end
