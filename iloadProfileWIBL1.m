function y0 = iloadProfileWIBL1(x, filename)
    load(filename,'y');

    y2 = reshape(y, length(y)/3, 3).';

    y0 = kron(y2, ones(length(x{1}),1));
end
