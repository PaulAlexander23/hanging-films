function y0 = iloadProfileWIBL1STF(x, filename)
    load(filename,'y');

    y2 = reshape(y, length(y)/2, 2).';

    y0 = kron(y2, ones(length(x{1}),1));
end
