function y0 = iloadInterp(xq, filename)

    load(filename,'x','h');
    
    y0 = periodicInterp2(x{1}, x{2}, h, xq{1}, xq{2}, 'spline');
end
