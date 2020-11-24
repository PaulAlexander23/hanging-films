function y0 = iloadInterpWIBL1(xq, filename)

    load(filename,'x','h');
    
    h0 = periodicInterp2(x{1}, x{2}, h(1:end/2,:), xq{1}, xq{2}, 'spline');
    f10 = periodicInterp2(x{1}, x{2}, h(1+end/2:end,:), xq{1}, xq{2}, 'spline');

    y0 = [h0; f10];
end
