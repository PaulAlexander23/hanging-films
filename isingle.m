function y0 = isingle(x)
    
    A = 1e-2;
    r = 1;
    y0 = 1 + A * (-r*cos(2*pi/x{1}(end) * x{1}) - cos(2*pi/x{2}(end) * x{2}'));
    
end