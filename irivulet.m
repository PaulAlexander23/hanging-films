function y0 = irivulet(x)
    
    A = 2e-1;
    r = 0.05;
    y0 = 1 + A * (-r*cos(2*pi/x{1}(end) * x{1}) - cos(2*pi/x{2}(end) * x{2}'));
    
end