function y0 = isingle(x,A)
    
    if nargin < 2, A = 1e-2; end
    y0 = 1 - A * (cos(2*pi/x{1}(end) * x{1}) + cos(2*pi/x{2}(end) * x{2}'));
    
end