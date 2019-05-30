function y0 = iload_with_pert(x,filename,alpha)
    if nargin < 3, alpha = 1e-2; end
    load(filename,'h');
    y0 = h + alpha * cos(2*pi*(x{1}/x{1}(end) + x{2}/x{2}(end)));
end
