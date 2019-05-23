function y0 = iload_with_pert(x,filename,alpha)
    if nargin < 3, alpha = 1e-2; end
    load(filename,'h');
    y0 = h + alpha * cos(x{1} + x{2});
end
