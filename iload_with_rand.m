function y0 = iload_with_rand(x,filename,alpha)
    if nargin < 3, alpha = 1e-2; end
    load(filename,'h');
    y0 = h + irand(x,alpha,5) - 1;
end
