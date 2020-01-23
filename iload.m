function y0 = iload(x, filename)
    load(filename,'h');
    
    if length(x{1}) == size(h,1) && length(x{2}) == size(h,2)
        y0 = h;
    else
        error('Interface loaded has incorrect size. Expected size: [%g, %g], Actual: [%g, %g]', ...
            length(x{1}), length(x{2}), size(h,1), size(h,2))
    end
end
