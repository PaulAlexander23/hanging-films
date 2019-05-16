function E = energy(x, y)

    dx = cellfun(@(x) x(1), x);
    xl = cellfun(@(x) x(end), x);
    
    E = sum((y - 1).^2*prod(dx)/prod(xl),[1,2]);
    
    E = squeeze(E);
end
