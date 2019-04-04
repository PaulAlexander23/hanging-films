function E = energy(x, y)

    dx = cellfun(@(x) x(1), x);
    
    E = sum((y - 1).^2*prod(dx),[1,2]);
    
    E = squeeze(E);
end
