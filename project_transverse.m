function project_transverse(t, x, y)
    
    y = reshape(y, [cellfun(@length, x)', length(t)]);
    
    fh = fft2(y);
    
    fh(2:end,:) = 0;
    
    fh = ifft2(fh);
    
    plot(t, energy(fh, cellfun(@length, x)))
    
end