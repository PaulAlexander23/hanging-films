function y = postprocess(method, domain, y)

    y = permute(y, [1, 3, 2]);
    y = domain.reshapeToDomain(y);

    % Optional postprocessing dependant on method.
    if method == "pseudo-spectral"
        y = domain.ifft(y);
    end
end
