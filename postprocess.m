function y = postprocess(method, domain, y)

    y = permute(y, [1, 3, 2]);
    y = domain.reshapeToDomain(y);

    if method == "pseudo-spectral"
        y = domain.ifft(y);
    end
end
