function y = postprocess(method, model, domain, y)

    y = permute(y, [1, 3, 2]);
    y = domain.reshapeToDomain(y);

    if method == "pseudo-spectral"
        if model == "benney"
            y = domain.ifft(y);
        elseif model == "wibl1"
            y = [domain.ifft(y(1:end/2,:,:)); ...
                domain.ifft(y(1+end/2:end,:,:))];
        end
    end
end
