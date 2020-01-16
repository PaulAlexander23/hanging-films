function y = postprocess(method, model, domain, y)
        if method == "pseudo-spectral"
            if model == "benney"
                y = domain.ifft(y);
            elseif model == "wibl1"
                y = [domain.ifft(y(1:end/2,:,:)); ...
                    domain.ifft(y(1+end/2:end,:,:))];
            end
        end
    
end