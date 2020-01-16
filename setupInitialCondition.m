function y0 = setupInitialCondition(model, domain, interface, method)
    if model == "benney"
        y0 = interface(domain.x);
    elseif model == "wibl1"
        y0 = interface(domain.x);
        F0 = 2/3 + 0*y0;
        y0 = [y0; F0];
    end
    
    if method == "pseudo-spectral"
        if model == "benney"
            y0 = domain.fft(y0);
        elseif model == "wibl1"
            y0 = [domain.fft(y0(1:end/2,:)); ...
                domain.fft(y0(1+end/2:end,:))];
        end
    end
end