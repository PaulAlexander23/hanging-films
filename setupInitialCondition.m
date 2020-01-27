function y0 = setupInitialCondition(model, domain, interface, method)
    
    y0 = interface(domain.x);
    
    if method == "pseudo-spectral"
        if model == "benney"
            y0 = domain.fft(y0);
        elseif model == "wibl1"
            y0 = [domain.fft(y0(1:end/2,:)); ...
                domain.fft(y0(1+end/2:end,:))];
        end
    end
end