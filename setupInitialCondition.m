function y0 = setupInitialCondition(domain, interface, method)
    
    y0 = interface(domain.x);
    
    if method == "pseudo-spectral"
        y0 = domain.fft(y0);
    end

    y0 = domain.reshapeToVector(y0);
end
