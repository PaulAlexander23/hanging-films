function odeFunction = setupODEFunction(model, domain, params)
    if model == "benney"
        odeFunction = @(t, y) fbenney2d(domain, y, params);
    elseif model == "wibl1"
        odeFunction = @(t, y) fwibl1(domain, y, params);
    end
end