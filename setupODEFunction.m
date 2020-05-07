function odeFunction = setupODEFunction(model, domain, params)
    if model == "benney"
        f = @fbenney2d;
    elseif model == "wibl1"
        f = @fwibl1;
    end

    odeFunction = @(t, y) f(domain, y, params);
end
