function odeFunction = setupODEFunction(model, domain, params)
    if model == "benney"
        f = @fbenney2d;
    elseif model == "wibl1"
        f = @fwibl1;
    end

    odeFunctionVector = matFuncToVecFunc(f);
    odeFunction = @(t, y) odeFunctionVector(domain, y, params);
end
