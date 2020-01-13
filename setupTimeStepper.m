function timeStepper = setupTimeStepper(method, model, domain, params, AbsTol, debug)
    
    odeopt = setupTimeStepperOptions(method, model, domain, params, AbsTol, debug);
    
    if method == "finite-difference"
        timeStepper = @(odefun, t, y0) ode15s(odefun, t, y0, odeopt);
    elseif method == "pseudo-spectral"
        timeStepper = @(odefun, t, y0) ode15s(odefun, t, y0, odeopt);
    end
end