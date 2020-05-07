function timeStepper = setupTimeStepper(args, odeoptIVP)
    
    odeopt = setupTimeStepperOptions(args.odeopt, odeoptIVP);
    
    timeStepper = @(odefun, t, y0) args.timeStepper(odefun, t, y0, odeopt);
end
