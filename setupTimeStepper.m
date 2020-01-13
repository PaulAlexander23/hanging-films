function timeStepper = setupTimeStepper(args, odeoptIVP, debug)
    
    odeopt = setupTimeStepperOptions(args.odeopt, odeoptIVP, debug);
    
    timeStepper = @(odefun, t, y0) args.timeStepper(odefun, t, y0, odeopt);
end