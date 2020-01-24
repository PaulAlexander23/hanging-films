function odeopt = setupTimeStepperOptions(odeoptDefault, outputOpt, odeoptIVP)
    
    odeopt = odeoptDefault;
    
    odeopt = odeset(outputOpt, odeopt);
    
    odeopt = odeset(odeoptIVP, odeopt);
end