function odeopt = setupTimeStepperOptions(odeoptDefault, odeoptIVP)
    
    odeopt = odeoptDefault;
    
    odeopt = odeset(odeoptIVP, odeopt);
end