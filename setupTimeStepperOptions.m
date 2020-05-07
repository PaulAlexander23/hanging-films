function odeopt = setupTimeStepperOptions(odeoptDefault, odeoptIVP)
    
    odeopt = odeoptDefault;
    
    if isfield(odeoptIVP, 'Jacobian')
        odeopt.Jacobian = odeoptIVP.Jacobian;
    end
end
