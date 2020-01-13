function odeopt = setupTimeStepperOptions(odeoptDefault, odeoptIVP, debug)
    
    odeopt = odeoptDefault;
    
    odeopt = odeset(odeoptIVP, odeopt);
    
    if debug
        odeopt = odeset(odeopt, ...
            'OutputFcn', 'odeprint', ...
            'OutputSel', 1, ...
            'Stats', 'on');
    end
end