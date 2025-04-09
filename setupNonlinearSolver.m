function nonlinearSolver = setupNonlinearSolver(args, odeoptBVP)
    
    odeopt = odeoptDefault;
    
    if isfield(odeoptBVP, 'Jacobian')
        odeopt.Jacobian = odeoptBVP.Jacobian;
    end
    
    nonlinearSolver = @(odefun, y0) args.nonlinearSolver(odefun, y0, odeopt);
end
