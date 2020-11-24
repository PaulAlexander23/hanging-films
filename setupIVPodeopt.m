function odeopt = setupIVPodeopt(odejac, domain, params, method)
    odeopt = odeset();
    
    if method == "finite-difference"
        odeopt = odeset(odeopt ...
            , 'Jacobian', @(t, y) odejac(domain, y, params) ...
            );
    end
end
