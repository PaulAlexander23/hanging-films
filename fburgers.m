function F = fburgers(domain,y,params)
    F = -domain.multiply(y, domain.diff(y, 1)) + ...
        params.nu * domain.diff(y, 2);
end
