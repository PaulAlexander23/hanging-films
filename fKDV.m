function F = fKDV(u,l,params)
    
    eta = params(1);

    F = diff_ps(u',3,l) - 6 * eta * diff_ps(u'.^2/2,1,l);
    
end