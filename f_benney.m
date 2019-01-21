function F = f_benney(H,L,params)
    
    delta = params(1);
    theta = params(2);
    Re = params(3);
    We = params(4);
    C = params(5);
    
    
    e1 = zeros(1,1,2);
    e1(1,1,1) = 1;
    Hx = diff_ps(H,1,L(1));
    
    P = H * cot(theta) - We * R(H, L) - 1/2 * 1/C * lap(H, L);
    F = div(...
        (2/3 * H.^3 + 8/15 * Re * delta * H .^ 6 .* Hx) .* e1...
        - 2/3 * delta * H.^3 .* grad(P, L),...
        L);
end