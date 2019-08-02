

% benney = @(a, b, params) -21i .* a + 8./15 .* params.Re .* a.^2 - ...
%     1./3 .* (a + b) .* (2 .* cot(params.beta) + 1./params.C .* (a.^2 + b.^2));
benney = @(a, b, params) -2.0.*(-4./15.*params.C.*params.Re.*a.^2 + ...
    1./3.*params.C.*a.^2.*cot(params.beta) + 1i.*params.C.*a + ...
    1./3.*params.C.*b.^2.*cot(params.beta) + 1./6.*a.^4 + 1./3.*a.^2.*b.^2 + ...
    1./6.*b.^4)./params.C;

wibl1_1 = @(a, b, params) -(68.*1i.*params.C.*params.Re.*a.*tan(params.beta) + 28.*params.C.*params.Re.*b.^2 + 105.*params.C.*tan(params.beta) + 14.*params.Re.*a.^2.*b.^2.*tan(params.beta) + 14.*params.Re.*b.^4.*tan(params.beta) + 1i.*sqrt(592.*params.C.^ ...
    2.*params.Re.^2.*a.^2 + 3808.*1i.*params.C.^2.*params.Re.^2.*a.*b.^2.*cot(params.beta) - 784.*params.C.^2.*params.Re.^2.*b.^4.*cot(params.beta).^2 + 11760.*params.C.^2.*params.Re.*a.^2.*cot(params.beta) + 21000.*1i.*params.C.^2....
    .*params.Re.*a + 5880.*params.C.^2.*params.Re.*b.^2.*cot(params.beta) - 11025.*params.C.^2 + 1904.*1i.*params.C.*params.Re.^2.*a.^3.*b.^2 - 784.*params.C.*params.Re.^2.*a.^2.*b.^4.*cot(params.beta) + 1904.*1i.*params.C.*params.Re.^2.*a.*b...
    .^4 - 784.*params.C.*params.Re.^2.*b.^6.*cot(params.beta) + 5880.*params.C.*params.Re.*a.^4 + 8820.*params.C.*params.Re.*a.^2.*b.^2 + 2940.*params.C.*params.Re.*b.^4 - 196.*params.Re.^2.*a.^4.*b.^4 - 392.*params.Re.^2.*a.^2.*b.^6 ...
    - 196.*params.Re.^2.*b.^8).*tan(params.beta))./(84.*params.C.*params.Re.*tan(params.beta));

wibl1_2 = @(a, b, params) -(68.*1i.*params.C.*params.Re.*a.*tan(params.beta) + 28.*params.C.*params.Re.*b.^2 + 105.*params.C.*tan(params.beta) + 14.*params.Re.*a.^2.*b.^2.*tan(params.beta) + 14.*params.Re.*b.^4.*tan(params.beta) - 1i.*sqrt(592.*params.C.^...
    2.*params.Re.^2.*a.^2 + 3808.*1i.*params.C.^2.*params.Re.^2.*a.*b.^2.*cot(params.beta) - 784.*params.C.^2.*params.Re.^2.*b.^4.*cot(params.beta).^2 + 11760.*params.C.^2.*params.Re.*a.^2.*cot(params.beta) + 21000.*1i.*params.C.^2....
    .*params.Re.*a + 5880.*params.C.^2.*params.Re.*b.^2.*cot(params.beta) - 11025.*params.C.^2 + 1904.*1i.*params.C.*params.Re.^2.*a.^3.*b.^2 - 784.*params.C.*params.Re.^2.*a.^2.*b.^4.*cot(params.beta) + 1904.*1i.*params.C.*params.Re.^2.*a.*b...
    .^4 - 784.*params.C.*params.Re.^2.*b.^6.*cot(params.beta) + 5880.*params.C.*params.Re.*a.^4 + 8820.*params.C.*params.Re.*a.^2.*b.^2 + 2940.*params.C.*params.Re.*b.^4 - 196.*params.Re.^2.*a.^4.*b.^4 - 392.*params.Re.^2.*a.^2.*b.^6 ...
    - 196.*params.Re.^2.*b.^8).*tan(params.beta))./(84.*params.C.*params.Re.*tan(params.beta));

wibl1A_1 = @(a, b, params) (68.*params.C.*params.Re.*a + 28.*1i.*params.C.*params.Re.*b.^2.*cot(params.beta) - 105.*1i.*params.C + 14.*1i.*params.Re.*a.^2.*b.^2 + 14.*1i.*params.Re.*b.^4 + sqrt(592.*params.C.^2.*params.Re.^2.*a.^2 + 3808.*1i.*params.C.^2.*...
    params.Re.^2.*a.*b.^2.*cot(params.beta) - 784.*params.C.^2.*params.Re.^2.*b.^4.*cot(params.beta).^2 + 11760.*params.C.^2.*params.Re.*a.^2.*cot(params.beta) + 21000.*1i.*params.C.^2.*params.Re.*a + 5880.*params.C.^2.*params.Re.*b.^2.*cot(params.beta) ...
    - 11025.*params.C.^2 + 1904.*1i.*params.C.*params.Re.^2.*a.^3.*b.^2 - 784.*params.C.*params.Re.^2.*a.^2.*b.^4.*cot(params.beta) + 1904.*1i.*params.C.*params.Re.^2.*a.*b.^4 - 784.*params.C.*params.Re.^2.*b.^6.*cot(...
    params.beta) + 5880.*params.C.*params.Re.*a.^4 + 8820.*params.C.*params.Re.*a.^2.*b.^2 + 2940.*params.C.*params.Re.*b.^4 - 196.*params.Re.^2.*a.^4.*b.^4 - 392.*params.Re.^2.*a.^2.*b.^6 - 196.*params.Re.^2.*b.^8))./(84.*params.C...
    .*params.Re.*a);


wibl1A_2 = @(a, b, params) (68.*params.C.*params.Re.*a + 28.*1i.*params.C.*params.Re.*b.^2.*cot(params.beta) - 105.*1i.*params.C + 14.*1i.*params.Re.*a.^2.*b.^2 + 14.*1i.*params.Re.*b.^4 - sqrt(592.*params.C.^2.*params.Re.^2.*a.^2 + 3808.*1i.*params.C.^2.*...
    params.Re.^2.*a.*b.^2.*cot(params.beta) - 784.*params.C.^2.*params.Re.^2.*b.^4.*cot(params.beta).^2 + 11760.*params.C.^2.*params.Re.*a.^2.*cot(params.beta) + 21000.*1i.*params.C.^2.*params.Re.*a + 5880.*params.C.^2.*params.Re.*b.^2.*cot(params.beta) ...
    - 11025.*params.C.^2 + 1904.*1i.*params.C.*params.Re.^2.*a.^3.*b.^2 - 784.*params.C.*params.Re.^2.*a.^2.*b.^4.*cot(params.beta) + 1904.*1i.*params.C.*params.Re.^2.*a.*b.^4 - 784.*params.C.*params.Re.^2.*b.^6.*cot(...
    params.beta) + 5880.*params.C.*params.Re.*a.^4 + 8820.*params.C.*params.Re.*a.^2.*b.^2 + 2940.*params.C.*params.Re.*b.^4 - 196.*params.Re.^2.*a.^4.*b.^4 - 392.*params.Re.^2.*a.^2.*b.^6 - 196.*params.Re.^2.*b.^8))./(84.*params.C...
    .*params.Re.*a);

fig1 = figure(1);
fig2 = figure(2);

myfunction(benney, 0.3, 0.3, fig1, fig2)
myfunction(wibl1_1, 0.3, 0.3, fig1, fig2)

function myfunction(func, a0, b0, fig1, fig2)
    params = struct();
    params.Re = 5;
    params.beta = 7*pi/8;
    params.C = 0.01;
    
    aCritical = fsolve(@(x) real(func(x, 0, params)), a0);
    bCritical = fsolve(@(x) real(func(0, x, params)), b0);
    
    b = 0;
    a = linspace(0,1);
    
    figure(fig1)
    hold on
    plot(a, real(func(a, b, params)));
    
    figure(fig2);
    hold on
    [A, B] = meshgrid(linspace(0,1),linspace(0,1));
    contour(A, B, real(func(A, B, params)),[0, 0]);
    
    fprintf("Critical alpha: %f\nCritical beta: %f\n", aCritical, bCritical)
end