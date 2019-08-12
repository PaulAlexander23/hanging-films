function dYdthat = fwibl1(domain, Yhat, params)
    yhat = Yhat(1:end/2, :, :);
    F1hat = Yhat(end/2+1:end, :, :);
    
    Phat = 2 * yhat * cot(params.theta) - 1 / params.C * (domain.diff(yhat, [2; 0]) + domain.diff(yhat, [0; 2]));

    F2hat = - domain.multiply(yhat, domain.diff(Phat, [0; 1]), [3, 1]) / 3;

    dydthat = -domain.diff(F1hat, [1; 0]) - domain.diff(F2hat, [0; 1]);

    dF1dthat =  domain.multiply( ...
        9 * domain.multiply( ...
        domain.multiply(domain.diff(yhat, [1; 0]), F1hat, [1, 2]), ...
        yhat, [1, -1]) - 17 * domain.multiply(F1hat, domain.diff(F1hat, [1; 0]), [1, 1]), ...
        (7 * yhat), [1, -1]) + ...
        (10 * yhat - 15 * domain.multiply(F1hat, yhat, [1, -2]) - ...
        5 * domain.multiply(yhat, domain.diff(Phat, [1; 0]), [1, 1])) / (6 * params.Re);

    dYdthat = cat(1, dydthat, dF1dthat);
end
