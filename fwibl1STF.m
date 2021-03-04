function dYdt = fwibl1STF(domain, Y, params)

    epsilon = 1;

    Y = domain.reshapeToDomain(Y);

    h = Y(1:end/2, :, :);
    F1 = Y(end/2+1:end, :, :);

    dhdx = domain.diff(h, [1; 0]);
    d2hdx2 = domain.diff(h, [2; 0]);
    d3hdx3 = domain.diff(h, [3; 0]);
    d2hdz2 = domain.diff(h, [0; 2]);
    d3hdxz2 = domain.diff(h, [1; 2]);

    dF1dx = domain.diff(F1, [1; 0]);
    
    P = 2 * h * cot(params.theta) ...
        - epsilon^2 * (d2hdx2 + d2hdz2) / params.C;

    F2 = - epsilon * m(h, domain.diff(P, [0; 1]), [3,1]) / 3;
    
    dF2dz = domain.diff(F2, [0; 1]);
    
    dhdt = - dF1dx - dF2dz;

    dF1dt = ...
        m(-17*m(F1,dF1dx),7*h,[1,-1]) ...
		+ m(9*m(F1,dhdx,[2,1]),h,[1,-2])/7 ...
		- 5*cot(params.theta)*m(h,dhdx)/(3*params.Re) ...
		+ 5*h/(3*params.Re*epsilon) ...
		- 5*m(F1,h,[1,-2])/(2*params.Re*epsilon) ...
		+ 5*epsilon^2*m(h,d3hdx3)/(6*params.C*params.Re) ...
		+ 5*epsilon^2*m(h,d3hdxz2)/(6*params.C*params.Re);

    dYdt = cat(1, dhdt, dF1dt);

    dYdt = domain.reshapeToVector(dYdt);

    function c = m(a,b,pow)
        if nargin<3, pow = [1,1]; end

        c = domain.multiply(a,b,pow);
    end
end
