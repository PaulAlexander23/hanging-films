function dYdt = fwibl2STF(domain, Y, params)

    epsilon = 1;

    Y = domain.reshapeToDomain(Y);

    h = Y(1:end/3, :, :);
    F1 = Y(end/3+1:2*end/3, :, :);
    F2 = Y(2*end/3+1:end, :, :);

    dhdx = domain.diff(h, [1; 0]);
    d2hdx2 = domain.diff(h, [2; 0]);
    d3hdx3 = domain.diff(h, [3; 0]);
    dhdz = domain.diff(h, [0; 1]);
    d2hdz2 = domain.diff(h, [0; 2]);
    d3hdz3 = domain.diff(h, [0; 3]);
    d2hdxdz = domain.diff(h, [1; 1]);
    d3hdx2z = domain.diff(h, [2; 1]);
    d3hdxz2 = domain.diff(h, [1; 2]);

    dF1dx = domain.diff(F1, [1; 0]);
    d2F1dx2 = domain.diff(F1, [2; 0]);
    dF1dz = domain.diff(F1, [0; 1]);
    d2F1dz2 = domain.diff(F1, [0; 2]);
    d2F1dxdz = domain.diff(F1, [1; 1]);
    
    dF2dx = domain.diff(F2, [1; 0]);
    d2F2dx2 = domain.diff(F2, [2; 0]);
    dF2dz = domain.diff(F2, [0; 1]);
    d2F2dz2 = domain.diff(F2, [0; 2]);
    d2F2dxdz = domain.diff(F2, [1; 1]);
    
    dhdt = - dF1dx - dF2dz;

    dF1dt = ...
        m(-17*m(F1,dF1dx) ...
		- 8*epsilon*m(F1,dF2dz) ...
		- 9*epsilon*m(F2,dF1dz),(7*h),[1,-1]) ...
		+ m(9*m(F1,dhdx,[2,1])/7 ...
		+ 9*epsilon*m(m(F1,F2),dhdz)/7,h,[1,-2]) ...
		+ 9*epsilon*d2F1dx2/(2*params.Re) ...
		+ epsilon*d2F1dz2/params.Re ...
		+ m(- 6*epsilon*m(F1,d2hdx2)/(params.Re) ...
		- 23*epsilon*m(F1,d2hdz2)/(16*params.Re) ...
		- 9*epsilon*m(dhdx,dF1dx)/(2*params.Re) ...
		- epsilon*m(dhdz,dF1dz)/(params.Re) ...
		+ m(4*epsilon*m(F1,dhdx,[1,2])/(params.Re) ...
		+ 3*epsilon*m(F1,dhdz,[1,2])/(4*params.Re) ...
		- 5*cot(params.theta)*m(h,dhdx)/(3*params.Re) ...
		+ 5*h/(3*params.Re*epsilon) ...
		- 5*m(F1,h,[1,-2])/(2*params.Re*epsilon) ...
		+ 5*epsilon^2*m(h,d3hdx3)/(6*params.C*params.Re) ...
		+ 5*epsilon^2*m(h,d3hdxz2)/(6*params.C*params.Re);

    dF2dt = ...
		m(-9*m(F1,dF2dx) ...
		- 8*m(F2,dF1dx) ...
		+ m(9*m(m(F1,F2),dhdx)/(7) ...
		+ 7*d2F1dxdz/(2*params.Re) ...
		+ m(- 73*m(F1,d2hdxdz)/(16*params.Re) ...
		- 43*m(dhdx,dF1dz)/(16*params.Re) ...
		- 13*m(dhdz,dF1dx)/(16*params.Re) ...
		+ m(13*m(m(F1,dhdx),dhdz)/(4*params.Re) ...
		- 5*cot(params.theta)*m(h,dhdz)/(3*params.Re*epsilon) ...
		- 5*m(F2,h,[1,-2])/(2*params.Re*epsilon^2) ...
		+ 5*epsilon*m(h,d3hdz3)/(6*params.C*params.Re) ...
		+ 5*epsilon*m(h,d3hdx2z)/(6*params.C*params.Re);

    dYdt = cat(1, dhdt, dF1dt, dF2dt);

    dYdt = domain.reshapeToVector(dYdt);

    function c = m(a,b,pow)
        if nargin<3, pow = [1,1]; end

        c = domain.multiply(a,b,pow);
    end
end