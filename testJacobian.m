function tests = testJacobian()
    tests = functiontests(localfunctions);
end

function testBenneyJacobian(testCase)
    addpath discretisationMethods

    domain = tFDDomain();
    y = 1 + 0.25 * cos(2*pi*domain.x{1}) + 0.25 * cos(2*pi*domain.x{2});
    params = struct('theta', 7/8*pi, 'Re', 1, 'C', 0.01);

    yVector = domain.reshapeToVector(y);

    F = @(u)fbenney2d(domain,u,params);

    expected = jacobianNumerical(F, yVector);
    actual = jbenney2d(domain, yVector, params);

   verifyEqual(testCase, actual, expected, 'RelTol', 1e-6)
end

function testWIBL1Jacobian(testCase)
    addpath discretisationMethods

    domain = tFDDomain();
    y = 1 + 0.25 * cos(2*pi*domain.x{1}) + 0.25 * cos(2*pi*domain.x{2});
    params = struct('theta', 7/8*pi, 'Re', 1, 'C', 0.01);

    y = [y; 2 / 3 * y];
    yVector = [domain.reshapeToVector(y(1:end/2, :, :)); ...
        domain.reshapeToVector(y(1+end/2:end, :, :))];

    F = @(u)fwibl1STF(domain, u, params);

    expected = jacobianNumerical(F, yVector);
    actual = jwibl1STF(domain, yVector, params);

    verifyEqual(testCase, actual, expected, 'AbsTol', 1e-3, 'RelTol', 1e-3)
end

function testHybridJacobian(testCase)
    addpath discretisationMethods

    domain = tFDDomain();
    y = 1 + 0.25 * cos(2*pi*domain.x{1}) + 0.25 * cos(2*pi*domain.x{2});
    params = struct('theta', 7/8*pi, 'Re', 1, 'C', 0.01, 'epsilon', 1, ...
        'delta', 1);

    y = [y; 2 / 3 * y];
    yVector = [domain.reshapeToVector(y(1:end/2, :, :)); ...
        domain.reshapeToVector(y(1+end/2:end, :, :))];

    F = @(u)fhybrid(domain, u, params);

    expected = jacobianNumerical(F, yVector);
    actual = jhybrid(domain, yVector, params);

    verifyEqual(testCase, actual, expected, 'AbsTol', 1e-3, 'RelTol', 1e-3)
end

function testHybridJacobianInplace(testCase)
    addpath discretisationMethods

    domain = tFDDomain();
    y = 1 + 0.25 * cos(2*pi*domain.x{1}) + 0.25 * cos(2*pi*domain.x{2});
    params = struct('theta', 7/8*pi, 'Re', 1, 'C', 0.01, 'epsilon', 1, ...
        'delta', 1);

    y = [y; 2 / 3 * y];
    yVector = [domain.reshapeToVector(y(1:end/2, :, :)); ...
        domain.reshapeToVector(y(1+end/2:end, :, :))];

    F = @(u)fhybrid(domain, u, params);

    expected = jacobianNumerical(F, yVector);
    [~, actual] = fhybrid(domain, yVector, params);

    verifyEqual(testCase, actual, expected, 'AbsTol', 1e-3, 'RelTol', 1e-3)
end


function fVector = Fbenney(domain, yVector, params)
    y = domain.reshapeToDomain(yVector);
    f = fbenney2d(domain, y, params);
    fVector = domain.reshapeToVector(f);
end

function fVector = FWIBL1(domain, yVector, params)
    y = [domain.reshapeToDomain(yVector(1:end/2)); ...
        domain.reshapeToDomain(yVector(1+end/2:end))];

    f = fwibl1STF(domain, y, params);

    fVector = [domain.reshapeToVector(f(1:end/2, :, :)); ...
        domain.reshapeToVector(f(1+end/2:end, :, :))];
end

function J = jacobianNumerical(F, u)
    N = length(u);
    epsilon = 1e-6;
    e = eye(N);
    J = sparse(N, N);

    Fu = F(u);

    for j = 1:N
        J(:, j) = (F(e(:, j)*epsilon+u) - Fu) / epsilon;
    end
end

function J = jacobianNumericalPS(F, u)
    N = length(u);
    epsilon = 1e-5;
    e = eye(N);
    J = zeros(N, N/2);

    Fu = F(u);

    for j = 1:N / 2
        epsilonRescaled = epsilon;
        J(:, j) = (F(e(:, j)*epsilonRescaled+u) - Fu) / epsilonRescaled;
    end
end

function domain = tFDDomain()

    diffDegrees = [1, 0; 0, 1; 2, 0; 0, 2; 4, 0; 0, 3; 0, 4; 2, 2; 2, 1; 1, 2; 3, 0]';
    domain = FDDomain(setupX(1, 1, 8, 8), diffDegrees, 2);
end
