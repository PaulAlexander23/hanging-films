function tests = testJacobian()
    tests = functiontests(localfunctions);
end

function testBenneyJacobian(testCase)
    addpath discretisationMethods

    diffDegrees = [1, 0; 0, 1; 2, 0; 0, 2; 4, 0; 0, 3; 0, 4; 2, 2; 2, 1; 1, 2; 3, 0]';
    domain = FDDomain(setupX(1, 1, 32, 32), diffDegrees, 4);
    y = 1 + 0.25 * cos(2*pi*domain.x{1}) + 0.25 * cos(2*pi*domain.x{2});
    params = struct('theta', 7/8*pi, 'Re', 1, 'C', 0.01);

    yVector = domain.reshapeToVector(y);

    F = @(u) Fbenney(domain, u, params);

    expected = jacobianNumerical(F, yVector);
    actual = jbenney2d(domain, yVector, params);

    verifyEqual(testCase, actual, expected, 'RelTol', 1e-7)
end

% function testBenneyJacobianPS(testCase)
%     addpath discretisationMethods
%
%     domain = PSDomain(setupX(1,1,2^6,2^6), nan, true, 2);
%     y = 1 + 0.25 * cos(2*pi*domain.x{1}) + 0.25 * cos(2*pi*domain.x{2}');
%     params = struct('theta', 7/8*pi, 'Re', 1, 'C', 0.01);
%     y = domain.fft(y);
%
%     yVector = domain.reshapeToVector(y);
%
%     F = @(u) Fbenney(domain, u, params);
%
%     expected = jacobianNumericalPS(F, yVector);
%     actual = jbenney2d(domain, yVector, params);
%
%     verifyEqual(testCase, actual, expected, 'RelTol', 1e-7)
% end

function testWIBL1Jacobian(testCase)
    addpath discretisationMethods

    diffDegrees = [1, 0; 0, 1; 2, 0; 0, 2; 4, 0; 0, 3; 0, 4; 2, 2; 2, 1; 1, 2; 3, 0]';
    domain = FDDomain(setupX(1, 1, 32, 32), diffDegrees, 4);
    y = 1 + 0.25 * cos(2*pi*domain.x{1}) + 0.25 * cos(2*pi*domain.x{2});
    params = struct('theta', 7/8*pi, 'Re', 1, 'C', 0.01);

    y = [y; 2 / 3 * y];
    yVector = [domain.reshapeToVector(y(1:end/2, :, :)); ...
        domain.reshapeToVector(y(1+end/2:end, :, :))];

    F = @(u) FWIBL1(domain, u, params);

    expected = jacobianNumerical(F, yVector);
    actual = jwibl1(domain, yVector, params);

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

    f = fwibl1(domain, y, params);

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
