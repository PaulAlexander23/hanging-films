
addpath('discretisationMethods/')

model = 'benney';
params = struct('theta', 7*pi/8, 'Re', 1, 'C', 0.01);
tFinal = 0.01;
interface = @(x)icos(x,0.1,0.1);
method = 'finite-difference';
AbsTol = 1e-6;

domainArguments = struct('xLength', 2*pi, 'yLength', 2*pi, 'xN', 2^5, ...
    'yN', 2^5, 'method', method);
ivpArguments = struct('domainArguments',domainArguments,'params',params,'method',method,'model',model,'interface',interface);
timePointsArguments = struct('tStep', 0.0000001, 'tFinal', tFinal);
odeoptDefault = odeset( ...
    ...'Vectorized', 'on', ...
    ...'BDF','on', ...
    'AbsTol', AbsTol ...
    ...'MaxStep', 5e-6 ...
    ...'InitialStep', 1e-3 ...
    );
odeoptOutput = odeset();
addpath ../../time-stepping-methods
timeStepperArguments = struct('timeStepper', @ab1, ...
    'odeopt', odeoptDefault, 'outputOpt', odeoptOutput);

solution = solveIVP(ivpArguments, timePointsArguments, timeStepperArguments);


%%

domain = setupDomain(domainArguments);

fbenney2dExplicitVec = matFuncToVecFunc(@fbenney2dExplicit);
fbenney2dImplicitVec = matFuncToVecFunc(@fbenney2dImplicit);
explicitOdefun = @(t, y) fbenney2dExplicitVec(domain, y, params);

t = linspace(0, solution.t(end), 201)';

y0 = domain.reshapeToVector(solution.y(:,:,1));

timeStepper = @bdf1si;

options = struct('optimmethod', @(fun, x0) fsolve(fun, x0, ...
    optimoptions('fsolve', 'Display', 'off', 'SpecifyObjectiveGradient', true)));

tic
y = timeStepper(explicitOdefun, @(t, y) implicitOdefun(t, y, domain, params), t, y0, options);
timeTaken = toc;

y = domain.reshapeToDomain(y);

semiImplicitSolution = struct('domain', domain, 't', t, 'y', y, ...
    'timeTaken', timeTaken);

%%

solution.timeTaken
semiImplicitSolution.timeTaken
expected = solution.y;
actual = semiImplicitSolution.y;

difference = actual(:,:,1:500:end) - expected;
max(abs(difference(:,:,end)),[],'all')
max(abs(difference(:,:,end)./expected(:,:,end)),[],'all')
figure;
surf(expected(:,:,end));
figure;
surf(actual(:,:,end));
figure;
surf(difference(:,:,end));

%% Functions

function [F, J] = implicitOdefun(t, y, domain, params)
    [f, J] = fbenney2dImplicit(domain, domain.reshapeToDomain(y),params);
    
    F = domain.reshapeToVector(f);
end
