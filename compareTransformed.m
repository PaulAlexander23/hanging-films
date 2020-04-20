
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

t = solution.t;

domain = setupDomain(domainArguments);

y0 = domain.reshapeToVector(diffusionTransformation(solution.y(:,:,1)));
odeFunction = @fbenney2dTransformedP;
odeFunction = matFuncToVecFunc(odeFunction);
odeFunction = @(t, y) odeFunction(domain, y, params);

odeoptIVP = odeset(...odeopt ...
...    ,'Jacobian', @(t, y) jbenney2d(domain, y, args.params) ...
);

timeStepper = setupTimeStepper(timeStepperArguments, odeoptIVP);

[y, t, timeTaken] = iterateTimeStepper(odeFunction, t, y0, timeStepper);

y = inverseDiffusionTransformation(y);
y = domain.reshapeToDomain(permute(y, [1, 3, 2]));

transformedSolution = struct('domain', domain, 't', t, 'y', y, ...
    'timeTaken', timeTaken);

%%

expected = solution.y;
actual = transformedSolution.y;

difference = actual - expected;
max(abs(difference(:,:,end)),[],'all')
max(abs(difference(:,:,end)./expected(:,:,end)),[],'all')
figure;
surf(expected(:,:,end));
figure;
surf(actual(:,:,end));
figure;
surf(difference(:,:,end));
