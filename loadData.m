function solution = loadData(ivpArguments, timePointsArguments, timeStepperArguments, suffix)
    if nargin < 4
        suffix = "";
    end
    
    filename = makeFilename("data", ivpArguments, timePointsArguments, ...
        timeStepperArguments) + suffix;
    
    load(filename + '.mat', 'solution');
end