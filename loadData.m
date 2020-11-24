function solution = loadData(ivpArguments, timePointsArguments, suffix)
    if nargin < 3
        suffix = "";
    end
    
    filename = makeFilename("data", ivpArguments, timePointsArguments) ...
        + suffix;
    
    load(filename + '.mat', 'solution');
end
