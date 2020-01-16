function saveData(solution, ivpArguments, timePointsArguments, timeStepperArguments)
    filename = makeFilename("data", ivpArguments, timePointsArguments, ...
        timeStepperArguments);
    
    filename = ensureUnique(filename);
    
    save(filename + '.mat', 'solution', 'ivpArguments', 'timePointsArguments', ...
        'timeStepperArguments');
end