function filename = saveData(solution, ivpArguments, timePointsArguments, timeStepperArguments)
    filename = makeFilename("data", ivpArguments, timePointsArguments);

    if strlength(filename) > 255
        filename = extractBefore(filename, 256);
    end
    
    filename = ensureUnique(filename);
    
    save(filename + '.mat', 'solution', 'ivpArguments', 'timePointsArguments', ...
        'timeStepperArguments', '-v7.3');
end
