function filename = makeFilename(prefix, ivpArguments, timePointsArguments, timeStepperArguments)
    
    timeStepperArguments = rmfield(timeStepperArguments, 'outputOpt');
    
    filename = prefix + struct2str(ivpArguments) + ...
        struct2str(timePointsArguments) + struct2str(timeStepperArguments);
end