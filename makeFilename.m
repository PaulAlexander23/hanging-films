function filename = makeFilename(prefix, ivpArguments, timePointsArguments, timeStepperArguments)
    
    filename = prefix + struct2str(ivpArguments) + ...
        struct2str(timePointsArguments);% + struct2str(timeStepperArguments);
end