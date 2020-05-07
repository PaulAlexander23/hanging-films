function filename = makeFilename(prefix, ivpArguments, timePointsArguments)
    
    filename = prefix + struct2str(ivpArguments) + ...
        struct2str(timePointsArguments);
end
