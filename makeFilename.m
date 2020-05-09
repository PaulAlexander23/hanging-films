function filename = makeFilename(prefix, ivpArguments, timePointsArguments)
    ivpArguments = rmfield(ivpArguments, 'interface');
    
    filename = prefix + struct2str(ivpArguments) + ...
        struct2str(timePointsArguments);
end
