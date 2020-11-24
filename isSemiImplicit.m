function bool = isSemiImplicit(timeStepper)
    semiImplicitMethods = {'ab1be', 'ab2be', 'ab3cn', 'bdf1si', 'bdf2si', ...
        'bdf3si','bdf4si','bdf5si','bdf6si'};
    bool = ismember(char(timeStepper), semiImplicitMethods);
end
