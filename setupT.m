function t = setupT(tFinal, tStep)
    t = 0:tStep:tFinal;
    if rem(tFinal, tStep) ~= 0
        t = [t, tFinal];
    end
end