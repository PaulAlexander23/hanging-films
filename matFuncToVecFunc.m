function F = matFuncToVecFunc(f)

    F = @(domain, y, params) domain.reshapeToVector( ...
        f(domain, domain.reshapeToDomain(y), params));
end
