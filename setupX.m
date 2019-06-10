function x = setupX(xLength, yLength, xN, yN)
    xL = [xLength, yLength];
    xN = [xN, yN];
    xS = xL ./ xN;
    dim = 2;
    x = cell(dim, 1);
    for n = 1:dim
        x{n} = linspace(xS(n), xL(n), xN(n))';
    end
end