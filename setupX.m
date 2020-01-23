function x = setupX(xLength, yLength, xN, yN)
    xL = [xLength, yLength];
    xN = [xN, yN];
    xS = xL ./ xN;
    
    x{1} = linspace(xS(1), xL(1), xN(1))';
    x{2} = linspace(xS(2), xL(2), xN(2));
end