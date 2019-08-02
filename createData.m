function createData(model, theta, Re, C, xLength, yLength, tFinal, interface, xN, yN, AbsTol, method)
    if nargin < 9, xN = 64; end
    if nargin < 10, yN = 64; end
    if nargin < 11, AbsTol = 1e-6; end
    if nargin < 12, method = "finite-difference"; end
    
    if model == "benney"
        create(theta, Re, C, xLength, yLength, tFinal, interface, xN, yN, AbsTol);
    elseif model == "wibl1"
        if method == "finite-difference"
            createWIBL1(theta, Re, C, xLength, yLength, tFinal, interface, xN, yN, AbsTol)
        elseif method == "pseudo-spectral"
            createWIBL1PseudoSpectral(theta, Re, C, xLength, yLength, tFinal, interface, xN, yN, AbsTol)
        end
    end
end