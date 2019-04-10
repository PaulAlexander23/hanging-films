function y0 = irand(x)
    %IRAND Interface on x with 5 random modes
    
    A = 0.1;
    V = zeros(1,1,25,2);
    V(1,1,:,:) = combvec(1:5,1:5)';
    R1 = rand(1,1,25);
    R2 = rand(1,1,25);
    R3 = rand(1,1,25);
    R4 = rand(1,1,25);
    y0 = 1 ...
        + sum(A * R1 .* cos(V(1,1,:,1) .* x{1} + V(1,1,:,2) .* x{2}' + 2*pi*R2),3)...
        + sum(A * R3 .* cos(V(1,1,:,1) .* x{1} - V(1,1,:,2) .* x{2}' + 2*pi*R4),3);
    
end