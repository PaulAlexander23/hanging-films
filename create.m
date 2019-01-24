function create(delta,theta,Re,We,C,xL,yL,T)
        
    % delta, theta, Re, We, C
    params = [delta,theta,Re,We,C];
    
    xN = 2^7;
    xS = xL/xN;
    x = linspace(xS,xL,xN)';
    
    yN = 2^7;
    yS = yL/yN;
    y = linspace(yS,yL,yN);
    
    V = zeros(1,1,25,2);
    V(1,1,:,:) = combvec(1:5,1:5)';
    
    R1 = rand(1,1,25);
    R2 = rand(1,1,25);
    R3 = rand(1,1,25);
    R4 = rand(1,1,25);
    
    H0 = 1 ...
        + sum(1e-4 * R1 .* cos(V(1,1,:,1) .* x + V(1,1,:,2) .* y + 2*pi*R2),3)...
        + sum(1e-4 * R3 .* cos(V(1,1,:,1) .* x - V(1,1,:,2) .* y + 2*pi*R4),3);

    tic
    [H, ~, t] = solver(@f_benney, params, H0, T, [xL, yL], [xN, yN], @(~,H) is_wetted(H), 1e-3);
    toc
    
    filename = replace(sprintf('data-d-%g-theta-%g-Re-%g-We-%g-C-%g-xL-%g-yL-%g-T-%g.mat',[params,xL,yL,T]),'.','_');
    save(filename,'H','params','t','x','y');
end
