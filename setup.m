%SETUP 

xN = 2^7;
xL = 10;
xS = xL/xN;
x = linspace(xS,xL,xN)';

yN = 2^7;
yL = 10;
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

% delta, theta, Re, We, C
params = [1,pi/2,1,0,1];

tFinal = 1;

%plot_surface(x,y,H0)
