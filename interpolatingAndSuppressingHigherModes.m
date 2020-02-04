

x = {linspace(32/64, 32, 64)', linspace(32/64, 32, 64)};
domain = FDDomain(x,[1,0;2,0;3,0;4,0;0,1;0,2;0,3;0,4]',2);
%%
fy = fft2(y);
figure; surf(log10(abs(fy)));
index = abs(fy)>1e-1;
fy = fy .* index;
figure;surf(log10(abs((fy))))
y2 = ifft2(fy);
figure;surf(y2)
laph = domain.diff(y2,[2,0]')+domain.diff(y2,[0,2]');
figure; surf(laph);

%%
fF = fft2(F);
figure;surf(log10(abs(fF)));
index = abs(fF)>10^(0.5);
fF = fF .* index;
figure;surf(log10(abs((fF))))
F2 = ifft2(fF);
figure;surf(F2)
lapF = domain.diff(F2,[2,0]')+domain.diff(F2,[0,2]');
figure; surf(lapF);