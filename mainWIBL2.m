main("wibl2", 7*pi/8, 1, 0.01, 32, 32, 1, @icosWIBL2, 32, 32, 1e-3, "finite-difference", @ode15s, 0.1, 1e-2, "00:10:00")
%main("wibl2", 7*pi/8, 1, 0.01, 32, 32, 1, @icosWIBL2, 64, 64, 1e-3, "pseudo-spectral", @ode45, 0.1, 1e-2, "00:10:00")
