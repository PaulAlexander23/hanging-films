function tests = testInterp()
    tests = functiontests(localfunctions);
end

function testPeriodicInterp2SplineVector(testCase)
    x = linspace(0.5,48,96)';
    y = linspace(0.5,32,64);
    xq = linspace(0.25,48,192)';
    yq = linspace(0.25,32,128);
    V = cos(pi/24 * x) + sin(pi/16*y);

    expected = cos(pi/24 * xq) + sin(pi/16*yq);
    actual = periodicInterp2(x,y,V,xq,yq,'spline');

    verifyEqual(testCase, actual, expected, 'AbsTol', 1e-3);
end

function testPeriodicInterp2SplineMeshgrid(testCase)
    [X, Y] = meshgrid(linspace(0.5,48,96),linspace(0.5,32,64));
    [Xq, Yq] = meshgrid(linspace(0.25,48,192),linspace(0.25,32,128));
    V = cos(pi/24 * X) + sin(pi/16*Y);

    expected = cos(pi/24 * Xq) + sin(pi/16*Yq);
    actual = periodicInterp2(X,Y,V,Xq,Yq,'spline');

    verifyEqual(testCase, actual, expected, 'AbsTol', 1e-3);
end
