# Hanging films project

The aim is to generate numerical solutions to the 3D Benney, WIBL1 and WIBL2 models of fluid flow down the underside of an inclined plate.
The programming language used is Matlab, tested for R2018b.
A collaboration with Dr. Ruben Tomlin and Prof. Demetrios Papageorgiou.

Run tests in Matlab: 

    runtests

Code coverage:

    profile on
    runtests
    profile off

Current Folder > Reports > Coverage Report

## Usage:

Include a ODE function in the folder with an interface:

    dydt = fODE(domain, y, params);

where 

    domain

is an object formed by initialising one of the domain classes within the discretisation methods folder

    y

is a vector of interface values, and

    params

is a structure of all the parameters contained within the ODE.


The solver is then called as,

    main(fODE, theta, Re, C, xLength, yLength, tFinal, interface, xN, yN, RelTol, method, timeStepper, timeStepOut, timeStep, timeout)
