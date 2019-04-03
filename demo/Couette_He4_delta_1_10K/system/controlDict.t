/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  1.6                                   |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      controlDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

application     dsmcFoamPlus;

startFrom       latestTime;

nTerminalOutputs    100;

startTime       0;

stopAt          endTime;

endTime         9.812263e-05;

deltaT          9.812263e-09;

writeControl    runTime;

writeInterval   9.812263e-06;

purgeWrite      1;

writeFormat     ascii;

writePrecision  6;

writeCompression off;

timeFormat      general;

timePrecision   6;

runTimeModifiable true;

adjustTimeStep  no;

functions
(
   
);

libs
(
   "libAIdsmcCollision.so"
);

// ************************************************************************* //
