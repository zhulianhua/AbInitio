## Quantum interatomic scattering in dsmcFoamPlus

This code 


_Ab initio_ potential of interatomic interaction


### How to comiler it?

NOTE: Require to install the [OpenFOAM-2.4.0-MNF](https://github.com/MicroNanoFlows/OpenFOAM-2.4.0-MNF) first.

This code should be compilerd to an user library by:

```bash
wmake
```

#### How to use it in your dsmc run case?

First add the following entry to the `system/controlDict`:

```txt
libs
(
   "libAIdsmcCollision.so"
);
```

Then in constant/



``` txt
// Binary Collision Model
// ~~~~~~~~~~~~~~~~~~~~~~

BinaryCollisionModel            AbInitio;


AbInitioCoeffs
{
    deflectionAngleTableFileName    "mmc2.csv";
    numRows                         800;
    numColumns                      100;

}
```




### Reference

