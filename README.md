## Quantum interatomic scattering in dsmcFoamPlus

This small code implements the _ab initio_ potentials based binary collision model in dsmcFoamPlus solver. 

Currently works only for <sup>3</sup>He, <sup>4</sup>He and Ne single-component gas.


### How to compile it?

NOTE: Require to install the [OpenFOAM-2.4.0-MNF](https://github.com/MicroNanoFlows/OpenFOAM-2.4.0-MNF) first.

This code should be compilerd to an user library by:

```bash
wmake
```

### How to use this collision model in your dsmcFoamPlus run case?

####  Step 1: change system/controlDict
 First add the following entry to the `system/controlDict` to reference to the user library compiled just now.

```txt
libs
(
   "libAIdsmcCollision.so"
);
```

####  Step 2: change to this binary collision model

Then in the constant/dsmcProperties


Use the following binary collision model specification:

``` txt
// Binary Collision Model
// ~~~~~~~~~~~~~~~~~~~~~~

BinaryCollisionModel            AbInitio;

AbInitioCoeffs
{
    He-He
    {
        deflectionAngleCosinTableFileName    "xiHe4.csv";
        numRows                         900;
        numColumns                      100;
        G                               400.0;
    }
    He-Ne
    {
        deflectionAngleCosinTableFileName    "xiHe4Ne.csv";
        numRows                         800;
        numColumns                      100;
        G                               400.0;
    }
    Ne-Ne
    {
        deflectionAngleCosinTableFileName    "xiNe.csv";
        numRows                         900;
        numColumns                      100;
        G                               200.0;
    }
}

```

The `deflectionAngleTableFileName` specifies the file stores the deflection angle's cosine values maxtrix (table) and the total cross section area (unit: 10<sup>-20</sup>). This file has `numColumns+2` columns and `numRows` rows of data. The meaning of data are explained in Ref.[1]. The coressponding files for <sup>3</sup>He, <sup>4</sup>He  has been given in the `deflectionAngleTables` directories, they are explained as follows:

* `'mmc1.csv'` for the <sup>3</sup>He that is used and give in Ref.[1], and is valid up to temperature of 3000K and down to 1K.  `nowRows`, `nowColumns` and `G` in the `AbInitioCoeffs` dictionary should be 800, 100 and 400, respectively.

* `'mmc2.csv'` for the <sup>4</sup>He that is used and give in Ref.[1], and valida up to temperature of 3000K and down to 1K. `nowRows`, `nowColumns` and `G` in the `AbInitioCoeffs` dictionary should be 800, 100 and 400, respectively.

* `'xiHe3'` for the <sup>3</sup>He and is valid up to temperature of 10000K (higher than `'mmc1.csv'`) and down to 1K. `nowRows`, `nowColumns` and `G` in the `AbInitioCoeffs` dictionary should be 800, 100 and 200, respectively.

* `'xiHe4'` for the <sup>4</sup>He and is valid up to temperature of 10000K (higher than `'mmc1.csv'`) and down to 1K. Along this file, `nowRows`, `nowColumns` and `G` in the `AbInitioCoeffs` dictionary should be 900, 100 and 400, respectively.

* `'xiNe'` for the Ne gas. Along this file, `nowRows`, `nowColumns` and `G` in the `AbInitioCoeffs` dictionary should be 900, 100 and 200,  respectively.

**NOTE**: This collision model will use only the `mass` entry in the `moleculeProperties` dictionary of `constant/dsmcProperties`.

####  Step 4: Run the case

Keep other setting the same as using other VHS/VSS models, and run the dsmcFoamPlus solver.

### Demostrational case

A demostration case is provided in the `demo` directory to reproduce the resuts in Ref. [1] using the `mmc2.csv` file.


### Limitations


### Reference

1. Modeling of transport phenomena in gases based on quantum scattering, _Physica A: Statistical Mechanics and its Applications
Volume 508, Pages 797-805_ https://doi.org/10.1016/j.physa.2018.05.129

### Acknowledgement

Professor Felix Sharipov [http://fisica.ufpr.br/sharipov/].
