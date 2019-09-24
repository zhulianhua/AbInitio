## Quantum interatomic scattering in dsmcFoamPlus

### Project Homepage
[AbInitio](https://github.com/zhulianhua/AbInitio/)

### What it is for?

This project implements the _ab initio_ potential based collision model in the dsmcFoamPlus and dsmcFoam (in the future ?) solver. 

The _ab initio_ potential DSMC with quantum scattering is introduced in Ref. [1, 2] etc. 


### How to compile it in OpenFOAM-2.4.0-MNF?

NOTE: Require to install the [OpenFOAM-2.4.0-MNF](https://github.com/MicroNanoFlows/OpenFOAM-2.4.0-MNF) first.

This code should be compilerd to a user library by typein

```bash
wmake
```
in the project directory.

### How to compile it in latest offical OpenFOAM releases?
[TO BE ADDED]

### How to use it?

####  Step 1: change system/controlDict
 Firstly, add the following user defined library entry to the `system/controlDict`, so the dsmcFoam solver can dynamically load the library compiled above.

```txt
libs ("libAIdsmcCollision.so");
```

####  Step 2: change to this binary collision model

Then in the `constant/dsmcProperties` file use the binary-collision model like follows:

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
        deflectionAngleCosinTableFileName    "xiHe4-Ne.csv";
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

The `deflectionAngleTableFileName` option specifies the file stores the deflection angles' cosine values and the total cross section area for a serial of discrete relative velocitie $g_i$. These files have **`numColumns+2` columns** and `numRows` rows of floating point number. Each row of the files stores:
*  the relative velocity at the last elsemnt;
*  the total cross-section (TCS, $\sigma_T$) in the unit of of $10^{-20} m^2$ at the last-but-one element;
*  the **numColumns**  equally probabel deflection angles' cosin values in the first numColumns elements;

The files used should be provided in the `constant` directory of the dsmcFoam case. The fllowing files are gathered from various literature have been given in the `deflectionAngleTables` directory of this repository and their specification are explained in the following table.

| matrix file | numRows | numColumns | G | comment |
|-------------|-------|---------|---|---------|
|  `xiHe3.csv` | 900 | 100 | 400.0  |   Quantum potential of <sup>3</sup>He-<sup>3</sup>He collision, valid up to tempeature of 15,000K, provided in Ref. [2] |
|  `xiHe4.csv` | 900 | 100 | 400.0  |   Quantum potential of <sup>4</sup>He-<sup>4</sup>He collision, valid up to tempeature of 15,000K, provided in Ref. [2] |
|  `xiNe.csv`  | 900 | 100 | 200.0  |   Quantum potential of Ne-Ne collision, valid up to tempeature of 15,000K, , provided in Ref. [2]|
|  `xiHe4-Ne.csv` | 800 | 100 | 400.0  |   Quantum potential of <sup>4</sup>He-Ne collision, valid up to tempeature of 15,000K, provided by Professor Felix Sharipov |

**NOTE**:
* This collision model will use only the `mass` entry in the `moleculeProperties` dictionary of `constant/dsmcProperties`.
* These files can be used universally for different case setups, i.e, no need to change.

####  Step 4: Run the case

Keep other setting the same as using other VHS/VSS models, and run the dsmcFoamPlus solver.

### Demostrational case

A demostration case is provided in the `demo` directory to reproduce the resuts in Ref. [1] using the `mmc2.csv` file.

### Reference

1. Sharipov, F., 2018. Modeling of transport phenomena in gases based on quantum scattering. Physica A: Statistical Mechanics and its Applications 508, 797–805. https://doi.org/10.1016/j.physa.2018.05.129

2. Sharipov, F., Dias, F.C., 2019. Temperature dependence of shock wave structure in helium and neon. Physics of Fluids 31, 037109. https://doi.org/10.1063/1.5088556


### Acknowledgement

The discusstions with Professor Felix Sharipov [http://fisica.ufpr.br/sharipov/] contributed to this project.
