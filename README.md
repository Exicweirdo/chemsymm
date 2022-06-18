# chemsymm
To decide point group symmetry of molecules

Support all point groups and SO(3) for single atom
## Example
```bash
(base) PS root-path> python .\src\main.py -f .\examples\molecules\BH3.gjf
 
path: .\examples\molecules\BH3.gjf
molecular formula: B1H3
coordinates:

        Atom    x               y               z
        B       2.1829e+00      -5.0000e-01     -1.7180e-01
        H       3.3629e+00      -5.0000e-01     -1.7180e-01
        H       1.5929e+00      5.2191e-01      -1.7180e-01
        H       1.5929e+00      -1.5219e+00     -1.7180e-01

symmetry information:
point group: D_3h
h: 12
rotation_type: symmetric
inertia moment: [2.08859999 2.0886     4.17719999]
i: False
standard_base: [[-8.52263495e-01  5.23112736e-01  0.00000000e+00]
 [-5.23112736e-01 -8.52263495e-01 -2.71588768e-23]
 [-1.42071544e-23 -2.31465193e-23  1.00000000e+00]]
---------------------------------------------------
```
