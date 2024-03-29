![teaser](https://user-images.githubusercontent.com/38473042/167678103-ce292271-4962-4726-98c0-ea6cf376eb98.png)

# QGP3D – Quantized Global Parametrization of Volumes

`QGP3D` is an implementation of [Volume Parametrization Quantization for Hexahedral Meshing \[Brückler et al. 2022\]](http://graphics.cs.uos.de/papers/Volume_Parametrization_Quantization-SIGGRAPH2022.pdf) (SIGGRAPH 2022) distributed under GPLv3.

If you make use of `QGP3D` in your scientific work, please cite our paper. For your convenience,
you can use the following bibtex snippet:

    @article{QGP3D,
        author     = {Hendrik Br{\"{u}}ckler and
                     David Bommes and
                     Marcel Campen},
        title      = {Volume Parametrization Quantization for Hexahedral Meshing},
        journal    = {ACM Trans. Graph.},
        volume     = {41},
        number     = {4},
        year       = {2022},
    }

***

## How does it work?

`QGP3D` makes use of the [3D Motorcycle Complex](https://github.com/HendrikBrueckler/MC3D) to partition a tetrahedral mesh, equipped with a suitable seamless map, into blocks.
It then constructs a valid quantization of the non-conforming partition, that is an assignment of integer lengths to the partition's edges, and deduces optimal integer spacings
between all crucial entities like singularities and boundaries from this.
We guarantee that these spacings, when used as constraints in reparametrization, permit the construction of a globally valid integer-grid-map parametrization
from which a hex mesh without defects can be extracted. Greedy rounding, the previous state-of-the-art approach for quantization, can not provide this guarantee and
therefore often enforces defects in the output.

![QGP3D](https://user-images.githubusercontent.com/38473042/167678775-2de020fd-1527-4e43-9350-fab9e6482ef2.png)


***

### Dependencies
- GMP (NOT included, must be installed on your system)
- GUROBI (NOT included, must be installed on your system), set GUROBI_BASE via cmake to point to correct folder
    * OPEN-SOURCE-ALTERNATIVE: coinor-bonmin, including all its dependencies, namely coinutils, Ipopt (+HSL MA27/coinhsl), Osi, Cgl, Cbc, Clp (NOT included, must be installed on your system)
    * This is automatically chosen if GUROBI could not be detected on your system
- [MC3D](https://github.com/HendrikBrueckler/MC3D) (Included as submodule, together with all subdependencies)

### Building
In root directory

    mkdir build
    cd build
    cmake -DGUROBI_BASE=<path/to/gurobi/> ..
    make

### Usage
An example command-line application is included that reads a tetrahedral mesh including a seamless parametrization from a file in .hexex-format, as used and documented in [libHexEx](https://www.graphics.rwth-aachen.de/software/libHexEx/).
It outputs the integer spacing constraints between critical vertices of the input mesh.

After building the CLI app can be found in ```build/Build/bin/cli``` .
For full information on its usage, execute

    qgp3d_cli --help

Example input can be found in folder ```extern/MC3D/tests/resources```.

### API
A simple interface is provided in ```QGP3D/Quantizer.hpp```.
It can be used like this:

```cpp
// Generation of tet mesh and seamless parametrization, optionally feature markers
// ...

// tetMesh: your mc3d::TetMesh
qgp3d::Quantizer quantizer(tetMesh);

// param: your property defining mc3d::Vec3d parameter values for each tet-vertex-pair
for (auto tet: tetMesh.cells())
    for (auto v: tetMesh.cell_vertices(tet))
        quantizer.setParam(tet, v, param[tet][v]);

// isFeatureF, isFeatureE, isFeatureV: your properties marking features
for (auto f: tetMesh.faces())
    if (isFeatureF[f])
        quantizer.setFeature(f, true);
for (auto e: tetMesh.edges())
    if (isFeatureE[e])
        quantizer.setFeature(e, true);
for (auto v: tetMesh.vertices())
    if (isFeatureV[v])
        quantizer.setFeature(v, true);

// This scaling factor is applied to the parametrization before quantization
double scaling = 1.0;
// This will contain integer spacing constraints that will enforce a valid quantization
std::vector<qgp3d::PathConstraint> constraints;
// This will hold the number of hexahedra implied by the quantization
int nHexahedra = 0;

// Compute a quantization and retrieve the corresponding spacing constraints
quantizer.quantize(scaling, constraints, nHexahedra);

// Generate a quantized seamless parametrization (integer-grid map)
// by enforcing the retreived integer spacing constraints during reparametrization
// ...

```

For more details on the API of the library, check the headers in ```include```, they are thoroughly documented. Apart from that, ```cli/main.cpp``` demonstrates usage of the entire pipeline for both simple and advanced usage.
