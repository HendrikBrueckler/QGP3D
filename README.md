![teaser](https://user-images.githubusercontent.com/38473042/167607191-79ecfd18-2e50-4a88-a0f3-76cd31808d39.png)

# quant3d – Quantization of Volumetric Parametrizations

`quant3d` is an implementation of [Volume Parametrization Quantization for Hexahedral Meshing \[Brückler et al. 2022\]](http://graphics.cs.uos.de/papers/Volume_Parametrization_Quantization-SIGGRAPH2022.pdf) (SIGGRAPH 2022) distributed under GPLv3.

If you make use of `quant3d` in your scientific work, please cite our paper. For your convenience,
you can use the following bibtex snippet:

    @article{quant3d,
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

`quant3d` makes use of the [3D Motorcycle Complex](https://github.com/HendrikBrueckler/MC3D) to partition a tetrahedral mesh, equipped with a suitable seamless map, into blocks. 
It then constructs a valid quantization of the non-conforming partition, that is an assignment of integer lengths to the partition's edges, and deduces optimal integer spacings 
between all crucial entities like singularities and boundaries from this.
We guarantee that these spacings, when used as constraints in reparametrization, permit the construction of a globally valid integer-grid-map parametrization 
from which a hex mesh without defects can be extracted. Greedy rounding, the previous state-of-the-art approach for quantization, can not provide this guarantee and 
therefore often enforces defects in the output.

![quant3d](https://user-images.githubusercontent.com/38473042/167615196-3d012d93-f403-48cb-ae16-40e9b7a4fe45.png)


***

### Dependencies
- GMP (NOT included, must be installed on your system)
- GUROBI (NOT included, must be installed on your system)
- [MC3D](https://github.com/HendrikBrueckler/MC3D) (Included as submodule, together with all subdependencies)
- [TrulySeamless3D](https://github.com/HendrikBrueckler/TrulySeamless3D) (Included as submodule, together with subdependency libHexEx)

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

    quant3d_cli --help

Example input can be found in folder ```extern/MC3D/tests/resources```.

### API
For details on the API of the library, check the headers in ```include```, they are thoroughly documented. Apart from that, ```cli/main.cpp``` demonstrates usage of the entire pipeline for both simple and advanced usage.
