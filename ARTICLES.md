# Reference Articles

## RTree Basics

- [blog reference](https://geoffboeing.com/2016/10/r-tree-spatial-index-python/)
- Datasets API: [`osmnx`](https://osmnx.readthedocs.io/en/stable/index.html) [1]
- rtree implementation from [`rtree 0.7.0`](https://toblerity.org/rtree/index.html#)

## Visibility operations

- shapely
- geopandas
- scikit-geometry (skgeom)
    [blog reference](https://wolfv.medium.com/introducing-scikit-geometry-ae1dccaad5fd)
- [Rotational sweep visibility](https://doc.cgal.org/latest/Visibility_2/classCGAL_1_1Rotational__sweep__visibility__2.html): Theoretical background for implementation

    ### skgeom Installation
    1. Install CGAL
    ```bash
    sudo apt-get install libcgal-dev
    ```
    2. Download the `scikit-geometry` git repo
    3. Inside the repo, run
    ```bash
    sudo python3 setup.py install
    ```

    ### CGAL Precondition error

    [NOTE]: `conda` does provide a stable python package for scikit-geometry which can be installed directly. For more info, kindly check their docs.
    This is an errors arising only with the newer versions of CGAL. Follow these steps to fix it.
    1. Comment that precondition line.
    2. Uninstall the `skgeom` python package and reinstall after cleaning the installation cache folders.
    (In case you have done a manual install after cloning the repo, delete the newly generated folders.)

## Papers:

[1] Boeing, G. 2017. OSMnx: New Methods for Acquiring, Constructing, Analyzing, and Visualizing Complex Street Networks. Computers, Environment and Urban Systems 65, 126-139. **[doi link](doi:10.1016/j.compenvurbsys.2017.05.004)**