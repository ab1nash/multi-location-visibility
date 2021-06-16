# Reference Articles

## RTree Basics

- https://geoffboeing.com/2016/10/r-tree-spatial-index-python/


## Visibility operations

- shapely
- geopandas

- scikit-geometry

    ### Installation
    1. Install CGAL
    ```bash
    sudo apt-get install libcgal-dev
    ```
    2. Download the `scikit-geometry` git repo
    3. Inside the repo, run
    ```bash
    sudo python3 setup.py install
    ```

    ### CGAL Precindition error
    1. Comment that precondition line.
    2. Uninstall the `skgeom` python package and reinstall after cleaning the installation cache folders.
    (In case you have done a manual install after cloning the repo, delete the newly generated folders.)