# Liggghts_VBP_VTKReader

A simple tool that parses VTK files generated by LIGGGHTS for Vibrating Packed Bed simulations and calculates particle counts in slices of the packed bed.

## How to build

You need to install vtk to build this project. It has been tested using vtk6.3.
This also depends on matplotlibcpp to plot graphs. A [fork](https://github.com/Cryoris/matplotlib-cpp) of the [original matplotlibcpp](https://github.com/lava/matplotlib-cpp) has been used. Therefore you also need to install matplotlib and the python2-dev libs.

```console
sudo apt-get install libvtk6.3 libvtk6-dev
sudo apt-get install python-matplotlib python2.7-dev
```

From root directory of the project:
```console
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```


## How to run

From build directory:
```console
./Liggghts_VBP_VTKReader <path to vtk file>
```
