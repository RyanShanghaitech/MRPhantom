# MRI Digital Phantom

## Introduction
A quantitative phantom of a slime with cardiac and respiratory motion.

Since a 3D + T dataset is so large, this package generates phantom on each tick on demand, according to the given shape parameters. This phantom has a parallel C++ backend for nearly real-time 3D phantom generation.

## Install
Install `gxx` via conda:
```bash
$ conda install gxx 
```
Using `pip`:
```bash
$ pip install mrphantom
```
Offline:
```bash
$ bash install.bash
```

## Usage
Please refer to the [Examples](https://github.com/RyanShanghaitech/MRPhantom/tree/main/example) page to learn how to use this package. You can find examples for 2D dynamic, 3D dynamic, static simulations, and also the method to obtain the quantitative parameters (e.g. T1, T2, off-resonance) of each pixel.