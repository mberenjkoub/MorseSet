# MorseSet

The goal of this project was to extract the Morse Sets defined as strongly connected component regions in the flow. For the first time, we propose a method to extract the Morse Sets in piece-wise linear vector field. The input of this software is the velocity fields data for 3D datasets and the output is the extracted Morse sets. 
This software implemented in C++ (windows) using two libraries of Qt and CUDA and OpenGL. 


## Description
This software has three features:
- Visualizing the streamlines.
- Extract the Morse sets in different levels.
- Visualizaing the Morse sets.

## Binding to Other libraries
I used 
- [CUDA7.5](https://developer.nvidia.com/cuda-75-downloads-archive)
- [Qt](https://www.qt.io/)



## Further Reading
[Morse Decomposition of 3D Piecewise Linear Vector Fields](http://www.ingentaconnect.com/content/ist/ei/2016/00002016/00000001/art00003?crawler=true&mimetype=application/pdf)

## install
In order to run this code, you need to 
- Install Qt5.8 library.
- Install CUDA7.5 library.
- Install OpenGL 2.1

and you can use the software.
