<script src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script> 
<div align="center"> <img src="./bin/src/softrobot.png" width="650"> </div> <br/>

Sorotoki is an open-source MATLAB toolkit that includes an array of modeling and designing tools for soft robotics. The toolkit aims to bring together the various scientific disciplines within the field of soft robotics, e.g., continuum mechanics, dynamical system- and control theory, and computer graphics. 

# Installation
Download the latest stable version ([.zip](https://github.com/BJCaasenbrood/SorotokiCode/zipball/master) or [.tar](https://github.com/BJCaasenbrood/SorotokiCode/tarball/master)) and unpack the compressed folder at any desired work directory. To install the toolkit, simply run the command below. That's it, the soft robotics toolkit is now ready-to-use.
```matlab
% installation command
sorotoki();
```

# Getting started
The main class of **ToPy** is 'Topology'. It defines the main constraints, grid and parameters of optimization -- but you don't really have to bother yourself with this if you just want to get some results.
```matlab
% list of classes
msh = Mesh();	 % meshing class
fem = Fem();   	 % finite elements class
obj = Gmodel();  % graphics model class
mdl = Model();   % dynamical model class
```

### Documentation and libary usage
* [**Mesh generation**](./bin/Mesh.md). 
* [**Nonlinear Finite Elements**](./bin/Fem.md).
* [**Computational Design**](./bin/Topo.md).
* [**Graphics and Implicit Modeling**](./bin/Gmodel.md).
* [**Dynamical Modeling**](./bin/Model.md).

# Citation
If you are using Sorotoki in your academic work, please consider to cite the toolkit:
```
@misc{Caasenbrood2018,
  author = {Caasenbrood, Brandon},
  title = {Sorotoki - A Soft Robotics Toolkit for MATLAB},
  year = {2018},
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\url{https://github.com/BJCaasenbrood/SorotokiCode}},
  }
```

# References
[1] B. Caasenbrood, A. Pogromsky, and H. Nijmeijer, *A Computational Design Framework for Pressure-driven Soft Robots through Nonlinear Topology Optimization*, RoboSoft 2019 - 2019 IEEE International Conference on Soft Robotics, 2020.

[2] N. Kim, *Introduction Analysis Finite Element to Nonlinear*. 2018.

[3] M. Bendsoe and O. Sigmund, *Topology Optimization. Theory, Methods and Applications*. 2003.


