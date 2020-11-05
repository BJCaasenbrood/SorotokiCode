<script src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script> 
<div align="center"> <img src="./bin/src/softrobot.png" width="650"> </div> <br/>

[**Sorotoki**](https://bjcaasenbrood.github.io/SorotokiCode/) is an open-source MATLAB toolkit that includes an array of modeling and designing tools for soft robotics. Due to its scientific diversity, exploring the field of soft robotics can be significantly challenging. With the aim to bridge this multi-disciplinary gap, Sorotoki aims to incorporate several layers of scientific disciplines in one toolkit. Examples include: continuum mechanics, dynamic system- and control theory, topology optimization, computer graphics, and much more! 

# Installation
Download the latest stable version ([.zip](https://github.com/BJCaasenbrood/SorotokiCode/zipball/master) or [.tar](https://github.com/BJCaasenbrood/SorotokiCode/tarball/master)), and unpack the compressed folder at any desired directory. Alternatively, you can directly clone the repository with the command:
```git
git clone https://github.com/BJCaasenbrood/SorotokiCode.git
```

To install the toolkit, simply run the command below. That's it, the soft robotics toolkit is now ready-to-use.
```ruby
% installation command
sorotoki();
```
The command above is also used to update the toolkit. It is recommended to run 'sorotoki.m' to check for updates occasionally. 

# Getting started
The Sorotoki toolkit includes a library of class-oriented functions for designing, modeling, and controlling a wide variety of soft robots. The toolkit is designed to solve particular sub-problem within soft robotics such that complexity can be easily decomposed in smaller sub-modules. For instance, the classes *Mesh()* and *Fem()* provide some numerical tools to model the continuum material behavior of soft structures, whose extrapolated material models can later be used in the dynamic modeling part with *Model()*. The main classes of Sorotoki are listed below:

```matlab
% list of classes
msh = Mesh();	    % meshing class
fem = Fem();   	    % finite elements class
mag = Mfem();       % magnetics finite elements class
obj = Gmodel();     % graphical model class
mdl = Model();      % dynamical model class
rig = Rig();        % IK rigging class
ctr = Control();    % (real-time) control class
```

To get started, type the following line in the command window:
```matlab
% show demos
sorotoki('demo');
```


# Documentation and libary usage
* [**Mesh generation**](./bin/Mesh.md). 
* [**Finite Element Method**](./bin/Fem.md).
* [**Topology Optimization**](./bin/Topo.md).
* [**Graphics and Rendering**](./bin/Gmodel.md).
* [**Dynamical Modeling**](./bin/Model.md).

## Citation
If you are using Sorotoki in your (academic) work, please consider citing the toolkit:
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

## Bibliography
[1] B. Caasenbrood, A. Pogromsky and H. Nijmeijer, *A Computational Design Framework for Pressure-driven Soft Robots through Nonlinear Topology Optimization,* 2020 3rd IEEE Inter. Conf. on Soft Robotics (RoboSoft), pp. 633-638, 2020.

[2] B. Caasenbrood, A. Pogromsky, and H. Nijmeijer, *Dynamic modeling of hyper-elastic soft robots using spatial curves,* IFAC World Congress, 2020.

[3] C. Talischi, G. H. Paulino, A. Pereira, and I. F. M. Menezes, *PolyMesher: A general-purpose mesh generator for polygonal elements written in Matlab*, Struct. Multidiscip. Optim., vol. 45, no. 3, pp. 309–328, 2012.

[4] N. Kim, *Introduction Analysis Finite Element to Nonlinear*. 2018.

[5] M. Bendsoe and O. Sigmund, *Topology Optimization. Theory, Methods and Applications*. 2003.


