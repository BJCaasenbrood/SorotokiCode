<div align="center">
	<img src="./bin/src/softrobot.png" width="700">
</div>


Sorotoki is an open-source toolkit for Soft Robotics. 

# Installation
Download the latest stable version ([.zip](https://github.com/BJCaasenbrood/SorotokiCode/zipball/master) or [.tar](https://github.com/BJCaasenbrood/SorotokiCode/tarball/master)) and unpack the compressed folder at your desired work directory. To install, simply execute the command below, and that's it, the toolkit is ready to use.
```matlab
% installation command
sorotoki
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

# References
[1] B. Caasenbrood, A. Pogromsky, and H. Nijmeijer. *A Computational Design Framework for Pressure-driven Soft Robots through Nonlinear Topology Optimization*, RoboSoft 2019 - 2019 IEEE International Conference on Soft Robotics, 2020.

# Citation
If you are using SOROTOKI in your academic work, please consider to cite the toolkit:
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

