<script src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script> 
<div align="center"> <img src="./bin/src/softrobot.png" width="650"> </div>

Sorotoki is an open-source toolkit for developing, modeling and controlling soft robots. 

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

# Finite Elements for Hyper-elastic Materials

$$\begin{aligned} 
\Psi & = c_1 (J_1 - 3) & (\text{Neo-Hookean})  \\ \vspace{2mm}
\Psi & = c_1 (J_1 - 3) + c_2 (J_2 - 3) & (\text{Mooney-Rivlin}) \\
\Psi & = \sum^3_{i=1} c_i (J_1 - 3)^i & (\text{Yeoh}) 
\end{aligned}
$$ 

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


