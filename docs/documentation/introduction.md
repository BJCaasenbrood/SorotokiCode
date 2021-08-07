---
layout: default
title: About SOROTOKI
parent: Documentation
nav_order: 1
---

# About

SOROTOKI is an open-source MATLAB toolkit for soft robotics that includes an array of tools for design, modeling, and control. Due to its scientific diversity, it can be challenging for researchers to quickly familiarize themselves with multiple scientific disciplines. With the aim to lower the threshold, Sorotoki aims to incorporate multiple layers of soft robotics research into one toolkit. Examples include: continuum mechanics, dynamic systems and control theory, topology optimization, computer graphics, and much more to come! The combination provides a highly flexible modeling environment and hopefully aids the development of novel soft robotic research.

# Object-oriented architecture
bla
```mermaid!
classDiagram
class Sdf {
	<<signed distance function>>
	eval()
  	show()
}

class Mesh {
	<<mesh generation>>
	Blender
  	generate()
}

class Gmodel {
	<<graphical models>>
	Material
  	render()
  	bake()
  	update()
}

class Fem {
	<<finite elements>>
	Material
  	solve()
  	optimize()
}

class Model {
	<<dynamic beam models>>
	Rig
	Shapes
  	build()
  	simulate()
}


class Bdog {
	<<Real-time interface via Balloondog>>
	SetPoint
  	connect()
 	execute()
  	shell()
  	close()
}

Sdf --> Mesh : Meshing domain
Sdf ..> Gmodel : Implicit modeling
Mesh ..> Gmodel
Gmodel ..> Model : Rendering

Mesh --> Fem : Import
Fem --> Mesh  : Export topology
Fem ..> Model : FEM-driven data
Model ..> Bdog : Model-based control
Fem ..> Bdog : FEM-based IK-control

``` 
