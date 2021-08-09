---
layout: default
title: Example 1 - Rendering the Stanford Bunny (.stl)
parent:  Computer Graphics
grand_parent: Examples
nav_order: 1
permalink: /docs/examples/3_graphics/renderstl/
---

#  Rendering the Stanford Bunny

<div align="center"> <img src="./img/rotateBunny.gif" width="450"> </div>
<div align="center"> Rendering the Stanford bunny in MATLAB with responsive textures [1].  </div>

---

#### Difficulty: `beginner`{: .fs-3 .text-green-200}
 - Required classes: `Gmodel`{: .text-purple-000}
 - Lines of code: `4 lines`{: .text-purple-000} (without comments)

---

Rendering a 
```matlab
obj = Gmodel('Bunny.stl');  
obj.Material = base;  
```

``` matlab
obj = obj.bake();
obj = obj.render();
```
Alternatively, we can rewrite the code above more compactly.
``` matlab
obj = obj.bake().render();
```

```matlab
view(20,10);
obj.update();
```

## Complete code (4 lines without comments)
```matlab
%% loading .stl file
obj = Gmodel('Bunny.stl');
     
%% asssigning texture 
obj.Texture = base;

%% baking texture and rendering the bunny
obj = obj.bake().render(); view(10,20);
obj.update();
```

## References 
[**[1]**](http://graphics.stanford.edu/data/3Dscanrep/) **The "Stanford Bunny"**, Stanford University Computer Graphics Laboratory -- (4.9 MB compressed, 22 MB uncompressed) [http://graphics.stanford.edu/data/3Dscanrep/](http://graphics.stanford.edu/data/3Dscanrep/)
{: .fs-3} 
