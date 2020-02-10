<div align="center"> <img src="./src/gmodel.png" width="650"> </div>

# Graphics and Implicit Modeling
[Go back to home page](https://bjcaasenbrood.github.io/SorotokiCode/)

## Rendering
Sorotoki is equipped with a wide range of material rendering options. 

### Material captures
```matlab
%% preview
obj = Gmodel(@(x) dSphere(x,0,0,0,1),[-1,1,-1,1,-1,1]);

%% material list
mat = {aniso, bluered, bump, chroma, chromium, ...
       clean, copper, egg, grey, hotmetal, jade, ...
       matcap, mateplastic, metal, metalclean, ...
       oldwax, orangeresin, planet, plastic, ...
       redshine, redwax, retro, rim, soft, ...
       skin, studioclay};

% loop materials
for ii = 1:length(mat)
  pause(0.1);
  obj.Texture = mat{ii};
  if ii == 1, obj.render;
  else, obj.update();
  end
end
```

<div align="center"> <img src="./src/matcap.gif" width="250"> </div>

### Ambient occlusion (AO)
<div align="center"> <img src="./src/bunny_AO.png" width="450"> </div>


### Sub-Surface Scattering (SSS)
<div align="center"> <img src="./src/bunny_AO.png" width="450"> </div>
