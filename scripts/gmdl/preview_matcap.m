clr;
%% preview
obj = Gmodel(@(x) dSphere(x,0,0,0,1),[-1,1,-1,1,-1,1]);

mat = {aniso, bluered, bump, chroma, chromium, ...
       clean, copper, egg, grey, hotmetal, jade, ...
       matcap, mateplastic, metal, metalclean, ...
       oldwax, prusa, planet, plastic, ...
       redshine, redwax, retro, rim, soft, ...
       skin, studioclay};
   
for ii = 1:length(mat)
pause(0.1);
obj.Texture = mat{ii};
if ii == 1, obj.render;
else, obj.update();
end
end
