clr; 
%% preview
obj = Gmodel(sSphere(1),'Quality',50);

mat = {aniso, bluered, bump, chroma, chromium, ...
       clean, copper, egg, grey, hotmetal, jade, ...
       matcap, mateplastic, metal, metalclean, ...
       oldwax, prusa, planet, plastic, ...
       redshine, redwax, retro, rim, soft, ...
       skin, studioclay,bubble};
   
f = figure(101);
obj.bake().render();

for ii = 1:numel(mat)
   obj.Texture = mat{ii};
   obj.update;
   view(10,20); drawnow;
   pause(0.2);
end

function updateSphere(src,~,mat)
    class = whoClasses('Gmodel');
    for i = 1:length(class)
        class{i}.Texture = mat{round(src.Value)};
        r = update(class{i},'tex');
    end
    if src.Value >= 10
        title(['Diffuse value = ', num2str(round(src.Value),2)]);
    else
        title(['Diffuse value = ', num2str(round(src.Value),1)]);
    end
end