%clr; 
X = 1;
%% preview
obj = Gmodel(@(x) dSphere(x,0,0,0,1),[-1,1,-1,1,-1,1]);

mat = {aniso, bluered, bump, chroma, chromium, ...
       clean, copper, egg, grey, hotmetal, jade, ...
       matcap, mateplastic, metal, metalclean, ...
       oldwax, prusa, planet, plastic, ...
       redshine, redwax, retro, rim, soft, ...
       skin, studioclay};
   
f = figure(103);
obj.set('TextureStretch',0.65);
obj.bake().render();


b = uicontrol('Parent',f,'Style','slider','Position',[81,24,419,23],...
              'value',X,'min',1, 'max',length(mat));
          
bl1 = uicontrol('Parent',f,'Style','text','Position',[50,24,23,23],...
                'String','0','BackgroundColor',gitpage);
            
bl2 = uicontrol('Parent',f,'Style','text','Position',[500,24,23,23],...
                'String',num2str(length(mat)),'BackgroundColor',gitpage);
          
b.Callback = @(s,e) updateSphere(s,e,mat);

title(['Material capture number: ', num2str(X,3)]);

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