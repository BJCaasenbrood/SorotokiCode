clr; 
X = 0;
%% preview
obj = Gmodel(@(x) dSphere(x,0,0,0,1),[-1,1,-1,1,-1,1]);
obj.Texture = diffuse(X);

f = figure(101);
obj.bake().render();

%% setup control menu
b = uicontrol('Parent',f,'Style','slider','Position',[81,24,419,23],...
              'value',X, 'min',0, 'max',1);
          
bl1 = uicontrol('Parent',f,'Style','text','Position',[50,24,23,23],...
                'String','0','BackgroundColor',gitpage);
            
bl2 = uicontrol('Parent',f,'Style','text','Position',[500,24,23,23],...
                'String','1','BackgroundColor',gitpage);
          
b.Callback = @updateSphere;

title(['Diffuse value = ', num2str(X,3)]);

function updateSphere(src,evnt)
    class = whoClasses('Gmodel');
    for i = 1:length(class)
        class{i}.Texture = diffuse(src.Value);
        r = update(class{i},'tex');
    end
    title(['Diffuse value = ', num2str(src.Value,3)]);
end