% clr
% cl = Console();
% 
% cl.println('Hello world');
% 
% headers = {'\pi', 'Values', 'Residuals'};
% data = [(1:20).', rand(20,1), rand(20,1)*1e-3];
% 
% cl = cl.set('progressSteps',200,'progressName','Rendering');
% 
% for ii = 1:200
%     cl = cl.progress(ii);
%     pause(1e-3);
% end
% 
% cl.pause(3);
% cl.table(headers,data,5);
% 
% A = 'SOROTOKI - a toolkit for design and control of soft robots';
% C1 = ['Current version: ', soropatch(1), ' (',soropatch(2),')'];
% C2 = 'Written by: b.j.caasenbrood - TU/e';
% C3 = ['Modified on: ', soropatch(3)];
% 
% cl.table({A},{C1;C2;C3})
clr;

M = 1;
N = 150;

Y = pccspace(N,M);

shp = Shapes(Y,[M,M,M,0,0,0],'Length',60);
shp = shp.setBase(roty(-pi/2));

fig(101);
view(60,30);
% 
% pointA = [0,-20,0];
% pointB = [0,20,0];
% pointC = [0,20,50];
% pointD = [0,-20,50];
% points=[pointA' pointB' pointC' pointD']; % using the data given in the question
% fill3(points(1,:),points(2,:),points(3,:),'b','LineStyle','none');
% grid on
% alpha(0.3)

for ii = 1:5

    shp = shp.render([25 * (ii-1),0,-25]*1e-3);
    shp.Gmodel = shp.Gmodel.set('Shading','Face');
    shp.Gmodel.ground([-20,20,-20,20,0,1]);
    shp.Gmodel.update();
    axis([-30,30,-30,30,0,50]);
    h = plotSE3(shp.Log.FK.g(:,:,end),'arrowsize',10,'linewidth',3);
    %pause();
    export_fig(['img',num2str(ii),'.png'],'-q100','-a4','-r300');
    delete(h{1});
    delete(h{2});
    delete(h{3});
end
