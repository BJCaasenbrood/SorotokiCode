function [SHAPES, GRIPPER] = sinatra_gripper(varargin)

    clr;
    [Y,~] = chebyspace(80,3);

    currentDir = fileparts(mfilename('fullpath'));
    stlPath1 = fullfile(currentDir, '..', 'assets', 'stl', 'sinatra_gripper_base.stl');
    stlPath2 = fullfile(currentDir, '..', 'assets', 'stl', 'sinatra_gripper_holder.stl');    

    obj1 = Gmodel(stlPath1,'Shading','Face');
    obj1.Texture = matcap_grey * 1.25;
    obj1.render();
    
    obj2 = Gmodel(stlPath2,'Shading','Face');
    obj2.Texture = matcap_whitebase * 0.65;
    obj2.render();

    sdf = sSphere(45,[0,0,-65]);
    obj = Gmodel(sdf);
    obj.render();

    GRIPPER = {obj1,obj2,obj};
    clearvars obj1 obj2 obj

    Q = initialStates;
    G = initialBase;

    mat = NeoHookean(0.01,0.33);
    mat.contact.TangentFriction = 1.5;

    for ii = 1:6
        shp = Shapes(Y,[0,0,3,0,0,0],'Length',130,'ContactDistance',0.85,...
            'Texture',matcap_diffuse(.45)*1.15,'TimeStep',1/60);
        shp = shp.setMaterial(mat);
        shp = shp.setRadius([6,0.85,0,1]);
        shp = shp.setBase(G{ii});
        shp = shp.addDrag(100e-12);
        shp = shp.rebuild();

        shp = shp.addContact(sdf);
        shp.solver.sol.x = Q(ii,:).';
        shp.solver.TimeHorizon = Inf;
        shp.solver.TimeStep = 1/30;
        shp.solver.MaxIteration = 15;

        SHAPES{ii} = showRenderShapes(shp,Q(ii,:));
    end

    axis tight;
    view(190,-10);
end

function Q = initialStates

    Q = zeros(6,3);
    Q(:,1) = -1;
    Q(:,2) = -10;
    Q(1,1) = 1;
    Q(1,2) = -15;
    Q(2,1) = 14;
    Q(2,2) = 0;
    Q(3,1) = 2;
    Q(3,3) = -8;
    Q(4,1) = 2;
    Q(4,3) = -8;
    Q(5,1) = 8;
    Q(5,2) = -8;
    Q(6,1) = 20;
    Q(6,2) = 1;

    Q = -Q * 1e-3;
end

function G = initialBase
    G{1} = SE3(rotx(pi + 0.175),[18,20,20]) * SE3(roty(-pi/2));
    G{2} = SE3(rotx(pi + 0.175),[-18,20,20]) * SE3(roty(-pi/2));
    G{3} = SE3(rotz(-pi * 0.5),[-39,0,20]) * SE3(rotx(+0.0)) * SE3(roty(pi/2));
    G{4} = SE3(rotz(pi * 0.5),[+39,0,20]) * SE3(rotx(-0.0)) * SE3(roty(pi/2));
    G{5} = SE3(rotx(-pi - 0.175),[18,-20,20]) * SE3(roty(-pi/2)) *  SE3(rotx(+pi));
    G{6} = SE3(rotx(-pi - 0.175),[-18,-20,20]) * SE3(roty(-pi/2)) *  SE3(rotx(+pi));
end