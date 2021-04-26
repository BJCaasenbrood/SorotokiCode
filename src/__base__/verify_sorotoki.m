function out = verify_sorotoki
%VERIFY_SOROTOKI Summary of this function goes here
%   Detailed explanation goes here
clr;
fprintf('* Starting verification check - SOROTOKI \n'); pause(0.2)
fprintf('* The verification check might take a couple of minutes. \n');
fprintf('  The analysis will perform a full check shortly... \n');
pause(2);

List = cell(0);

% check plotting within SOROTOKI
List = verifySorotoki(List,@(x) checkPlot,'Plotting tools');
pause(1); close all;

List = verifySorotoki(List,@(x) check2,'Check 2');

out = List;
end

function checkPlot
figure(101);
subplot(2,2,1);
colorshow();

subplot(2,2,2);
[X,Y,Z] = peaks(40);
surf(X,Y,Z,'linestyle','none')
colormap(turbo); drawnow; pause(0.1);
colormap(barney); drawnow; pause(0.1);
colormap(blackwhite); drawnow; pause(0.1);
colormap(bluesea); drawnow; pause(0.1);
colormap(bounce); drawnow; pause(0.1);
colormap(inferno); drawnow; pause(0.1);
colormap(metro); drawnow; pause(0.1);
colormap(noir); drawnow; pause(0.1);
colormap(turbo); drawnow; pause(0.1);
colormap(viridis); drawnow; pause(0.1); 

subplot(2,2,3);
showColormap(turbo); drawnow; pause(0.1);
showColormap(turbo(-1)); drawnow; pause(0.1);
showColormap(turbo(0)); drawnow; pause(0.1);
showColormap(turbo(-100)); drawnow; pause(0.1);
showColormap(turbo(100)); drawnow; pause(0.1);

subplot(2,2,4);
imshow('FEMbot.png');
end

function list = verifySorotoki(list,exec,label)
list = verifyFunction(list,label,0);
try
    exec(0);
    list = verifyFunction(list,label,1);
catch e
    list = verifyFunction(list,e,-1);
end

end

function List = verifyFunction(List,input,index)
if index == 0
    cout('text','* ');
    cout('text','Assesing: ');
    str = char(input);
    cout('key',str);
    cout('key','... ');
    List{end+1,1} = str;
elseif index == 1
    List{end,2} = empty;
    cout(' ');
    cout('green','completed!\n');
else
    List{end,2} = input;
    cout(' ');
    cout('red','unsuccesfull!\n');
end

end