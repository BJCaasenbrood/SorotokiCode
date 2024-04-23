clr;
%% assign free DOF
port = 8888;
usr = 'pi';
pwd = 'softroboticsSA';
ip  = '192.168.4.1';

brd = Control(usr,ip,pwd,'autoConnect',false,...
    'Port',port,'NVeab',2);

%% get data
brd = brd.set('Frequency',60);

Pset = 15;

%% execute control loop
while brd.loop(55)
    
    % read dat
    t = brd.t;
    Pd = zeros(1,4);   

    y = clamp(sign(sin(t+pi/2 +pi/2)),0,1);
    p1 = 70*(y*1.1 - 0.5);
    p2 = Pset*sigmoid(t)*sin(t - 0*pi + pi/2);
    p3 = Pset*sigmoid(t)*sin(t - pi + pi/2);

    Pd(1) = p2;
    Pd(2) = p3;
    Pd(4) = p1;

    brd.setInput(Pd);
end

brd.resetInput();
pause(0.1);

brd.disconnect();

%%
function [p1, p2, p3] = reference1(t)
    Pset = 15;

end