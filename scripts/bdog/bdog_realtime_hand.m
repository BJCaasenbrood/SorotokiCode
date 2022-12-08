clr;
pause(3);
%% assign free DOF
usr = 'pi';
pwd = 'softroboticsSA';
ip  = '192.168.0.2';

brd = Bdog(usr,ip,pwd,'NVeab',3);

%% get data
brd = brd.set('Frequency',30);

%% execute control loop
phi = @(k) (k-1)*pi/6;

w = 4; a = -1;
while brd.loop(300)
    
    t = brd.t;
    
    Pd = zeros(1,6);
    %Pd(1) = -5;
    Pd(1) = -3;
    Pd(2) = -25;
    Pd(5) = 10;
    %Pd(1) = Pset*sign(sin(w*t - phi(1)));
    %Pd(2) = Pset*sign(sin(w*t - phi(2)));
    %Pd(3) = Pset*sign(sin(w*t - phi(3)));
    %Pd(4) = Pset*sign(sin(w*t - phi(4)));
    %Pd(5) = Pset*sign(sin(w*t - phi(5)));
    
    Pd = Pd*smoothstep(t);

    brd.setInput(Pd);
    fprintf(['t=',num2str(t,3),'\n']);
end

brd.resetInput();
brd.disconnect();




