clr;
%% construct dynamic
mdl = Model([0,1,1,1,0,0],'NModal',1);
mdl = mdl.generate();

%% set inputs
mdl.q0 = [0,2,-0.1];
mdl.Pressure = [2;0;0];

mdl = mdl.csolve(); 

mdl.showModel();
