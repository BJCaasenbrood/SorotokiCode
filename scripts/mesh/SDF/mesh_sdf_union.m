clr;
%% generate SDF class
sdf1 = sCircle(0, 0, 0.5);
sdf2 = sCircle(0.5, 0, 0.5);

%% math operations
S = sdf1+sdf2;
S.show()
