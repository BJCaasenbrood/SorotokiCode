clr;

% load variables
entrypoint_FKF;
NumberOfCalls = 1e3;

%% regular
time0 = cputime;clc

for ii = 1:NumberOfCalls
    
[g,J] = computeForwardKinematicsFast(x,... % states
    ds,...      % spatial steps
    p0,...      % position zero
    Phi0,...    % phi zeroclc
    xia0,...    % intrinsic strain vector
    Th,...      % evaluated Theta matrix
    Ba);

end

time0 = cputime - time0;
fprintf('Regular computation time = %2.3f\n',time0);

%% accelerated MEX
time1 = cputime;

for ii = 1:NumberOfCalls
    
[g,J] = computeForwardKinematicsFast_mex(x,... % states
    ds,...      % spatial steps
    p0,...      % position zero
    Phi0,...    % phi zeroclc
    xia0,...    % intrinsic strain vector
    Th,...      % evaluated Theta matrix
    Ba);
end

time1 = cputime - time1;
fprintf('Accelerated computation time = %2.3f\n',time1);