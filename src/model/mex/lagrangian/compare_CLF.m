clr;

% load variables
entrypoint_CLF;
NumberOfCalls = 1e3;

%% regular
time0 = cputime;

for ii = 1:NumberOfCalls
    
computeLagrangianFast(x,dx,... % states
    ds,...      % spatial steps
    p0,...      % position zero
    Phi0,...    % phi zero
    xia0,...    % intrinsic strain vector
    Th,...      % evaluated Theta matrix
    Ba,...      % state to strain matrix
    Ktt,...     % geometric stiffness
    Mtt,...     % geometric inertia
    Zeta,...
    [0;0;0]);

end

time0 = cputime - time0;
fprintf('Regular computation time = %2.3f\n',time0);

%% accelerated MEX
time1 = cputime;

for ii = 1:NumberOfCalls
    
computeLagrangianFast_mex(x,dx,... % states
    ds,...      % spatial steps
    p0,...      % position zero
    Phi0,...    % phi zero
    xia0,...    % intrinsic strain vector
    Th,...      % evaluated Theta matrix
    Ba,...      % state to strain matrix
    Ktt,...     % geometric stiffness
    Mtt,...     % geometric inertia
    Zeta,...
    [0;0;0]);
end

time1 = cputime - time1;
fprintf('Accelerated computation time = %2.3f\n',time1);