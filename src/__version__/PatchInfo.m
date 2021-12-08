%% ------------------------------------------------------------------------
%% SOROTOKI - Early Alpha - 1.12.06 - Dec 6 - 2021
% - Fixed some numerical issues with Fem.simulate() function. Stability is 
%   now further improved.
%
% - Fixed some numerical issues with Fem.Contact. A wider range of
%   TimeSteps result now in stable solutions. Still the timestep size
%   should be taken with care. Also fixed the issues of induced
%   (unstable) oscillations due to fast impact. Contact is now based on the
%   intial Modulus of the hyper-elastic material,i.e., Material.Emod().
%
% + Added time-base function_handle for the external pressures in dynamic
%   finite-element simulations. They can added using:
%       3) Fem.AddConstraint('Pressure',id, @(t) sin(t));
%
% + Dynamic simulations now record the potential and kinetic energies.
%
% ? Future implementation will have Load and Tendon-based dynamic forces
%
% ? The potential energy of the extenal load is still missing...
%% ------------------------------------------------------------------------
%% SOROTOKI - Early Alpha - 1.12.02 - Dec 2 - 2021
% + Added Fem.simulate() function to the Nonlinear Finite Element. The
%   function is a standard Newmark-Beta dynamic solver routine for NLFEM.
%
%   Alternatively, fem.simulate() can be used as an alternative for
%   fem.solve() by setting the Fem.Material.Zeta large (overdamped).
%
% - Fixed the stability issue of Fem.Contact. Contact can now be added with
%   Fem.AddConstraint('Contact',@(x) sdf(x), [x,y]) (x and y move the SDF). 
%   for best stability, use Fem.simulate() with a TimeStep < 1/75.
%
% ? Open issue: Fem.simulate does not use dynamic external forces, rather
%   the forces are scaled with a sigmoid fucntion f_ext = f*sigmoid(t);
%   This will be removed in the future with fext = @(x) f* .... - an
%   explicit function of time that can be evaluated at time t.
%       Possible implementations: Displace, Load, Pressure, Gravity.

function x = PatchInfo(Request)
switch(Request)
    case(1); x = '1.12.02';
    case(2); x = '2 December 2020';
    case(3); x = '30 August 2018';
end
end

