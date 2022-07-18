%% ------------------------------------------------------------------------
% % SOROTOKI - Alpha - 2.05.19 - May 19 - 2022
% **Fem.m**
%
% - @muyi-max: Fixed issue with Fem.Log. If a nonlinear residual increment takes less
%   than 2 iterations, data was not properly stored. To solve this, any
%   nonlinear time increment must take atleast two steps, as to ensure a
%   new displacement field is introduced to the stress-strain calculations.
%
%% ------------------------------------------------------------------------
% % SOROTOKI - Alpha - 2.05.17 - May 17 - 2022
% + Added MEX compiler for Model, Gmodel, and Fem classes. This signicantly
%   improves the computation times of the Sorotoki toolkit. When calling
%   sorotoki.m, the installer askes to compile the mex files. *Note* we do
%   not support the orginal Matlab functions anymore, simple because of
%   computational speed. Hence, the Matlab Coder Toolkit is a prerequisite 
%   for Sorotoki.
% 
% + Updated verifySorotoki to include also Model.
% 
% **Shapes.m**
% + Updated the Fem -> Shapes converter that extracts the geometric modes
%   from dynamic or quasi-static FEM data. This can be done using shp.Fem =
%   Fem followed by shp.reconstruct, shp.rebuild. Then a dynamic model can
%   be generated using mdl = Model(shp). We call this new strain basis the
%   **Geometry-Informed Variable Strain** basis as the geoemtry of the soft
%   robot is preserved into the functional basis.
% 
% + Added Shapes.Material. Now hyper-elastic materials can be loaded into
%   the stiffness and interia tensor construction. The work is experimental
%   since more research is required to properly include hyp-elast materials
%   into the Cosserat beam models.
% 
% **Model.m**
% + Added an auxiliary flow function to the Model class: dx = f(x,t). This
%   can be used for integrator actions, like PI controllers. or adaptive
%   controllers where unknown parameters evolutions can be included.
% 
% % ------------------------------------------------------------------------
% % SOROTOKI - Alpha - 2.01.28 - Jan 28 - 2022
% + Improved stability of Newmark solver.
% 
%% ------------------------------------------------------------------------
% % SOROTOKI - Alpha - 2.01.24 - Jan 24 - 2022
% - Fixed broken installer. `vernum.m` file was missing on Repo. It has been
%   replaced with `soropatch.m` which also includes the patch notes.
%
% + **Sdf.m**:
%   - Added a new function `Sdf.showcontour()`. It will show the contour of 2D
%   signed distance functions. Currently not implemented for 3D Sdfs.
% 
% + **Open issue (2.24.01)**:    
%   - @martijnschouten Missing DOI for citation, and long-term support/access.
%% ------------------------------------------------------------------------
% % SOROTOKI - Alpha - 2.01.13 - Jan 13 - 2022
% + Moved SOROTOKI from early alpha to alpha (prepping for official release).
% 
% + Signifcant update to the class Shapes. Shapes now requests a Fem class
%   , a Matlab function_handle of functionals (e.g., ), or a evaluated matrix of
%   shape functions (); to construct the class. Then, it serves as an input
%   for the class Model (e.g., dynamic modeling of soft robots), or as an
%   Inverse Kinematic solver for beam models. 
% 
% + Added a function Shapes.jointEstimate(Fem). This will produce an
%   optimal set of modal coefficient to reconstruct the beam model from
%   FEM-driven data (acccessed by Fem.Log.(...) ). 
% 
% - Current Shapes.jointEstimate() exploits the fact that the shapes are
%   orthonormal w.r.t. their integral over the spatial domain [0,L]. If  
%   shapes do not satisfy this conditions, inaccurate estimates can occur. 
% 
% + Added Y = gsog_poly(X) which ensure the columns of Y are mutually
%   orthonormal derived from the matrix X. X must be full-rank. This
%   function is simply the gramm-smith method for the innerproduct space
%   int_C y_i y_j ds on the domain C := [0,L].
%% ------------------------------------------------------------------------
% % SOROTOKI - Early Alpha - 1.12.06 - Dec 6 - 2021
% - Fixed some numerical issues with Fem.simulate() function. Stability is 
%   now further improved for larger timesteps.
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
% % SOROTOKI - Early Alpha - 1.12.02 - Dec 2 - 2021
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

function x = soropatch(varargin)
if isempty(varargin)
   varargin{1} = 1; 
end

switch(varargin{1})
    case(1); x = '2.06.02';
    case(2); x = 'May 19, 2022';
    case(3); x = 'August 30, 2018';
end
end

