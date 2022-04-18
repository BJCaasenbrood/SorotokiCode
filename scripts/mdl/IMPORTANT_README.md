To efficiently run the Model.m class, make sure you have compiled the mex files.

1) Type 'coder' in the command window, check if it is installed!
2) Nagivate to SorotokiCode/src/mdl/mex
3) You find two folders: forwardkin and lagrangian

A: forwardkin
A1) open forwardkin folder, make sure the directory is:
    XX/SorotokiCode/src/model/mex/forwardkin
A2) type 'coder' in command window, a new window opens
A3) type 'ComputeForwardKinematicsFast' in Generate code, hit enter
A4) Click 'reopen' in the yellow box
A5) A new window opens, hit enter again
A6) It should say 'entrypoint_FKF' under the Define Input Types menu, hit next
A7) Matlab will now check the code, this can take a while, then hit next
A8) Last step, hit generate; either C or C++, it doesnt matter significantly.
A9) Matlab will build mex, and it should say Build succeeded
A10) Hit next, and close window; that it done!

B: lagrangian
B1) open lagrangian folder, make sure the directory is:
    XX/SorotokiCode/src/model/mex/lagrangian
A2) type 'coder' in command window, a new window opens
A3) type 'computeLagrangianFast' in Generate code, hit enter
A4) Click 'reopen' in the yellow box
A5) A new window opens, hit enter again
A6) It should say 'entrypoint_CLF' under the Define Input Types menu, hit next
A7) Matlab will now check the code, this can take a while, then hit next
A8) Last step, hit generate; either C or C++, it doesnt matter significantly.
A9) Matlab will build mex, and it should say Build succeeded
A10) Hit next, and close window; that it done!

Last step:
4) Open Model.M under XX/SorotokiCode/src/model/
5) set 'MexSolver = true' if it was false. Model.m will now use compiled code
