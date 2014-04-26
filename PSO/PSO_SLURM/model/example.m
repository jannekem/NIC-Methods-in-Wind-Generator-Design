% An example of how to run design_PMSM_generator -function.
% example_params.mat contains an example set of input parameters for the
% design algorithm. 

% The function design_PMSM_generator gives efficiency and torque density as
% outputs. If a constraint is violated inside the design_PMSM_generator
% funciton, it returns -1 for both outputs.

load example_params.mat
[objects, constraints] = design_PMSM_generator(params);


