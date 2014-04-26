function kv = kinematic_viscosity_of_air(T)
% This function calculates the kinematic viscosity of air in temperature T.
% The used formula is valid for temperatures between 100 and 1600 K.
kv = -1.1555 .* 10^(-14) .* T.^3 + 9.5728 .* 10^(-11) .* T.^2 + 3.7604 .* 10^(-8) .* T - 3.4484 .* 10^(-6);

