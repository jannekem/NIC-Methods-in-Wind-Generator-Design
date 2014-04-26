function [ rho_T_op ] = resistivity_at_T_op( T_op, T_ref, rho_ref, alpha_0 )
% calculates the resistivity of a material in a certain temperature
rho_T_op = rho_ref*(1 + alpha_0*(T_op-T_ref));
