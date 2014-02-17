function [ R_Fey ] = yoke_reluctance( w_y, mu_Fe, mu_0, L, h_y )
%This function evaluates the reluctance of the stator yoke
%   The flux in the yoke flows in the circumferential direction. Therefore,
%   for a rectangular shape, the reluctance is the width, w_y, of the yoke
%   divided by the product of its permeability, length and height, h_y.
R_Fey = w_y./(mu_Fe*mu_0*L*h_y);

end

