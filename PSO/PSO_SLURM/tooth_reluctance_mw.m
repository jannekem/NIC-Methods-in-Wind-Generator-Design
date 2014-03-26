function [ R_Fet ] = tooth_reluctance_mw( h_s, mu_Fe, mu_0, L, zmin,zmax)
%This function evaluates the reluctance of a tooth
%   The flux in the teeth flows in the radial direction. Therefore, for a
%   rectangular shape, the reluctance is the height, h_s, of the teeth
%   divided by the product of their permeability, length and width.
%R_Fet = h_s./(mu_Fe*mu_0*L*z);

R_Fet = 1/(mu_Fe*mu_0*L)*(h_s./(zmax-zmin)).*log(zmax./zmin);

end

