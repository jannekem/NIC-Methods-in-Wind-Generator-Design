% function [ k_w, k_wfund, k_w2, nu, nu_fund ] = ...
function [ k_w, k_wfund, nu, nu_fund ] = ...
    integral_winding_factors( tau_p, c_w, l_w, c, m, q, Q, p )
%This function evaluates the winding factors for an integral slot winding
%   Detailed explanation goes here
if ( l_w == 1 )
    %disp('Single layer winding')
    % Short pitching of single layer winding not possible. Different
    % overhangs, i.e. the connections of coil sides, however, can be
    % employed such as diamond, concentric...
    W = 1.*tau_p;
else
    %disp('Double layer winding')
    W = c_w.*tau_p;
end
g = -c:1:c;
nu = 2.*m.*g + 1;
nu_fund = 2.*m.*q.*g + 1;
% k_w2 = sin(nu.*W./tau_p.*pi()./2).*sin(nu.*pi()./(2.*m))...
%     ./(q.*sin(nu.*pi()./(2.*m.*q)));
k_wfund = 2.*sin(nu_fund.*pi().*W./(2.*tau_p)).*sin(nu_fund.*pi()./(m.*2))...
    ./(Q./(m.*p).*sin(nu_fund.*pi().*p./Q));
k_w = 2.*sin(nu.*pi().*W./(2.*tau_p)).*sin(nu.*pi()./(m.*2))...
    ./(Q./(m.*p).*sin(nu.*pi().*p./Q));

end