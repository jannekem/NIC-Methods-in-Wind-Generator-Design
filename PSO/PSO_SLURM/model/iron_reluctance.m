function [ R_Fe ] = iron_reluctance( R_Fet, R_Fey )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
R_Fe = R_Fet + R_Fey.*R_Fey./(R_Fey + R_Fey);

end

