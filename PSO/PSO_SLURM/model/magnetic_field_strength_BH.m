function [ H_v ] = magnetic_field_strength_BH( B_data, H_data, B_v )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
H_v = spline(B_data,H_data,B_v);

end

