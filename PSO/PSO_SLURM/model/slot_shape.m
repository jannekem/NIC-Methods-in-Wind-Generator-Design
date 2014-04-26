function [h1,h2,h3,h_tooth, h_slot, b1,b2,b3,bs1,bs2] = slot_shape(Q, tau_u, bs1_rel, bs2_rel, h1_rel, h2_rel, D, A_slot)
% slot shape construction

bs1 = tau_u*bs1_rel;
b1 = tau_u*(1-bs1_rel);
bs2 = bs1*bs2_rel;

h3 = A_slot/bs2;
h1 = h1_rel*h3;
h2 = h2_rel*h3;
h_tooth = h1+h2+h3;
h_slot = h_tooth;

tau_u2 = (D + 2*(h1+h2))*pi/Q;
b2 = tau_u2 - bs2;

tau_u_3 = (D + 2*h_tooth)*pi/Q;
b3 = tau_u_3 - bs2;


% aa = pi/Q + pi^3/(2*Q^2);
% bb = bs2+pi^2/(2*Q)*bs2;
% cc = pi/8*bs2^2 - A_slot;
% 
% x3_1 = (-bb + sqrt(bb^2 - 4*aa*cc))/(2*aa);
% x3_2 = (-bb - sqrt(bb^2 - 4*aa*cc))/(2*aa);
% 
% h3 = max(x3_1,x3_2);
% bs3 = 2*h3*pi/Q + bs2;
% b2 = (D + 2*(h1_rel + h2_rel)*h3)*pi/Q - bs2;
% h4 = bs3/2;
% h_tooth = h1+h2+h3+h4;
% h_slot = h_tooth;

