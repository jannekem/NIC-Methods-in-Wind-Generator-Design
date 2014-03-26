function [Q,Nph, zQ_final, A_final, a_final, k_wnu,nu,W] = winding_design(p,m,c_ws,tau_p,A_a,D,Is,l_ws,q)
% Winding design for PMSM, simple method
% vakoluku
% q = Q/(2*p*m);
% c_ws = 1;
W = c_ws*tau_p;

% number of coil turns in a phase
% Nph = ceil(A_a*pi*D/(2*Is*m));
% zQ = Nph*a/(p*q);
% l_ws
% q
%A_app = 10000:50:500000;

% Nph = ceil(A_app*pi*D/(2*Is*m));
% zQ_mod = mod(Nph*a,p*q);
% 
% Nph_pos = Nph(zQ_mod==0);
% zQ_pos = Nph_pos*a./(p*q);
% A_pos = 2*Is*Nph_pos*m/(pi*D);
% 

% 
% A_real = A_pos(abs(A_pos - A_a) == min(abs(A_pos - A_a)));
% Nph_real = Nph_pos(abs(A_pos - A_a) == min(abs(A_pos - A_a)));
% zQ_real = zQ_pos(abs(A_pos - A_a) == min(abs(A_pos - A_a)));
% 
% A_final = A_real(1); 
% Nph_final = Nph_real(1);
% zQ_final = zQ_real(1);
% Q = 2*p*m*q;

Q = 2*m*p*q;
% k = 1:10;
% a = 1:(p*q);
% Nph = 
% zQ_mod = mod(Nph*a,p*q);
% a_pos = a(zQ_mod==0);
% a_final = a_pos(abs(a_init-a_pos) == min(abs(a_init-a_pos)));
% zQ_final = Nph*a_final/(p*q);
% A_final = 2*m*Is*Nph/(pi*D);
a_app = 1:ceil(p/2);
zQ_app = ceil(A_a*pi*D*a_app/(Q*Is));
nph_mod = mod(zQ_app*p*q,a_app);
a_pos = a_app(nph_mod == 0);
zQ_pos = zQ_app(nph_mod == 0);
A_pos = Q*Is*zQ_pos./(pi*D*a_pos);
if isempty(a_pos)
    zQ_final = -1;
    Nph = -1;
    A_final = -1;
    k_wnu = -1;
    nu = -1;
    Q = -1;
    a_final = -1;
    W = -1;
    return
end



A_final = A_pos(abs(A_pos-A_a) == min(abs(A_pos-A_a)));
zQ_final = zQ_pos(abs(A_pos-A_a) == min(abs(A_pos-A_a)));
a_final = a_pos(abs(A_pos-A_a) == min(abs(A_pos-A_a)));
A_final = A_final(1);
zQ_final = zQ_final(1);
a_final = a_final(1);
Nph = zQ_final*p*q/a_final;


% winding factor
[ k_wnu, k_wfund, nu, nu_fund ] = ...
    integral_winding_factors( tau_p, c_ws, l_ws, 50, m, q, Q, p );


