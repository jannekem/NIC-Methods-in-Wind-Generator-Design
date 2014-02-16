function Lm = magnetizingInductance(m, alpha_i, p,tau_p, delta_eff,lef,k_w1,Nph)
% in case of a salient pole machine, this function can be used to 
% calculate the magnetizing inductance for both d- and q-axis.
mu0 = 4*pi*10^-7;
Lm = (m/2)*(4/pi)*alpha_i*mu0*(1/(2*p))*(tau_p/delta_eff)*lef*(k_w1*Nph)^2;