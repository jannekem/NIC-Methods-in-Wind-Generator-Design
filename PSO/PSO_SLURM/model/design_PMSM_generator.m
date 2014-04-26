function [objectives, constraints]  = design_PMSM_generator(optim_params)
P = 3e6; % output power W
U_ph = 690/sqrt(3); % phase voltage V
n = 16.98/60; % rotation speed (1/s)
%% Design parameters ( these can be used for optimisation )
% pole pair number
p = optim_params(1); 
% Magnetic flux density in the air-gap (peak) T
A_a = optim_params(2);
%Magnetic flux densities in the stator teeth and yoke 
B_zinit = 1.5;%optim_params(3);
B_yinit = 1.2;%optim_params(4);
% Current density in the stator winding A/m^2
J = optim_params(3); 
% Air gap length (m)
delta = optim_params(4);
% Machine length/air-gap diameter
% Dr/lstk
khii = optim_params(5);
% relative magnet width (= magnet width / pole pitch)
lp_rel = optim_params(6);
% Machine constant (mehcanical) 
% tangential stress
sigma_Ftan = optim_params(7);
% magnetic flux density in the rotor yoke
B_yrmax = optim_params(8);%optim_params(1.3 - 1.6 T); 
% number of parallel paths
%a_init = optim_params(9);
% number of slots per pole per phase
q = optim_params(9);
% number of conductors in a slot
Dse_rel = 1/optim_params(10);
%
bs1_rel = optim_params(11);
bs2_rel = 1/optim_params(12);
h1_rel = optim_params(13);
h2_rel = optim_params(14);
%sigma_Ftan = optim_params(16);
c_ws = 1;
l_ws = 1;
%% BH-curve for the stator and rotor iron
load BH_SB_M600_50A.mat

% h2 = 3e-3;%% Constants and fixed initial values

% for the efficiency and power factor (these are calculated later)
eta = 0.94; % efficiency (desired)
cosphi = 0.95; % power factor
Pel = P; % electrical power
Ps = P/eta;
% number of phases
m = 3;
%iron space factor 
fr = 0.95;
%copper space factor
kCu = 0.6; 
% ventilation ducts 
nv = 0; % number of ventilation ducts
bv = 6e-3; % width of ventilation ducts
c_EMF = 1.05;
% % slot dimensions
% h1 = 3.5e-3;

% h3 = 0;
% b1rel = 0.75; % relative width of the slot opening

%% material properties
mu0 = 4*pi*10^-7; %permeability of vacuum
% magnet properties
Hc = 800000; %coersive field strength 
muPM = 1.05; % relative permeability of the magnet material
rhoPM = 7400; % density of the permanent magnet material
% density of iron
rho_Fe = 7870;
% density of copper
rho_Cu = 8960; 
% endwinding leakage inductance coefficients
%c_w = 1; tämä on jänteistys 
lambda_W = 0.166;
lambda_lew = 0.371;

%coefficients for iron losses
k_Fey = 1.5;
k_Fez = 2;
p10 = 1.75;
p15 = 6.6;
k_rho = 10; % windage and ventilation loss factor
mu = 0.003; %friction coefficient
%% Material prices
Cu_price = 12; % Copper price (€/kg)
Fe_price = 4; % Electrical steel price (€/kg)
Magnet_price = 60; % Price of NdFeB magnets (€/kg)

kLoss = 2; % Loss price (€/W)

% FREQUENCY AND ROTATION SPEED
f = p*n; % supply frequency
omega = 2*pi*f;
Omega = 2*pi*n;
T_est = Ps/Omega;
Vr_est = T_est/(2*sigma_Ftan);
l = nthroot(4*Vr_est/(khii^2*pi),3);
% STATOR LENGTH & INNER DIAMETER    
% stator inner diameter = air-gap diameter
kappa = (bv/delta)/(5+bv/delta); % coefficient for calculating the effective 
% ventilation duct width bve
bve = kappa * bv; 
% % real machine length
% l = lef + nv*bve - 2*delta; 
lef = l - nv*bve + 2*delta;
Dr_out = l*khii;
D = Dr_out + 2*delta;
%% winding design
% WINDING
% pole pitch
tau_p = pi*D/(2*p);
% estimated stator current
Is = Pel/(m*U_ph*cosphi);
% conductor area
A_cond = Is/J;
%A_a = sigma_Ftan*sqrt(2)/(B*cosphi);
[Q,Nph, zQ, A,a, k_wnu,nu,W] = winding_design(p,m,c_ws,tau_p,A_a,D,Is,l_ws,q);
% if Nph == -1
%     %cost = -1;
%     efficiency = -1;
%     torque_density = -1;
%     return
% end
k_w1 = k_wnu(nu == 1);
tau_u = pi*D/Q;
% if tau_u < 7e-3
%     %cost = -1;
%     efficiency = -1;
%     torque_density = -1;
%     return
% end
%constr1 = tau_u; % slot pitch has to be at least 7 mm



% SLOT AND TOOTH DIMENSIONS
% slot area
A_slot = A_cond*zQ/kCu;
[h1,h2,h3,h_tooth, h_slot, b1,b2,b3,bs1,bs2] = slot_shape(Q, tau_u, bs1_rel, bs2_rel, h1_rel, h2_rel, D, A_slot);
% stator yoke inner diameter
Dsy = D + 2*h_tooth;
% stator yoke
Dse = Dse_rel*Dsy;
% stator outer diameter
hy = (Dse - Dsy)/2;
by = m*q*pi*(D+2*h_tooth+hy)/Q; %yoke width 

% MAGNET
E_PM_req = c_EMF*U_ph;
B_est = sqrt(2)*E_PM_req/(2*pi*f*k_w1*Nph*lef*tau_p);% ilmavälin B
% h_mag = delta*B/(mu0*(Hc - B/lp_rel/(muPM*mu0*lp_rel)));
% B = 0.8;
%B_est = sigma_Ftan*sqrt(2)/(A*cosphi);
%lm = delta*B_est/(mu0*Hc - B_est/(muPM*lp_rel));
lm = delta*lp_rel*B_est/(mu0*(Hc - B_est/(mu0*muPM)));
%constr2 = Hc - B_est/(mu0*muPM); % > 0
delta_lm = lm;
if Hc < B_est/(mu0*muPM) % magnet material not strong enough to produce required air-gap flux density (possibly not enough coil turns)
    %cost = -1;
%     efficiency = -1;
%     torque_density = -1;
    lm = 0.5;
    delta_lm = 1e-3;

end

% h_mag = 8e-3;
%D_r = D-2*delta-2*lm; % rotor yoke outer diameter
lp = tau_p*lp_rel;
%% air-gap flux and flux density, iteration


%lm = 0.5e-3;
E_PM = 0;
lm = lm - delta_lm;
% magnetic field strengths (maximum values)
H_zmax = magnetic_field_strength_BH(B_data,H_data,B_zinit);
H_ymax =  magnetic_field_strength_BH(B_data,H_data,B_yinit );
% relative permeability of iron, teeth and yoke
muFez = (B_zinit/H_zmax)/mu0;
muFey = (B_yinit/H_ymax)/mu0;
while abs(E_PM - c_EMF*U_ph) > 0.05*c_EMF*U_ph
    lm = lm + delta_lm;    
    ThetaPM = Hc*lm;
    %calculation of reluctances
    R_delta = delta/(mu0*tau_p*lef);
    R_delta2 = lm/(mu0*(tau_p-lp)*lef);
    R_PM = lm/(mu0*muPM*lp*l);
    %Phi_deltaPM = Bp_delta(:,k)*alpha_PM*tau_p*lef;
    B_delta_prev = -1;
    Bp_delta = 2;
    ii_loop_stop = 0;
    while abs(Bp_delta - B_delta_prev)/Bp_delta > 0.01
        ii_loop_stop = ii_loop_stop + 1;
        B_delta_prev = Bp_delta;
        %disp('loop')
        % new iron reluctance
        R_Fey = yoke_reluctance(by, muFey, mu0, l, hy );
        R_Fez = tooth_reluctance_mw(h_tooth, muFez, mu0, l, b2,b3);
        R_Fe = iron_reluctance(R_Fez, R_Fey/2);
        % new flux
        Phi_deltaPM = airGapFlux(ThetaPM, R_PM, R_Fe, R_delta, R_delta2);
        %Phi_deltaPM = (4/pi)*Phi_deltaPM*sin(lp*pi/(2*tau_p)); %lisätty huippuarvon estimointi
        Bp_delta = Phi_deltaPM/(lp_rel*tau_p*lef);
        B_z = Phi_deltaPM/(3*q*(b2+b3)/2*l);
        B_y = Phi_deltaPM/(2*l*by);
        H_z = magnetic_field_strength_BH(B_data,H_data,B_z);
        H_y = magnetic_field_strength_BH(B_data,H_data,B_y);
%         H_z = HB_fitted_curve(B_z);
%         H_y = HB_fitted_curve(B_y);
        muFez = B_z/H_z/mu0;
        muFey = B_y/H_y/mu0;
        if ii_loop_stop > 100
            break
        end
    end
    E_PM_prev = E_PM;
    E_PM = 2*pi*f*k_w1*Nph*Phi_deltaPM/sqrt(2);
    if isnan(E_PM)
        E_PM = 0;
    end
    if E_PM_prev < 0.95*c_EMF*U_ph && E_PM > 1.05*c_EMF*U_ph
       lm = lm - 2*delta_lm;
       delta_lm = 0.5*delta_lm;
       if delta_lm < 1e-10
           break;
       end
    end
    if lm > 0.5
        break;
    end
end

    
%     if isnan(E_PM)
%         cost = -1;
%         temperatures = -1;
%         return
%     end
%     ind = find(abs(E_PM - c_EMF*U_ph) == min(abs(E_PM - c_EMF*U_ph)));
%     E_PM = E_PM(ind);
%     if abs(E_PM - c_EMF*U_ph) > 0.05*c_EMF*U_ph
%         %cost = -2;
%         efficiency = -1;
%         torque_density = -1;
%         return
%     end
%constr3 = abs(E_PM - c_EMF*U_ph)/(c_EMF*U_ph);

Dr = D - 2*delta - 2*lm;
hyr = (Phi_deltaPM/B_yrmax)/l;
Dri = Dr - 2*hyr;
if Dri > Dr/3
    Dri = Dr/3;
end
h_H = Dse/2 + 0.02; % frame radius
l_E = 0.2;


%% Inductances


kappa = (b1/delta)/(5+(b1/delta));
kC = tau_u/(tau_u-kappa*b1);
%kC = 1.07;
delta_eff = kC*delta + lm/muPM;
Lmd = magnetizingInductance(m, lp_rel, p,tau_p, delta_eff,lef,k_w1,Nph);
Lmq = magnetizingInductance(m, lp_rel, p,tau_p, delta_eff,lef,k_w1,Nph);

%Leakage inductances
% air-gap leakage inductance

% c = [-50:-1 1:50];
% nu = 1+2*c*m;
% k_wnu = zeros(1,length(nu));
% for j = 1:length(nu)
%     k_wnu(1,j) = windingFactor(nu(1,j),W,tau_p,p,Q,m);
% end
k_wnu_1 = k_wnu(nu~=1);
sigma_delta = sum((k_wnu_1./(nu(nu~=1).*k_w1)).^2);
L_delta = sigma_delta*Lmd;

% slot leakage inductance
lambda_u = h3/(3*bs2)+ 0/bs2 + h1/bs1 + h2/(bs2-bs1)*log(bs2/bs1);
Lu = (4*m/Q)*mu0*lef*Nph^2*lambda_u; 

% tooth tip leakage inductance
epsilon = 1 - (W/tau_p);
k2 = 1 - epsilon;
lambda_tt = k2*(5*(delta/b1)/(5+4*(delta/b1))); 
Ltt = (4*m/Q)*mu0*lef*lambda_tt*Nph^2;

% end winding leakage inductance

Wav = tau_p*c_ws; % average width of endwinding
lw = (2.8*Wav - 2*Wav + 0.4)/2; % average length of end winding
taup = pi*(D+h_tooth)/(2*p);
tauu = pi*(D+h_tooth)/Q;
WeW = taup - tauu;
lew = 0.5*(lw-WeW); % axial length of endwinding
lwlambdaw = 2*lew*lambda_lew + WeW*lambda_W;
Lw = (4*m/Q)*q*Nph^2*mu0*lwlambdaw;

Lstr = L_delta + Ltt + Lu + Lw;


Ld = Lmd+ Lstr;
Lq = Lmq+ Lstr;

% reactances
Xq = omega*Lq;
Xd = omega*Ld;

 %% Resistances

lav = 2*l + 2.8*Wav + 0.4; %average length of a coil turn (large machines with prefabircated winding)

lc = lav*Nph;
sigma_Cu = 57e6; %conductivity of commercial copper wire (S/m) (room temperature 20 degrees)
alpha_Cu = 3.81e-3; % temperature coefficient for resistivity of copper (1/K);
rho_Cu_100 = resistivity_at_T_op(100,20,1/sigma_Cu,alpha_Cu); 
sigma_c = 1/rho_Cu_100;
R_DC = lc/(sigma_c*a*A_cond); % DC-resistance

hc0 = A_cond/bs2;
bc0 = bs2;
za = 1;
zt = zQ/za;
xi = hc0*sqrt(0.5*omega*mu0*sigma_c*bc0/bs2);
kR = fi(xi)+ ((zt^2-1)/3)*psi1(xi);
%kR = 1 + ((zt^2-0.2)/9)*xi^4;
R_AC = kR*R_DC; % AC-resistance


%% Load angle
% if 3*U_ph*E_PM/Xd < Pel
%     %cost = -3;
%     efficiency = -1;
%     torque_density = -1;
%     return
% end
Pmax = 3*U_ph*E_PM/Xd;
initial_value = Pel - 3*(U_ph*E_PM/(omega*Ld)*sin(pi/2)+U_ph^2*(Ld-Lq)/(2*omega*Ld*Lq)*sin(2*pi/2));
if isnan(initial_value) || isinf(initial_value) || ~isreal(initial_value)
%     %cost = -3;
%     efficiency = -1;
%     torque_density = -1;
%     return     
    delta_a = pi/2;
elseif Pmax < Pel
    delta_a = pi/2;
else
    delta_a = fzero(@(delta_a) Pel - 3*(U_ph*E_PM/(omega*Ld)*sin(delta_a)+U_ph^2*(Ld-Lq)/(2*omega*Ld*Lq)...
    *sin(2*delta_a)),pi/2,optimset('Display','off'));
end
% 3*U_ph*E_PM/Xd
% delta_as = 0:0.001:2*pi;
% P_el_delta_as = 3*(U_ph*E_PM/(omega*Ld)*sin(delta_as)+U_ph^2*(Ld-Lq)/(2*omega*Ld*Lq)*sin(2*delta_as));
% % min(abs(Pel - P_el_delta_as))
% delta_a = delta_as(abs(Pel - P_el_delta_as) == min(abs(Pel - P_el_delta_as)));

%% Currents

Id = (U_ph*(omega*Lq*cos(delta_a)-R_AC*sin(delta_a)) - E_PM*omega*Lq)/(omega^2*Lq*Ld + R_AC^2);
Iq = (U_ph*(R_AC*cos(delta_a) + omega*Ld*sin(delta_a)) - E_PM*R_AC)/(omega^2*Lq*Ld+R_AC^2);

In = sqrt(Id^2+Iq^2);

Pout = m*U_ph*(Iq*cos(delta_a) - Id*sin(delta_a)); 
%constr5 = abs(Pout-Pel)/Pel;
% if abs(Pout-Pel) > Pel*0.05 || isnan(delta_a)
%     %cost = -3;
%     
%     efficiency = -1;
%     torque_density = -1;
%     return
% end
S_x = 3*U_ph*In;
cos_phi = Pout/S_x;


J = In/(A_cond*a);
A = m*Nph*In/(pi*D/2);


%% Losses and efficiency
% resistive losses
P_Cu = m*R_AC*In.^2;


% teeth iron losses
m_Fez = Q*l*b2*h_tooth*rho_Fe*fr; % mass of teeth
P_z = m_Fez*k_Fez*p15*(B_zinit/1.5)^2*(f/50)^(3/2);
% yoke iron losses
m_Fey = (pi/4)*(Dse^2-(D+2*h_tooth)^2)*l*rho_Fe*fr; %mass of yoke
P_y = m_Fey*k_Fey*p15*(B_yinit/1.5)^2*(f/50)^(3/2);
% total iron losses
P_Fe = P_y + P_z;

% mechanical losses
m_rotor = pi*(Dr/2)^2*l*rho_Fe;
F_bearing = m_rotor*9.81;
P_Br = 2*0.5*Omega*mu*F_bearing*Dri;

% windage and ventilation losses
P_rho = k_rho*Dr*(l + 0.6*tau_p)*(2*pi*n*Dr/2)^2;

% stray losses
P_Str = 0.0075*P;

% winding losses
P_w = P_Cu*2*l/lav;
% end winding losses
P_ew = P_Cu*2*lw/lav;
% eddy current losses
delta_PMec = delta + lm/(2*muPM);
omega_nu = 2*pi*n*Q;
sigma_PM = 670000; % conductivity of the permanent magnet material
alpha_PM = lp/tau_p; % relative permanent magnet width
k_nu = sqrt(omega_nu*muPM*mu0*sigma_PM/2);
beta_nu = omega_nu/(pi*D*n);
alpha_nu = (1 + 1i)*k_nu;
a_Rnu = (1/sqrt(2))*sqrt(sqrt(4 + (beta_nu/k_nu)^4) + (beta_nu/k_nu)^2);

kappa = (b1/delta_PMec)/(5 + b1/delta_PMec);
be = kappa * b1;
kC = tau_u/(tau_u - kappa*b1);
u = b1/(2*delta_PMec) + sqrt(1 + (b1/(2*delta_PMec))^2);
beta = (1 + u^2 - 2*u)/(2*(1+u^2));

B0 = beta*Bp_delta;

S_PM = pi*D*alpha_PM*l;
P_ec = (1/2)*a_Rnu*(1 + tau_u/(2*l))*abs(alpha_nu)^2/beta_nu^2*(B0/(mu0*muPM))^2*k_nu/sigma_PM*S_PM;


% total losses
P_loss = P_Cu + P_Fe + P_rho + P_Str + P_Br + P_ec;

Ps = Pout + P_loss; %shaft power
% efficiency 
eta2 = Pout/Ps;

% torque
T = Ps/Omega;

thermal_analysis_2;
temperatures = Theta - 273.15; % temperature in Celcius
% if temperatures(6) > 100
%     %cost = -4;
%     
%     efficiency = -1;
%     torque_density = -1;
%     return
% end
% if temperatures(5) > 130
%     %cost = -4;
%     
%     efficiency = -1;
%     torque_density = -1;
%     return
% end
% if temperatures(4) > 130
%     %cost = -4;
%     
%     efficiency = -1;
%     torque_density = -1;
%     return
% end
%mass of iron
mFes = m_Fey + m_Fez; % stator
mFer = (pi*(Dr/2)^2*l-pi*(Dri/2)^2*l)*rho_Fe; % rotor
mFe = mFes + mFer;
%mass of copper
mCu = lav*A_cond*Nph*m*rho_Cu;

% mass of permanet magnets
alphaMag = lp/((Dr/2+lm));
VMags = (alphaMag/2)*(Dr*lm+lm^2)*l;
mMagnets = 2*p*VMags*rhoPM;

cost = kLoss*P_loss + Cu_price*mCu + Fe_price*mFe + Magnet_price*mMagnets;
efficiency = eta2;
torque_density = T/(mFe+mMagnets+mCu);
mass = mFe+mMagnets+mCu;

objectives = [Pout, torque_density, mass, efficiency, cos_phi, cost];
% constraints
slot_constr = tau_u - 7e-3; % > 0
P_out_constr = Pout/Pel; % 1 <= Pout/Pel < 1.05
T_constr = temperatures(6); % < 100
constraints = [slot_constr P_out_constr T_constr];

end