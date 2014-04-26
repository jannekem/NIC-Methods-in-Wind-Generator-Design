bz = b2+b3/2;
%% Thermal resistances of different parts of the machine

% FRAME AND YOKE
lambda_Al = 209;% thermal conductivity of aluminium (electrotechnical)
hfr = h_H - Dse/2; % frame height = ??
Rthfr = hfr/(pi*l*lambda_Al*(Dse+hfr));

lambda_Fe = 74.7;% thermal conductivity of iron (pure)
Rthy = (log(Dse/2)-log(D/2+h_tooth))/(2*pi*l*fr*lambda_Fe);

% lambda_g = 0.025;%thermal conductivity of the gap
% ge = 40e-3; %equivalent air gap
% Ag = 1; %equivalent gap area = ????
% Rthcy = ge/(lambda_g*Ag); % contact resistance

% Approximation for the contact resistance
Rthcy = Rthy; % approximation


% STATOR TEETH

% approximate thermal resistance of the teeth, the real resistance is
% obtained by integrating because the area of the tooth is not constant
Rthz = h_tooth/(lambda_Fe*Q*fr*l*bz);

% INTERNAL AIR
% between internal air and frame
vr = Omega*Dr/2;
alpha1 = 15 + 6.75^0.65 + vr^0.65;
Aes = pi*(Dse/2)^2; %internal end shield area
Aif = Dse*pi*l; % internal frame area
A1 = 2*Aes + Aif;
R1 = 1/(alpha1*A1);

% BETWEEN INTERNAL AIR AND ROTOR
r_delta = D/2;
alpha2 = 16.5^0.65*vr^0.65;
bfin = 0; % jäähdytysrivat
hfin = 0;
nfin = 0;
A2 = 2*bfin*hfin*nfin + pi*r_delta^2; 
R2 = 1/(alpha2*A2);

% between internal air and end windings
lov = lav - l;
d1 = D+2*h_tooth; % outer diameter of the endwinding wreath
d2 = D;% inner diameter of the end winding wreath
Aew = pi*lov*(d1+d2)/2;
alpha3 = 6.5 + 5.25^0.65*vr^0.6;
R3 = 1/(alpha3*Aew);

% WINDING
% per unit-length resistances
% b_eqv and h_eqv are equivalent rectangular slot dimensions
dI = 0.001; % slot insulation thickness
dA = 0.001; % air film thickness
d = dI + dA; % equivalent insulation thickness
% b_eqv = (bs2 + bs1)/2 - 2*d;
% h_eqv = 2*Su/(bs1 + bs2) - 2*d;
b_eqv = bs2 - 2*d;
h_eqv = h_tooth - 2*d;


%lambdaQ = 0.2*kCu; % thermal conductivity of the slot 
lambdaQ = 3.5*0.2; % thermal conductivity of the slot 
lambdaA = 0.025; %thermal conductivity of air (stagnant)
lambdaI = 0.2; % thermal conductivity of insulation (typical insulation system)
Rx0 = bs2/(h_eqv*lambdaQ);
Ry0 = h_tooth/(b_eqv*lambdaQ);
Rix = dI/(h_eqv*lambdaI) + dA/(h_eqv*lambdaA);
Riy = dI/(b_eqv*lambdaI) + dA/(h_eqv*lambdaA);

Rx = (1/2)*(Rix + Rx0/6); 
Ry = (1/2)*(Riy + Ry0/6);

Rth4 = Rx*Ry/(Q*l*fr*(Rx+Ry))*(1- Rx0*Ry0/(720*(Rx0+Ry0)));

% between coil side and endwinding
l_av = lav/2; % average conductor length of a half turn
A_Cu = A_cond; % cross sectional area of a conductor
lambda_Cu = 394; % thermal conductivity of copper (electrotechnical)
Rthw = l_av/(3*A_Cu*lambda_Cu);
Rth5 = Rthw/(2*Q);

% AIR-GAP
r_delta = D/2;
% kinematic viscosity of air, actually this is temperature dependent, this
% value is for about 60 degrees Celcius. It is just a guess for testing.
T_airgap = 60;
nuA = kinematic_viscosity_of_air(T_airgap);
Tam = (Omega^2 * r_delta * delta^3)/nuA^2;
%Nusselt number according to Lindström, (according to Pyrhönen, Jokinen and
%Hrabovcova, the equation is slightly different)
if Tam < 1740
    Nu = 2;
else
    Nu = 0.409*Tam^0.241 - 137*Tam^-0.75; 
end

alpha_delta = (Nu*lambdaA/(2*delta));
Rthdelta = 1/(alpha_delta*2*pi*r_delta*l);

% ROTOR YOKE

Rthyr = (log(Dr/2) - log(Dri/2))/(2*pi*fr*l*lambda_Fe);

% PERMANENT MAGNETS
lambda_PM = 9;
%RthPM = (log(Dr/2 + lm) - log(Dr/2))/(2*pi*l*lambda_Fe)*(2*pi)/(2*p*lp);
RthPM = (log(Dr/2 + lm) - log(Dr/2))/(2*pi*l*lambda_PM)*(tau_p/lp);

% RETAINING SLEEVE
Rthsl = 0;
% insulation between the magnets and the rotor 
dIr = 0.5e-3;
lambda_Ir = 0.2; %thermal conductivity of the insulation (typical insulation system)
%gamma_m = 
%Rthins = dIr/(pi*Dr*l*lambda_Ir)*(2*pi/(2*p*lp));
Rthins = dIr/(pi*Dr*l*lambda_Ir)*(tau_p/lp);

Rth9 = Rthdelta + (1/2)*Rthz + Rthsl + (1/2)*RthPM;

% SHAFT
lambda_sh = 45; % thermal conductivity of the shaft (steel, carbon steel 0.5 %)

lbb = 1.4*l;
Rthsh = lbb/(pi*(Dri/2)^2*lambda_sh);

% BEARING
% The thermal resistance of the bearings depends of the rotation speed.
% This equation is actually valid for rotation speeds that are close to the ones
% used in some experiments made by Kylander. However,
% Lindström has used this equation and neclegted the rotation speed
% dependence. 


k1 = 0.45;
k2 = 1;
k3 = 1;
D_Br_ave = Dri+Dri*0.5; % average diameter of bearing

kR = lambda_Al;
RR1 = 1/(4*pi*kR*0.05)*(1-(Dri/2)^2*log(D_Br_ave/Dri)/((D_Br_ave/2)^2 - (Dri/2)^2)); 
RR2 = 1/(4*pi*kR*0.05)*((D_Br_ave/2)^2*log(D_Br_ave/Dri)/((D_Br_ave/2)^2 - (Dri/2)^2) - 1);
Rthb = RR1 + RR2;
Rth11 = (1/4)*Rthb;

% BETWEEN THE SHAFT AND THE ROTOR CORE
Rther = Rthyr;

Rth10 = (1/2)*RthPM + Rthins + Rthyr + Rther + (1/2)*Rthsh + (1/4)*Rthb;

Rth1 = ((1/2)*Rthfr);
%% Thermal conductances of the network

Gth1 = 1/((1/2)*Rthfr);
Gth2 = 1/((1/2)*Rthfr + (1/2)*Rthy + (1/2)*Rthcy);
Gth3 = 1/((1/2)*(Rthy + Rthz));
Gth4 = 1/Rth4;
Gth5 = 1/(Rthw/(2*Q));
Gth6 = 1/(R1 + R2 + R1*R2/R3);
Gth7 = 1/(R2 + R3 + R2*R3/R1);
Gth8 = 1/(R1 + R3 + R1*R3/R2);
Gth9 = 1/Rth9;
Gth10 = 1/Rth10;
Gth11 = 1/Rth11;



Gth = [Gth1+Gth2+Gth11+Gth6+Gth8 -Gth2 0 0 -Gth8 -Gth6 -Gth11;
    -Gth2 Gth2+Gth3 -Gth3 0 0 0 0;
    0 -Gth3 Gth3+Gth4+Gth9 -Gth4 0 -Gth9 0;
    0 0 -Gth4 Gth4+Gth5 -Gth5 0 0;
    -Gth8 0 0 -Gth5 Gth5+Gth7+Gth8 -Gth7 0;
    -Gth6 0 -Gth9 0 -Gth7 Gth7+Gth6+Gth10+Gth9 -Gth10;
    -Gth11 0 0 0 0 -Gth10 Gth10+Gth11];

%% Ambient temperature


T0 = 17 + 273.15;
P_loss_vec = [T0/Rth1 P_y P_z P_w P_ew P_ec P_Br]';

Theta = Gth\P_loss_vec;

  