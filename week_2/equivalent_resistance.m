% File created to generate the transistor equivalent resistance 
% Written on Sep. 20, 2015 by S. Rakheja
% For ECE 6473, Fall 2015

clear all
clc
%% define constants for the technology [Check chapter 3 from Jan M. Rabaey]

Vdd = 0.5:0.01:2.5; % supply voltage for the technology

% NFET process parameters
kn = 115e-6; % [A/V^2] product of mu_n*Cox 
VTn0 = 0.4; % [V] threshold voltage 
gamma_n = 0.4; %[sqrt(V)] body-effect coefficient 
lambda_n = 0.06; % [1/V] channel-length modulation 
% PFET process parameters
kp = 30e-6; % [A/V^2] product of mu_p*Cox 
VTp0 = -0.4; % [V] threshold voltage 
gamma_p = -0.4; %[sqrt(V)] body-effect coefficient 
lambda_p = -0.1; % [1/V] channel-length modulation

phif = 0.3; % Fermi potential in [V]. Example 3.5, Chapter 3 from Jan M. Rabaey

%% define design parameters
WL_ratio_NFET = 1; % ratio W/L for NFET
WL_ratio_PFET = 1; % ratio W/L for PFET

%% I-V equations (unified model)
%% NFET equations
% define voltages for NFET
Vgs = Vdd;
Vds = Vdd;
Vbs = 0;

VTn = VTn0+gamma_n*(sqrt(2*phif-Vbs)-sqrt(2*phif));
Vmin = min(Vds, Vgs-VTn); 
if (Vgs < VTn) 
    Id_nfet = 0;
else
    Id_nfet = kn*WL_ratio_NFET*((Vgs-VTn).*Vmin-Vmin.^2/2).*(1+lambda_n*Vmin);
end
Idsatn = 1/2*kn*WL_ratio_NFET*(Vgs-VTn).^2;
Reqn = 3/4*Vdd.*(1-5/6*lambda_n*Vdd)./Idsatn; % equivalent resistance when drain goes from Vdd to Vdd/2 
% assumption is that the transistor remains in saturation


% PFET equations
% define voltages for PFET
Vsg = Vdd;
Vsd = Vdd;
Vsb = 0;
VTp =  VTp0+gamma_p*(sqrt(2*phif-Vsb)-sqrt(2*phif));
Vmin = min(Vsd, Vsg-abs(VTp));

if (Vsg < abs(VTp))
    Id_pfet = 0;
else
    Id_pfet = kp*WL_ratio_PFET*((Vsg-abs(VTp)).*Vmin-Vmin.^2/2).*(1-lambda_p*Vmin);
end
Idsatp = 1/2*kp*WL_ratio_PFET*(Vsg-abs(VTp)).^2;
Reqp = 3/4*Vdd.*(1-5/6*abs(lambda_p)*Vdd)./Idsatp; % equivalent resistance when drain goes from 0 to Vdd/2
% assumption is that the transistor remains in saturation

plot(Vdd, Reqn, Vdd, Reqp)






