% _______________________________________________________________________
% main.m
% version 1 (November, 2024)
% subroutines required: propolar.m, prospect_DB.m,
%                      calctav.m, dataSpec_PDB.m, data_external.m
% _______________________________________________________________________

% This script allows to simulate multi-angular photometric and polarimetric properties of leaf reflection (BRF, BPRF, DOLP) and
% to retrieves both leaf biochemical and surface structural parameters from BRF.
%
% Non-polarized component (400-2500 nm) with 1 nm step is calculated using PROSPECT-DB (Feret et al. 2017).
% The matlab codes of PROSPECT-DB are available on  website of
% http://teledetection.ipgp.jussieu.fr/prosail/.
%
% Polarized component is calculated based on a three-parameter function
% (linear coefficient, refractive index, and roughness of leaf surface).
% _______________________________________________________________________
%
% Xiao Li, Zhongqiu Sun, Shan Lu, Kenji Omasa,2024
% A radiative transfer model for characterizing photometric and polarimetric properties of leaf reflection: 
% combination of PROSPECT and a polarized reflection function
%
% Feret, J.-B., Gitelson, A.A., Noble, S.D., & Jacquemoud, S., 2017. 
% PROSPECT-D: towards modeling leaf optical properties through a complete lifecycle. 
%
% Bousquet, L., Lachérade, S., Jacquemoud, S., & Moya, I., 2005.
% A radiative transfer model for characterizing photometric and polarimetric properties of leaf reflection:
% combination of PROSPECT and a polarized reflection function
% 
% Litvinov, P., Hasekamp, O., & Cairns, B., 2011. 
% Models for surface reflection of radiance and polarized radiance: Comparison with airborne multi-angle 
% photopolarimetric measurements and implications for modeling top-of-atmosphere measurements. 
% _______________________________________________________________________

clear
clc
%% **************  Spectral simulation  ***********************

% load  structure and  biochemistry parameters
load('leaf_parameter.txt');%leaf parameter
% aerfa     = leaf_parameter(1);% linear coefficient
% rough    = leaf_parameter(2);% roughness of leaf surface
% n     = leaf_parameter(3);% refractive index factor
% N     = leaf_parameter(4);% leaf structure parameter
% Cab   = leaf_parameter(5);% chlorophyll a+b content in μg/cm2
% Car   = leaf_parameter(6);% carotenoids content in  μg/cm2
% Anth  = leaf_parameter(7); % Anthocyanin content in nmol/cm2
% Cbrown= leaf_parameter(8);% brown pigments content in arbitrary units
% Cw    = leaf_parameter(9);% equivalent water thickness in g/cm2 or cm
% Cm    = leaf_parameter(10);% dry matter content in g/cm2

% load illumination-viewing geometry
geo=load('geometry.txt');%leaf parameter
% SZA    = geo(1,:);% Source zenith angle,degree
% VZA    = geo(2,:);% Viewing zenith angle,degree
% VAA    = geo(3,:);% Viewing azimuth angle,degree
geometry=deg2rad(geo);% transfer degree into radian

% identification of wavelength
waveo=400:2500;
wave=waveo-399;% from 400 nm

% choose leaf optical properties
choose=999; %1 BRF; 2 BPRF;3 DOLP; 999 ALL

% simulate
opticalout=propolar(leaf_parameter,wave,geometry,choose);%simulate BRF
% outBRF=[geo;opticalout{1}]; % Correspond geometric parameters to BRF
% outBPRF=[geo;opticalout{2}]; % Correspond geometric parameters to BPRF
% outDOLP=[geo;opticalout{3}]; % Correspond geometric parameters to DOLP

% get simulated leaf optical properties
if(choose==1)
    outBRF=[geo;opticalout{1}]; % Correspond geometric parameters to BRF
elseif (choose==2)
    outBPRF=[geo;opticalout{2}]; % Correspond geometric parameters to BPRF
elseif(choose==3)
    outDOLP=[geo;opticalout{3}]; % Correspond geometric parameters to DOLP
elseif(choose==999)
    outBRF=[geo;opticalout{1}]; % Correspond geometric parameters to BRF
    outBPRF=[geo;opticalout{2}]; % Correspond geometric parameters to BPRF
    outDOLP=[geo;opticalout{3}]; % Correspond geometric parameters to DOLP
end

% output as excel file
if(choose==1)
    xlswrite('leaf_spectrum.xlsx',outBRF,'BRF')% output BRF
elseif(choose==2)
    xlswrite('leaf_spectrum.xlsx',outBPRF,'BPRF')% output BPRF
elseif(choose==3)
    xlswrite('leaf_spectrum.xlsx',outDOLP,'DOLP')% output DOLP
elseif(choose==999)
    xlswrite('leaf_spectrum.xlsx',outBRF,'BRF')% output BRF
    xlswrite('leaf_spectrum.xlsx',outBPRF,'BPRF')% output BPRF
    xlswrite('leaf_spectrum.xlsx',outDOLP,'DOLP')% output DOLP
end
%% *****************  Model inversion  *************************
choose=1;% invert based on BRF

% identify parameter boundary and initial value
% P0=[aerfa n rough N  Cab  Car  Anth  Cbrown  Cw  Cm ]
P0=[20  0.3 1.1  1.5	40	10	0.1	0	0.01	0.01];% initial value
LB=[0.1  0.01 1.001 1	 0	0	0	0	0	0]; %lower boundary
UB=[50  1  1.2 3.5	120	30	40	1	0.06	0.05];%upper boundary

options = optimset('Algorithm','trust-region-reflective','TolFun',1e-10);
sol= lsqcurvefit(@propolar,P0,wave,opticalout{1},LB,UB,options,geometry,choose);% invert based on BRF,geometry parameters input as prior information

% reconstruct spectral
choose=999;% reconstruct BRF, BPRF and DOLP
mnout=propolar(sol,wave,geometry,choose);% simulate specular and diffuse component with inverted parameters
mnBRF=mnout{1};% mnout{1}:simulated BRF
mnBPRF=mnout{2};% mnout{2}:simulated BPRF
mnDOLP=mnout{3};% mnout{1}:simulated DOLP

%  output
xlswrite('leaf_spectrum.xlsx',mnBRF,'Simulated_BRF')% output simulated BRF
xlswrite('leaf_spectrum.xlsx',mnBPRF,'Simulated_BPRF')% output simulated BRF
xlswrite('leaf_spectrum.xlsx',mnDOLP,'Simulated_DOLP')% output simulated BRF

