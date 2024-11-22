% _______________________________________________________________________
% prospecular.m (November, 2024)
%
% This script used to simulate multi-angular photometric and polarimetric properties of leaf reflection (BRF, BPRF, DOLP).
% _______________________________________________________________________
% ***********************************************************************
% Xiao Li, Zhongqiu Sun, Shan Lu, Kenji Omasa,2024
% A radiative transfer model for characterizing photometric and polarimetric properties of leaf reflection: 
% combination of PROSPECT and a polarized reflection function
% ***********************************************************************
% This script allows to simulate multi-angular photometric and polarimetric properties based on the PROSPECT-DB model (Féret et al. 2021)
% and a three-parameter polarized function (linear coefficient, refractive index, and roughness of leaf surface) (Bousquet et al. 2005)

% ***********************************************************************
% Leaf optical properties are calculated from 400 nm to 2500 nm (1 nm step) 

function out=propolar(P0,wave,xr,choose)
% P0: surface structure and biochemical parameters of PROSPECULAR model;
%       - aerfa = linear coefficient 
%       - rough = roughness of leaf surface
%       - n   = refractive index factor
%       - N   = leaf structure parameter
%       - Cab = chlorophyll a+b content in μg/cm2
%       - Car = carotenoids content in  μg/cm2
%       - Anth = Anthocyanin content in nmol/cm2
%       - Cbrown= brown pigments content in arbitrary units
%       - Cw  = equivalent water thickness in g/cm2 or cm
%       - Cm  = dry matter content in g/cm2
%
% xr: source zenith angle, viewing zenith angle and viewing azimuth angle, radian;
% wave: wavelenth (nm)


%%  Non-polarized component calculation
N=P0(4);
Cab=P0(5);
Ccx=P0(6);
Can=P0(7);
Cbp=P0(8);
Cw=P0(9);
Cdm=P0(10);
RT=prospect_DB(N,Cab,Ccx,Can,Cbp,Cw,Cdm);
R=RT(wave,2); %non-polarized component

%%  Polarized component caliculation
A= P0(1); %α linear coefiicient factor
cu= P0(2); % roughness

% refractive index of leaf surface of 2009 Stucken
data    = data_external;
waven = data(wave,2);
n= P0(3).*waven;% leaf surfafce refractive index, wavelength-dependent

% illumination-viewing geometry
zenithi = xr(1,:);% source zenith angle, radian
zenithr = xr(2,:);% viewing zenith angle, radian
azimuthi = zeros(size(zenithi,1),size(zenithi,2));% source azimuth angle, radian
azimuthr = xr(3,:);% viewing azimuth angle, radian

czi=cos(zenithi);%cosθi
szi=sin(zenithi);%sinθi
czr=cos(zenithr);%cosθr
szr=sin(zenithr);%sinθr
fir=cos(abs(azimuthr-azimuthi));%cos δφ
sir=sin(azimuthr-azimuthi);% sin δφ

aerfa=0.5.*acos(czi.*czr+szi.*szr.*fir);% half of the phase angle
comg=-cos(2.*aerfa);%cos(scattering angle).There is a complementary relationship between scattering angle and phase angle
aerfat=asin(sin(aerfa)./n);% refractive angle
beita=acos((czr+czi)./(2.*cos(aerfa)));%
caerfa=cos(aerfa);%cos half of the phase angle
caerfat=cos(aerfat);%cos refractive angle

% a three-parameter polarized function (linear coefficient, refractive index, and roughness of leaf surface)
fm=8.*pi.*czi.*czr.*cos(beita);%
fen1=(n.*caerfat-caerfa)./(n.*caerfat+caerfa);
fen2=(n.*caerfa-caerfat)./(n.*caerfa+caerfat);
Fp=0.5.*(fen1.^2-fen2.^2);
cu2=2.*cu.^2;
b1=(A.*pi.*Fp)./(4.*cos(beita).*(czi+czr));
b21=exp(-(tan(beita).^2)./cu2);
b22=pi.*(cos(beita).^3).*cu2;
b2=b21./b22;
BPDF=b1.*b2;
BPRF=BPDF.*pi;

% Fresnel correction coefficient 
sf11=aerfat-aerfa;
sf12=aerfat+aerfa;
sf21=sin(sf11).^2;
sf22=sin(sf12).^2;
sf23=tan(sf11).^2;
sf24=tan(sf12).^2;
sfp=(sf21.*sf24-sf22.*sf23)./(sf21.*sf24+sf22.*sf23);

BPRF=BPRF./sfp;

if(choose~=2)
BRF=R+BPRF;% BRF=non-polarized component+polarized component
DOLP=BPRF./BRF; %DOLP=polarized component /non-polarized component+polarized component

if(choose==1)% output iprf
    out=BRF;
elseif(choose==2)% output bprf
    out=BPRF;
elseif(choose==3)% output dolp
    out=DOLP;
elseif(choose==999)
    out={BRF,BPRF,DOLP, R};% output all the optical prooerties
end
end