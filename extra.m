% code to compute azimuth and elevation

[x,y,z]=lla2ecef_AB(latitude*(2*pi/360),-(longitude*(pi/180)),H); % this function takes east longitude. change if need be.
ecef_ref=[x/1000;y/1000;z/1000];
eci_ref=ECEFtoECI(JD2,ecef_ref); % reference (ground station) location in ECI
eci_ref=transpose(eci_ref);
GMST = JD2GMST(JD2);
LST = GMST-longitude;
[range,az,el]=eci2rangeazel(r.*1000,eci_ref.*1000 ,LST, latitude);

S.a = 'test';
S.b = 1;
S.c = rand(2);
fNames = fieldnames(S);
for n = 1:length(fNames)
    disp(S.(fNames{n}))
end