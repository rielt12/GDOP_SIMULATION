clear all
% Goal of simulation is to take constellation TLEs and:
% 1) simulate doppler shifts from them
% 2) calculate a navigation fix from the shifts
% 3) calculate GDOP


% constant decleration
R_e=6378.1; % Equatorial radius
f=0.00335; % Oblatness of the earth
lambda = 3e8/1626.5e6; % wavelength iridium

% We first propagate a TLE

%% 
% example: We are trying to gett the doppler shift at Jun 26, 05:50:05am (local time), 9:50:05 GMT/UTC with
% the Iridium 70 satellite. At latitude: 38.9649°N and longitude:
% 77.3790°W 
longitude=77.3790;
latitude=38.9649; 
H=0; 


% get TLE info
fid = fopen("test_tle.TLE");
tline = fgets(fid);
longstr1=char(tline);
tline=fgetl(fid);
longstr2=char(tline);

satrec = twoline2rvMOD(longstr1,longstr2);

JD1=satrec.jdsatepoch; % Julian Date of epoch of TLE
JD1=double(vpa(JD1));



Y=2015; 
M=6; 
D=26; 
HR=9; 
MN=50; 
S=05;


JD2=jday(Y,M,D,HR,MN,S); % We get the Julian date of doppler shift
elapsedtime=(JD2-JD1)*24*60; % minutes

doppler_shift =[];
elevation = [];

i = 1;
for prop_time = elapsedtime:elapsedtime+60
[satrec, r, v] = sgp4(satrec,prop_time);
time = [2015; 6; 26; 9; 50; 5]';
location = [latitude,-longitude,H];
aer = eci2aer(r.*1000,time, location);
elevation(i,1) = aer(1,2);
%if(elevation >0)
[x,y,z]=lla2ecef_AB(latitude*(2*pi/360),(2*pi)-(longitude*(pi/180)),H); % this function takes east longitude. change if need be.
ecef_ref=[x/1000;y/1000;z/1000];
eci_ref=ECEFtoECI(JD2,ecef_ref); % reference (ground station) location in ECI
eci_ref=transpose(eci_ref);


r1_topo=r-eci_ref; % topocentric location of satellite in ECI
% we project the eci velocity onto the eci position
v_c = dot(v,r1_topo./norm(r1_topo));% closing velocity
doppler_shift(i,1) = -norm(v_c)/lambda;
%end
i= i+1;
end





figure(1)
plot(doppler_shift)
figure(2)
plot(elevation)





