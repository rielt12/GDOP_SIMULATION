function exist=predict(fname,latitude,longitude,H,Y,M,D,HR,MN,S,Mirror)

%%
% inputs:
% fname is the name of the tle file
% longitude is the west longitude of the ground station 
% latitude is the north latitude of the ground station
% H is the altitude
% Y is the year of the predicted flare
% M is the month of the predicted flare
% D is the day of the predicted flare
% HR is the hour of the predicted flare
% MN is the minute of the predicted flare
% S is the second of the predicted flare
% Mirror is which mirror is causing the flare 
% Input 1 for the front mirror
% Input 2 for the right mirror 
% Input 3 for the left mirror.


%%
% We are trying to confirm the prediction that there will be an iridium
% flare on Y,M,D,HR,MN,S with the the time being in UTC and with the ground
% station at latitude and longitude and H.
%%


%% 
% example: We are trying to confirm the prediction that there will be an
% iridium flare on Jun 26, 05:50:05am (local time), 9:50:05 GMT/UTC with
% the Iridium 70 satellite. At latitude: 38.9649°N and longitude:
% 77.3790°W 
% longitude=77.3790;
% latitude=38.9649; 
% H=0; 
% altitude=0; 
% Y=2015; 
% M=6; 
% D=26; 
% HR=9; 
% MN=50; 
% S=05;





R_e=6378.1; % Equatorial radius
f=0.00335; % Oblatness of the earth




%%
% We figure out where the sun is at the predicted time.
%%
 
JD2=jday(Y,M,D,HR,MN,S); % We get the Julian date of the predicted flare
global suncoef
suncoef=1;
[rasc, decl, rsun] = sun2 (JD2);
[SUN]=rsun;
 

%%
% We need to find the position of the satellite in ECI coordinates.
% We do this with the SGP4 algorithm.
%%


fid = fopen(fname);
tline = fgets(fid);
longstr1=char(tline);
tline=fgetl(fid);
longstr2=char(tline);




satrec = twoline2rvMOD(longstr1,longstr2);

JD1=satrec.jdsatepoch; % Julian Date of epoch of TLE
JD1=double(vpa(JD1));
elapsedtime=(JD2-JD1)*24*60;



[satrec, r, v] = sgp4(satrec,elapsedtime);



%angular_difference=iridium_flare(EQ_X,EQ_Y,EQ_Z,V_X, V_Y, V_Z,TOPO_X,TOPO_Y,TOPO_Z,SUN_ALPHA, SUN_DELTA,ROT_1, ROT_2)
EQ_X=r(1,1);
EQ_Y=r(1,2);
EQ_Z=r(1,3);

V_X=v(1,1);
V_Y=v(1,2);
V_Z=v(1,3);





% convert to topocentric coordinates
    

[x,y,z]=lla2ecef_AB(latitude*(2*pi/360),(2*pi)-(longitude*(pi/180)),H); % this function takes east longitude. change if need be.
ecef_ref=[x/1000;y/1000;z/1000];
eci_ref=ECEFtoECI(JD2,ecef_ref); % reference (ground station) location in ECI
eci_ref=transpose(eci_ref);


[r1]=r-eci_ref; % topocentric location of satellite in ECI





 TOPO_X=r1(1,1);
 TOPO_Y=r1(1,2);
 TOPO_Z=r1(1,3);

 
 if Mirror==1
     ROT_2=0;
 elseif Mirror==2
     ROT_2=120;
 elseif Mirror==3
     ROT_2=-120;
 end
 
 
[angular_difference]=iridium_flare(EQ_X,EQ_Y,EQ_Z,V_X, V_Y, V_Z,TOPO_X,TOPO_Y,TOPO_Z,SUN,-40, ROT_2);


if angular_difference < 2
    exist=1;
else
    exist=0;
end 

end