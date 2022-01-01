clear all
close all
% Goal of simulation is to take constellation TLEs and:
% 1) simulate doppler shifts from them
% 2) calculate a navigation fix from the shifts
% 3) calculate GDOP


% constant decleration
R_e=6378.1; % Equatorial radius
f=0.00335; % Oblatness of the earth
lambda = 3e8/1626.5e6; % wavelength iridium
lambda = 3e8/11.325e9; % wavelength starlink

% We first propagate a TLE

%% 
% example: We are trying to gett the doppler shift at Jun 26, 05:50:05am (local time), 9:50:05 GMT/UTC with
% the Iridium 70 satellite. At latitude: 38.9649°N and longitude:
% 77.3790°W 
longitude=77.3790;
latitude=38.9649; 
H=0; 
location = [latitude,-1*longitude,H];



% get TLE info
fid = fopen("starlink_set.TLE");
n = linecount(fid);
fclose(fid);

fopen("starlink_set.TLE");
for i=1:3:n-3
tline = fgetl(fid);
longstr1=char(tline);
tline=fgetl(fid);
longstr2=char(tline);
tline=fgetl(fid);
longstr3=char(tline);
sat_name = convertCharsToStrings(longstr1);
sat_name = strrep(sat_name, ' ', '');
sat_name = strrep(sat_name, '-', '');
sat_name = strrep(sat_name, '(', '');
sat_name = strrep(sat_name, ')', '');
satrec.(sat_name)= twoline2rvMOD(longstr2,longstr3);
end


fNames = fieldnames(satrec);
for n = 1:length(fNames)
    disp(fNames(n))
end

% we start on the 26th of december 2021 at 0:00 and propagate for 24 hours
% minute by minute.

Y=2022; 
M=01; 
D=01; 
HR=16; 
MN=01; 
S=03;
JD_init= jday(Y,M,D,HR,MN,S);

end_t =10;
for t = 1:end_t
    t
JD_prop_to= (JD_init*24*60*60 + t)/(24*60*60);
for n = 1:length(fNames)
    elapsedtime=(JD_prop_to-satrec.(fNames{n}).jdsatepoch)*24*60; % minutes
    [satrec.(fNames{n}), r, v] = sgp4(satrec.(fNames{n}),elapsedtime);
    JD_prop_to_dt =datetime(JD_init,'convertfrom','juliandate');
    [y,mo,d] = ymd(JD_prop_to_dt);
    [h,m,s] = hms(JD_prop_to_dt);
    time  =[y,mo,d,h,m,s];
    %elevation calc
    aer = eci2aer(r.*1000,time, location) ;
    satrec.(fNames{n}).elevation(t,1) = aer(1,2);
    
    % doppler shift calc
    [x,y,z]=lla2ecef_AB(latitude*(2*pi/360),(2*pi)-(longitude*(pi/180)),H); % this function takes east longitude. change if need be.
    ecef_ref=[x/1000;y/1000;z/1000];
    eci_ref=ECEFtoECI(JD_prop_to,ecef_ref); % reference (ground station) location in ECI
    eci_ref=transpose(eci_ref);
    r1_topo=r-eci_ref; % topocentric location of satellite in ECI
    % we project the eci velocity onto the eci position
    v_c = dot(v,r1_topo./norm(r1_topo));% closing velocity
    satrec.(fNames{n}).doppler_shift(t,1) = -norm(v_c)/lambda;
        
end

end

% 
% for n = 1:length(fNames)
%     figure(n)
%     plot(1:end_t, satrec.(fNames{n}).doppler_shift)
%     legend(fNames{n})
%     xlabel('time (seconds)')
%     ylabel('dopper shift')
%     title('doppler shift vs. time')
%     
%     figure(n+1000)
%     plot(1:end_t, satrec.(fNames{n}).elevation)
%     legend(fNames{n})
%     xlabel('time (seconds)')
%     ylabel('elevation')
%     title('elevation vs. time')
% end

save('starlink_ephemeris_altitude')

