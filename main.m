clear all
close all
% Goal of simulation is to take constellation TLEs and:
% 1) simulate doppler shifts from them
% 2) calculate a navigation fix from the shifts
% 3) calculate GDOP


% constant decleration
R_e=6378.1; % Equatorial radius
f=0.00335; % Oblatness of the earth
c=3e8;
w_e = 7.2921150e-5; 
lambda = 3e8/1626.5e6; % wavelength iridium
lambda = 3e8/11.325e9; % wavelength starlink

% We first propagate a TLE

%% 
% example: We are trying to gett the doppler shift at Jun 26, 05:50:05am (local time), 9:50:05 GMT/UTC with
% the Iridium 70 satellite. At latitude: 38.9649°N and longitude:
% 77.3790°W 
longitude=-77.3790;
latitude=38.9649; 
H=0; 
location = [latitude,longitude,H];



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

% we start on the 1 of January 2022 at 18:01.13 (16:01:13 in UTC) and propagate for 24 hours
% second by second.

Y=2022; 
M=01; 
D=01; 
HR=16; 
MN=01; 
S=12;
JD_init= jday(Y,M,D,HR,MN,S);

end_t =1;
for t = 1:end_t
    t
JD_prop_to= (JD_init*24*60*60 + t)/(24*60*60);
for n = 1:length(fNames)
    elapsedtime=(JD_prop_to-satrec.(fNames{n}).jdsatepoch)*24*60; % minutes
    [satrec.(fNames{n}), r, v] = sgp4(satrec.(fNames{n}),elapsedtime);
    r= r.*1000; % convert to meters
    v = v.*1000; % convert to meters
    JD_prop_to_dt =datetime(JD_prop_to,'convertfrom','juliandate');
    [y,mo,d] = ymd(JD_prop_to_dt);
    [h,m,s] = hms(JD_prop_to_dt);
    time  =[y,mo,d,h,m,s];
    %elevation calc
    aer = eci2aer(r,time, location); 
    satrec.(fNames{n}).elevation(t,1) = aer(1,2);
    
    %convert to ecef
   [r_ecef,v_ecef]=eci2ecef(r',v',time(1),time(2),time(3),time(4),time(5),time(6)); %meters
    %[az,elev,slantRange]  =
    %ecef2aer(r_ecef(1),r_ecef(2),r_ecef(3),location(1), location(2), location(3), wgs84Ellipsoid, 'degrees');

    
    % doppler shift calc
    [x,y,z]=lla2ecef_AB(latitude*(2*pi/360),longitude*(pi/180),H); % this function takes east longitude. change if need be.
    ecef_ref=[x;y;z]; % meters
    eci_ref=ECEFtoECI(JD_prop_to,ecef_ref); % reference (ground station) location in ECI % meters
    eci_ref=transpose(eci_ref); % meters
    r1_topo=(r-eci_ref); % topocentric location of satellite in ECI meters
    r1_topo_ecef=(r_ecef-ecef_ref); % topocentric location of satellite in ECEF (meters)
   
    % we project the eci velocity onto the eci position
    v_c = dot(v,r1_topo./norm(r1_topo));% closing velocity
    
    % we project the ecef velocity onto the ecef position
    v_c_ecef =dot(v_ecef, r1_topo_ecef./norm(r1_topo_ecef));
   
%     satrec.(fNames{n}).doppler_shift(t,1) = -norm(v_c_ecef)/lambda;
%     satrec.(fNames{n}).pos(t,:) = r_ecef'; %meters
%     satrec.(fNames{n}).vel(t,:) = v_ecef'; % meters

    % equation (7) attempt
    % we assume a non stationary reciever 
    v_rec = [0;0;0];
    % we simulate 0 clock bias in the satellites
    del_j = 0; % seconds
    % we simulate 0 clock bias rates in the satellites
    del_j_dot = 0; % seconds/second
    % we simulate 0 reciever clock bias
    del_r = 0; % second0s
    % we simulate 0 reciever clock bias rase
    del_r_dot =0; % seconds/second

    t_prop = norm(r1_topo_ecef)/c; % approximation (?)
    [satrec.(fNames{n}), r_prop_eci, v_prop_eci] = sgp4(satrec.(fNames{n}),elapsedtime-del_r/60-t_prop/60);
    r_prop_eci = r_prop_eci.*1000;
    v_prop_eci = v_prop_eci.*1000; % convert to meters
    
    %convert to ecef
    [r_prop_ecef,v_prop_ecef]=eci2ecef(r_prop_eci',v_prop_eci',time(1),time(2),time(3),time(4),time(5),time(6));
    v_j = earth_rot_mat(t_prop)*v_prop_ecef;
    rho_j = (ecef_ref - earth_rot_mat(t_prop)*r_prop_ecef)./norm(ecef_ref - earth_rot_mat(t_prop)*r_prop_ecef);
    alpha_j = dot(-rho_j, cross([0;0;w_e],earth_rot_mat(t_prop)*r_prop_ecef));
% 
%     satrec.(fNames{n}).doppler_shift(t,1)=((dot(rho_j,v_rec-v_j)*((1+del_j_dot)/(1+(1/c)*(alpha_j-dot(rho_j,v_j))))+c*del_r_dot-c*del_j_dot)/(1+del_r_dot))/-lambda;
%     satrec.(fNames{n}).pos(t,:) = r_ecef'; %meters
%     satrec.(fNames{n}).vel(t,:) = v_ecef'; % meters
    
    %equation 10 attempt
    satrec.(fNames{n}).doppler_shift(t,1)= (acculumulated_delta_range_derivative(ecef_ref, del_r,elapsedtime, v_rec, del_r_dot, r_ecef,lambda,satrec.(fNames{n}), JD_prop_to));
    satrec.(fNames{n}).pos(t,:) = r_ecef'; %meters
    satrec.(fNames{n}).vel(t,:) = v_ecef'; % meters




    



   
  
    
        
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


