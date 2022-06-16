clear all

% March 27th 2022_03_27_16_34_55 local israel time
% day 86 of year
% in UTC it is 133455 pm 
% Week 2203
% GPS seconds of week around 48913
% LeapSeconds 18
% ephemeris time is 50400


% checking SV5
% UTC time is 1648388118.71606731414794921875
addpath '..\'
addpath '..\..\'

a=datetime(1648388118.71606731414794921875,'convertfrom','epochtime','Epoch','1970-01-01');
[T, D, W, R, LS] = GPSdatetime(a,'LeapSecond',18);



[nav leapsec cnt] = ReadNav("ABPO00MDG_R_20220860000_01D_GN.rnx");
%[nav leapsec cnt] = ReadNav("ABPO00MDG_R_20220861300_01H_GN.rnx");
%[nav leapsec cnt] = ReadNav("ABPO00MDG_R_20220861400_01H_GN.rnx");

for i=1:cnt
eph(1,i)= nav(i).prn;
eph(2,i)= nav(i).afc(3,1);
eph(3,i) = nav(i).M0;
eph(4,i) = nav(i).art;
eph(5,i) = nav(i).dn;
eph(6,i) = nav(i).ecc;
eph(7,i) = nav(i).w;
eph(8,i) = nav(i).cuc;
eph(9,i) = nav(i).cus;
eph(10,i) = nav(i).crc;
eph(11,i) =nav(i).crs;
eph(12,i) = nav(i).i0;
eph(13,i) = nav(i).idt;
eph(14,i) = nav(i).cic;
eph(15,i) = nav(i).cis;
eph(16,i) = nav(i).om0;
eph(17,i) = nav(i).omd;
eph(18,i) = nav(i).toe;
eph(19,i) = nav(i).afc(1,1);
eph(20,1) = nav(i).afc(2,1);
eph(21,1) = nav(i).transmit;
end

del_t = 0.1;
svid = 31;
satpos = get_satpos(T,svid,eph,2);
satpos1 = get_satpos(T+del_t,svid,eph,2);
satpos2 = get_satpos(T-del_t,svid,eph,2);
satvel = (satpos1-satpos2)./(2*del_t);
satvel = satvel(1:3,1);







% Aeronautics location
office_lat = 31.900487; 
office_lon = 34.7360919;
office_alt = 11;
[x_office,y_office,z_office] = lla2ecef_AB(office_lat*(2*pi/360),office_lon*(2*pi/360),office_alt);
v = [0;0;0];

r = [x_office;y_office;z_office];
pos = [satpos(1,1); satpos(2,1); satpos(3,1)];
del_R = 0;
del_R_rate = 0;
t_R = T;
lambda=3e8/1.57542003E9;




unit_vec_gps = (pos -[x_office;y_office;z_office])./norm(pos -[x_office;y_office;z_office]);
pseudorange_rate = dot(satvel,unit_vec_gps)

del_ADR = acculumulated_delta_range_derivative_gps(r, del_R,t_R, v, del_R_rate, pos,lambda,eph,svid)
