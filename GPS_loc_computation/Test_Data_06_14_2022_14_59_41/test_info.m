clear all

% June 14th 2022_06_14_14_59_41 local israel time
% day 165 of year
% in UTC it is 145941 pm 
% Week 2214
% GPS seconds of week around 226799
% LeapSeconds 18
% ephemeris time of day is 53999

%checking svid 21
%UTC time is 1655218812.5175187587738037109375
% day of week 2


addpath '..\'
addpath '..\..\'
UTC_time = 1655218812.5175187587738037109375;
%UTC_time  = 1655220010.518310070037841796875;

a=datetime(UTC_time,'convertfrom','epochtime','Epoch','1970-01-01');
[T, D, W, R, LS] = GPSdatetime(a,'LeapSecond',18);

[nav leapsec cnt] = ReadNav("brdc1650.22n");
%[nav leapsec cnt] =ReadNav("AMC400USA_R_20221651400_01H_GN.rnx")

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
 

svids = [21;
1
8
22
3
32
27
10
];


bias = [ 
106.56481186797157079126918688416
  106.106679178774356842041015625
108.28082670271470533407409675419
108.43689345071717866630933713168
106.06789046277521038064151071012
107.72690038631384368272847495973
106.26553490012906877382192760706
 105.1043247431520626378187444061];

for i =1:length(svids)
    svid = svids(i)
    satpos = get_satpos(T,svid,eph,2);
satpos1 = get_satpos(T+del_t,svid,eph,2);
satpos2 = get_satpos(T-del_t,svid,eph,2);
satvel = (satpos1-satpos2)./(2*del_t);
satvel = satvel(1:3,1);

% ramat aviv location
%32.1163502,34.7984651,33.000000


office_lat = 32.1163502; 
office_lon = 34.7984651;
office_alt = 33.000000;
[x_office,y_office,z_office] = lla2ecef_AB(office_lat*(2*pi/360),office_lon*(2*pi/360),office_alt);
v = [0;0;0];

r = [x_office;y_office;z_office];
pos = [satpos(1,1); satpos(2,1); satpos(3,1)];
del_R = 0;
del_R_rate = 0;
t_R = T;
lambda=3e8/1.57542003E9;

pos
unit_vec_gps = (pos -[x_office;y_office;z_office])./norm(pos -[x_office;y_office;z_office]);
pseudorange_rate = dot(satvel,unit_vec_gps)

del_ADR = acculumulated_delta_range_derivative_gps(r, del_R,t_R, v, del_R_rate, pos,lambda,eph,svid)-bias(i)
end
