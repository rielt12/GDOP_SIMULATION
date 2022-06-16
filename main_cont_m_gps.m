close all
clear all

lambda=3e8/1.57542003E9; %GPS L1 Frequency
longitude=34.7984651;
latitude=32.1163502; 
H=33; 

location = [latitude,longitude,H];

addpath(genpath('GPS_loc_computation'))
% We use heavens above to find all the satellites visible and we collect
% the doppler shift
% We collect the relevant shifts

data =xlsread("Test_Data_06_14_2022_14_59_41\raw2.xlsx");
shift(1:8,1) = data(1:8,5);
%shift(1:9,1) = data(1:9,5)





meas_ADR =[-269.41410048554342893112334422767, -517.876475118100643157958984375, 297.06787504255771636962890625, 119.83127115915218041664047632366, -422.42423600206774381149443797767, 378.67154149959486630905303172767, 634.413301013410091400146484375, 427.68595926463603973388671875];



svids = [21;
1
8
22
3
32
27
10
];
 
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


UTC_time = 1655218812.5175187587738037109375;
a=datetime(UTC_time,'convertfrom','epochtime','Epoch','1970-01-01');
[T, D, W, R, LS] = GPSdatetime(a,'LeapSecond',18);

p = get_satpos(T, svids(1),eph,2);
pos(1,:) = p(1:3,1);
p = get_satpos(T, svids(2),eph,2);
pos(2,:) = p(1:3,1);
p = get_satpos(T, svids(3),eph,2);
pos(3,:) = p(1:3,1);
p = get_satpos(T, svids(4),eph,2);
pos(4,:) = p(1:3,1);
p = get_satpos(T, svids(5),eph,2);
pos(5,:) = p(1:3,1);
p = get_satpos(T, svids(6),eph,2);
pos(6,:) = p(1:3,1);
p = get_satpos(T, svids(7),eph,2);
pos(7,:) = p(1:3,1);
p = get_satpos(T, svids(8),eph,2);
pos(8,:) = p(1:3,1);
% p = get_satpos(UTC_time, svids(9),eph,2);
% pos(9,:) = p(1:3,1);


del_t = 0.1;
v = ((get_satpos(T+del_t, svids(1),eph,2)-get_satpos(T-del_t, svids(1),eph,2))./(2*del_t))';
vel(1,:) = v(1,1:3);
v = ((get_satpos(T+del_t, svids(2),eph,2)-get_satpos(T-del_t, svids(2),eph,2))./(2*del_t))';
vel(2,:) = v(1,1:3);
v = ((get_satpos(T+del_t, svids(3),eph,2)-get_satpos(T-del_t, svids(3),eph,2))./(2*del_t))';
vel(3,:) = v(1,1:3);
v = ((get_satpos(T+del_t, svids(4),eph,2)-get_satpos(T-del_t, svids(4),eph,2))./(2*del_t))';
vel(4,:) = v(1,1:3);
v = ((get_satpos(T+del_t, svids(5),eph,2)-get_satpos(T-del_t, svids(5),eph,2))./(2*del_t))';
vel(5,:) = v(1,1:3);
v = ((get_satpos(T+del_t, svids(6),eph,2)-get_satpos(T-del_t, svids(6),eph,2))./(2*del_t))';
vel(6,:) = v(1,1:3);
v = ((get_satpos(T+del_t, svids(7),eph,2)-get_satpos(T-del_t, svids(7),eph,2))./(2*del_t))';
vel(7,:) = v(1,1:3);
v = ((get_satpos(T+del_t, svids(8),eph,2)-get_satpos(T-del_t, svids(8),eph,2))./(2*del_t))';
vel(8,:) = v(1,1:3);
% v = ((get_satpos(T+del_t, svids(9),eph,2)-get_satpos(UTC_time-del_t, svids(9),eph,2))./(2*del_t))';
% vel(9,:) = v(1,1:3);









% Goal to calculate GDOP and position
current_time = T;



%plot_cost(shift, pos, lambda, current_time, ephemeris, JD_prop_to)

[state ~, Error, delta_y_0,error_0,A_0, del_ADR] = doppler_shift_positioning_m_gps(shift, pos, lambda, current_time, eph, svids);
vpa(state)

figure(1000)
plot(Error)
xlabel('iterations')
ylabel('error')
title('error vs. iterations')

figure(700)
clf
num=1:8;
scatter(num,shift,'MarkerFaceColor','blue')
hold on
scatter(num,del_ADR,'MarkerFaceColor','red')
xlabel('satellite number')
ylabel('ADR Derivatives')
legend('measured','estimated')


meas_sim_err = norm(del_ADR'-shift);
