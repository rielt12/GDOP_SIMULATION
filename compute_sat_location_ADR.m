clear all

addpath '\\storm\users\arielb\Documents\GDOP_SIMULATION\GPS_loc_computation\'
addpath '\\storm\users\arielb\Documents\GDOP_SIMULATION\GPS_loc_computation\Test Data'

[nav leapsec cnt] = ReadNav("ABPO00MDG_R_20220861300_01H_GN.rnx");

for i=1:17
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






satpos = get_satpos(50400,5,eph,2);
vpa(satpos)