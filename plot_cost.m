function plot_cost(shift, pos, lambda, current_time, ephemeris, JD_prop_to)

close all
% inputs
% shifts are observed doppler shifts
% current time is the time in minutes since epoch (elapsed time) that all sats are in view.
% julian date is JD when sats are in view needed  to convert ECI to ECEF


%constants
c= 3e8;
longitude=-77.3790;
latitude=38.9649; 
H=0; 
%1085027.1960869177710264921188354
%-4845789.2537157507613301277160645
%3989288.0202742042019963264465332
[x,y,z]=lla2ecef_AB(latitude*(2*pi/360),longitude*(pi/180),H); % this function takes east longitude. change if need be.
rec_pos(1,1) =  x;
rec_pos(2,1) =  y;
rec_pos(3,1) =  z;
%rec_pos(1,1) = 0;
%rec_pos(2,1) = 0;
%rec_pos(3,1) = 0;
rec_clock_bias =0  ;
rec_vel = [0;0;0];
rec_clock_bias_rate= 0;

y_i = [rec_pos; rec_clock_bias; rec_vel; rec_clock_bias_rate];
y_0 =y_i;

iter = 1;
for g = 0.1e-3:1e-4:1000e-3
% rec_pos(1,1) = rec_pos(1,1)+g;
% X_change(iter,1)= rec_pos(1,1);
% rec_pos(iter,1) = rec_pos(2,1)+g;
% Y_change(iter,1)= rec_pos(2,1);
% rec_pos(3,1) = rec_pos(3,1)+g;
% Z_change(iter,1)= rec_pos(3,1);
% rec_clock_bias = rec_clock_bias + g;
% C_BIAS_change(iter,1) = rec_clock_bias;
% rec_vel(1,1) = rec_vel(1,1)+g;
% VX_change(iter,1)= rec_vel(1,1);
% rec_vel(2,1) = rec_vel(2,1)+g;
% VY_change(iter,1)= rec_vel(2,1);
% rec_vel(3,1) = rec_vel(3,1)+g;
% VZ_change(iter,1)= rec_vel(3,1);
  rec_clock_bias_rate = rec_clock_bias_rate + g;
 C_BIAS_rate_change(iter,1) = rec_clock_bias_rate;

   y_i = [rec_pos; rec_clock_bias; rec_vel; rec_clock_bias_rate];
Kost(iter,1) = compute_cost_doppler(y_i , shift, pos, lambda,ephemeris,JD_prop_to);
iter = iter + 1;
end

plot(C_BIAS_rate_change,Kost)


end