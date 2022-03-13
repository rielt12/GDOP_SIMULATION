function [state rec_GDOP, Error, delta_y_0,error_0,A_0] = doppler_shift_positioning_m(shift, pos, lambda, current_time, ephemeris, JD_prop_to)


max_iter = 1500;

% inputs
% shifts are observed doppler shifts
% current time is the time in minutes since epoch (elapsed time) that all sats are in view.
% julian date is JD when sats are in view needed  to convert ECI to ECEF


%constants
c= 3e8;
longitude=77.3790;
latitude=38.9649; 
H=0; 
%1085027.1960869177710264921188354
%-4845789.2537157507613301277160645
%3989288.0202742042019963264465332
[x,y,z]=lla2ecef_AB(latitude*(2*pi/360),(2*pi)-(longitude*(pi/180)),H); % this function takes east longitude. change if need be.
rec_pos(1,1) =  x;
rec_pos(2,1) =  y+100;
rec_pos(3,1) =  z;
%rec_pos(1,1) = 0;
%rec_pos(2,1) = 0;
%rec_pos(3,1) = 0;
c_rec_clock_bias =0;
rec_vel = [1;1;1];
c_rec_clock_bias_rate=0;

y_i = [rec_pos; c_rec_clock_bias; rec_vel; c_rec_clock_bias_rate];
y_0 =y_i;


options = optimoptions('lsqnonlin', 'FunctionTolerance', 1e-15, 'FiniteDifferenceStepSize',1e-15);
options.Algorithm = 'trust-region-reflective';

g = @(y_i) find_error_1(y_i,shift, pos, lambda, ephemeris, JD_prop_to);
state = lsqnonlin(g,y_0,[],[], options);
vpa(state)

rec_GDOP = 4;
Error  = 4;
delta_y_0  = 4;
error_0 = 4;
A_0  =4;




end