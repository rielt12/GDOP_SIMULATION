function [state rec_GDOP, Error, delta_y_0,error_0,A_0, del_ADR] = doppler_shift_positioning_m(shift, pos, lambda, current_time, ephemeris, JD_prop_to)


max_iter = 1500;

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
c_rec_clock_bias =c*101e-3;
rec_vel = [0;0;0];
c_rec_clock_bias_rate=0;

y_i = [rec_pos; c_rec_clock_bias; rec_vel; c_rec_clock_bias_rate];
y_0 =y_i;


options = optimoptions('lsqnonlin', 'FunctionTolerance', 1e-30, 'FiniteDifferenceStepSize',1e-5, 'Display','iter');
options.Algorithm = 'trust-region-reflective';
options.MaxFunctionEvaluations = 8000;
options.StepTolerance= 1e-20;
options.MaxIterations = 1500;
options.OptimalityTolerance = 1e-15;

options

g = @(y_i) compute_cost_doppler(y_i,shift, pos, lambda, ephemeris, JD_prop_to);
[state,resnorm,residual,exitflag,output,~,jacobian] = lsqnonlin(g,y_0,[],[], options);
vpa(resnorm)
vpa(0.5*norm(residual)^2)
%vpa(state)

A_0  =jacobian;
vpa(A_0)


rec_GDOP = 4;
Error  = 4;
delta_y_0  = 4;
error_0 = 4;




rec_pos = y_i(1:3,1);
rec_clock_bias = y_i(4,1)/c;
rec_vel = y_i(5:7,1);
rec_clock_bias_rate = y_i(8)/c;
for i=1:length(shift)
    elapsedtime=(JD_prop_to-ephemeris(i).jdsatepoch)*24*60;
    del_ADR(i) = acculumulated_delta_range_derivative(rec_pos, rec_clock_bias,elapsedtime, rec_vel, rec_clock_bias_rate, pos(i,:)',lambda,ephemeris(i,1), JD_prop_to);
    R(i,i) = (0.01)^2;
    error(i,1)= lambda*(shift(i)-del_ADR(i));
end
error


inv(A_0)*error

end