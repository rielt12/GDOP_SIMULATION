function [state rec_GDOP, Error, delta_y_0,error_0,A_0, del_ADR] = doppler_shift_positioning_m_gps(shift, pos, lambda, current_time, ephemeris, svids)

max_iter = 15e3;

% inputs
% shifts are observed doppler shifts multiplied by wavelength;


%constants
c= 3e8;
longitude=34.7984651;
latitude=32.1163502; 
H=33; 

% longitude=34.7984651;
% latitude=32.1163502; 
% H=33.000000; 
%4440231.992452892474830150604248, 
%3085865.2224882706068456172943115,
%3371383.3935853922739624977111816
[x,y,z]=lla2ecef_AB(latitude*(2*pi/360),longitude*(pi/180),H); % this function takes east longitude. change if need be.
rec_pos(1,1) =  x;
rec_pos(2,1) =  y;
rec_pos(3,1) =  z;
%rec_pos(1,1) = 0;
%rec_pos(2,1) = 0;
%rec_pos(3,1) = 0;
c_rec_clock_bias =0;
rec_vel = [0;0;0];
c_rec_clock_bias_rate= 0;



y_i = [rec_pos; c_rec_clock_bias; rec_vel; c_rec_clock_bias_rate];
y_0 = y_i;


options = optimoptions('lsqnonlin', 'FunctionTolerance', 1e-30, 'FiniteDifferenceStepSize',1e-5, 'Display','iter');
options.Algorithm = 'trust-region-reflective';
options.MaxFunctionEvaluations = 8000;
options.StepTolerance= 1e-3;
options.MaxIterations = 1e5;
options.OptimalityTolerance = 1e-20;

options

%here
g = @(y_i) compute_cost_doppler_gps(y_i , shift, pos, lambda, current_time, ephemeris, svids);
[state,resnorm,residual,exitflag,output,~,jacobian] = lsqnonlin(g,y_0,[],[], options);
vpa(resnorm)
vpa(0.5*norm(residual)^2)
vpa(state)

A_0  =jacobian;
vpa(A_0)


rec_GDOP = 4;
Error  = 4;
delta_y_0  = 4;
error_0 = 4;




rec_pos = y_i(1:3,1);
rec_clock_bias = y_i(4,1);
rec_vel = y_i(5:7,1);
rec_clock_bias_rate = y_i(8);
bias = [ 
106.56481186797157079126918688416
  106.106679178774356842041015625
108.28082670271470533407409675419
108.43689345071717866630933713168
106.06789046277521038064151071012
107.72690038631384368272847495973
106.26553490012906877382192760706
 105.1043247431520626378187444061];
for i=1:length(shift)
    del_ADR(i) = acculumulated_delta_range_derivative_gps(rec_pos, rec_clock_bias,current_time, rec_vel, rec_clock_bias_rate, pos(i,:)',lambda,ephemeris, svids(i))-bias(i);
    R(i,i) = (0.01)^2;
    error(i,1)= lambda*(shift(i)-del_ADR(i));
end
error



%here
%vpa(inv(A_0)*error)



end