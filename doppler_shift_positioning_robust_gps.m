function [state rec_GDOP, Error,y_0,delta_y_0,error_0,A_0, del_ADR] = doppler_shift_positioning_robust_gps(shift, pos, lambda, current_time, ephemeris, svids)

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
rec_pos(1,1) =  x+1000;
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

R  =zeros(length(shift),length(shift));

for i=1:length(shift)
   R(i,i) = (0.01)^2;
end



error = compute_cost_doppler_gps(y_i , shift, pos, lambda, current_time, ephemeris, svids);
error_0 = error;


for i=1:length(shift)
A(i,:) = Jacobian_Psiaki_Row_Numerical_Central_O2_gps(rec_pos, c_rec_clock_bias,rec_vel, c_rec_clock_bias_rate, shift(i), pos(i), lambda, current_time, ephemeris, svids(i));
end
A_0 = A;

vpa(A_0)



delta_y_0 =inv(A_0)*(error');





iter = 1;
nwd = 0;
Error(iter,1) = 0.5*norm(error)^2;
while(0.5*norm(error)^2 >1e-8 && iter < max_iter && nwd~=1)

iter  
0.5*norm(error)^2
Error(iter,1) = 0.5*norm(error)^2;

delta_y =inv((A'*R*A))*A'*R*error';

%delta_y


N_g  =1000;
[tau,nwd] = line_search_doppler_gps(y_i, delta_y, N_g ,shift, pos, lambda, current_time, ephemeris,svids);
%tau = 1e-3;


delta_y(4,1) = delta_y(4,1)/c; % helps with convergence
delta_y(8,1) = delta_y(8,1);
%y_i = y_i -tau*delta_y;
y_i(1) = y_i(1)-tau*delta_y(1,1);
y_i(2) = y_i(2)-tau*delta_y(2,1);
y_i(3) = y_i(3)-tau*delta_y(3,1);
  y_i(4) = y_i(4)-tau*delta_y(4,1);
 y_i(5) = y_i(5)-tau*delta_y(5,1);
 y_i(6) = y_i(6)-tau*delta_y(6,1);
 y_i(7) = y_i(7)-tau*delta_y(7,1);
  y_i(8) = y_i(8)-tau*delta_y(8,1);
%y_i(4:8,1)  = [0;0;0;0;0]; % static condition % and no clockbias
%y_i(5:7) = [0;0;0]; % static condition


state_record(iter,:) = y_i';




rec_pos =y_i(1:3,1);
c_rec_clock_bias =y_i(4,1)/c;
rec_vel =y_i(5:7,1);
c_rec_clock_bias_rate= y_i(8,1)/c;


t_R =current_time; % gps tow


for i=1:length(shift)
   R(i,i) = (0.01)^2;
end



error = compute_cost_doppler_gps(y_i , shift, pos, lambda, current_time, ephemeris, svids);
c=3e8;


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
    del_ADR(i) = acculumulated_delta_range_derivative_gps(rec_pos, c_rec_clock_bias/c,current_time, rec_vel, c_rec_clock_bias_rate/c, pos(i,:)',lambda,ephemeris, svids(i))-bias(i);
end


for i=1:length(shift)
A(i,:) = Jacobian_Psiaki_Row_Numerical_Central_O2_gps(rec_pos, c_rec_clock_bias/c,rec_vel, c_rec_clock_bias_rate/c, shift(i), pos(i), lambda, current_time, ephemeris, svids(i));
end


% plotting shifts
figure(400)
clf
num=1:8;
scatter(num,shift,'MarkerFaceColor','blue')
hold on
scatter(num,del_ADR,'MarkerFaceColor','red')
xlabel('satellite number')
ylabel('ADR Derivatives')
legend('measured','estimated')

iter = iter+1;



end

rec_pos = y_i(1:3,1);
c_rec_clock_bias = y_i(4,1)/c;
rec_vel = y_i(5:7,1);
c_rec_clock_bias_rate = y_i(8,1)/c;



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
    del_ADR(i) = acculumulated_delta_range_derivative_gps(rec_pos, c_rec_clock_bias/c,current_time, rec_vel, c_rec_clock_bias_rate/c, pos(i,:)',lambda,ephemeris, svids(i))-bias(i);
end



rec_GDOP = sqrt(trace(inv(A'*A)));
state = y_i;


end