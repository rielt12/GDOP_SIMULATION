function [state rec_GDOP, Error, delta_y_0,error_0,A_0, del_ADR] = doppler_shift_positioning_robust(shift, pos, lambda, current_time, ephemeris, JD_prop_to)

max_iter = 15e3;

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
rec_pos(1,1) =  x+100e4;
rec_pos(2,1) =  y+100e4;
rec_pos(3,1) =  z+100e4;
%rec_pos(1,1) = 0;
%rec_pos(2,1) = 0;
%rec_pos(3,1) = 0;
c_rec_clock_bias =1;
rec_vel = [1;1;1];
c_rec_clock_bias_rate= 1;

y_i = [rec_pos; c_rec_clock_bias; rec_vel; c_rec_clock_bias_rate];

R  =zeros(length(shift),length(shift));

for i=1:length(shift)
   R(i,i) = (0.01)^2;
end



error = compute_cost_doppler(y_i , shift, pos, lambda,ephemeris,JD_prop_to);
error_0 = error;
0.5*norm(error_0)^2

for i=1:length(shift)
elapsedtime=(JD_prop_to-ephemeris(i).jdsatepoch)*24*60;
A(i,:) = Jacobian_Psiaki_Row_Numerical_Central_O2(rec_pos, c_rec_clock_bias,elapsedtime, rec_vel, c_rec_clock_bias_rate, pos(i), shift(i),lambda,ephemeris(i), JD_prop_to);
end
A_0 = A;






delta_y_0 =inv(A_0)*(error')





iter = 1;
nwd = 0;
Error(iter,1) = 0.5*norm(error)^2;
while(0.5*norm(error)^2 >1e-8 && iter < max_iter && nwd~=1)

iter  
0.5*norm(error)^2
Error(iter,1) = 0.5*norm(error)^2;

delta_y =inv((A'*R*A))*A'*R*error';


 N_g  =1000;



%[tau, nwd]  =line_search_doppler(y_i, delta_y, N_g ,shift, pos, lambda,ephemeris,JD_prop_to);
tau = 0.5;


delta_y(4,1) = delta_y(4,1)/c;
delta_y(8,1) = delta_y(8,1);
y_i = y_i -tau*delta_y;

%delta_y
state_record(iter,:) = y_i';
tau



rec_pos =y_i(1:3,1);
c_rec_clock_bias =y_i(4,1);
rec_vel =y_i(5:7,1);
c_rec_clock_bias_rate= y_i(8,1);


t_R =current_time; % minutes


for i=1:length(shift)
   R(i,i) = (0.01)^2;
end



error = compute_cost_doppler(y_i , shift, pos, lambda,ephemeris,JD_prop_to);
c=3e8;

rec_pos = y_i(1:3,1);
rec_clock_bias = y_i(4,1)/c;
rec_vel = y_i(5:7,1);
rec_clock_bias_rate = y_i(8,1)/c;


for i=1:length(shift)
    elapsedtime=(JD_prop_to-ephemeris(i).jdsatepoch)*24*60;
    del_ADR(i) = acculumulated_delta_range_derivative(rec_pos, rec_clock_bias,elapsedtime, rec_vel, rec_clock_bias_rate, pos(i,:)',lambda,ephemeris(i,1), JD_prop_to);
end


for i=1:length(shift)
A(i,:) = Jacobian_Psiaki_Row_Numerical_Central_O2(rec_pos, c_rec_clock_bias,elapsedtime, rec_vel, c_rec_clock_bias_rate, pos(i), shift(i),lambda,ephemeris(i), JD_prop_to);
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
rec_clock_bias = y_i(4,1)/c;
rec_vel = y_i(5:7,1);
rec_clock_bias_rate = y_i(8,1)/c;

for i=1:length(shift)
    elapsedtime=(JD_prop_to-ephemeris(i).jdsatepoch)*24*60;
    del_ADR(i) = acculumulated_delta_range_derivative(rec_pos, rec_clock_bias,elapsedtime, rec_vel, rec_clock_bias_rate, pos(i,:)',lambda,ephemeris(i,1), JD_prop_to);
end



rec_GDOP = sqrt(trace(inv(A'*A)));
state = y_i;


end