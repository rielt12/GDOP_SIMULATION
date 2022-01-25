function [state rec_GDOP, Error] = doppler_shift_positioning(shift, pos, lambda, current_time, ephemeris, JD_prop_to)

% inputs
% shifts are observed doppler shifts
% current time is the time in minutes since epoch (elapsed time) that all sats are in view.
% julian date is JD when sats are in view needed  to convert ECI to ECEF


%constants
c= 3e8;
longitude=77.3790;
latitude=38.9649; 
H=0; 
[x,y,z]=lla2ecef_AB(latitude*(2*pi/360),(2*pi)-(longitude*(pi/180)),H); % this function takes east longitude. change if need be.
rec_pos(1,1) =  x+100;
rec_pos(2,1) =  y+100;
rec_pos(3,1) =  z+100;
%rec_pos(1,1) = 0;
%rec_pos(2,1) = 0;
%rec_pos(3,1) = 0;
rec_clock_bias =0;
rec_vel = [100;0;0];
rec_clock_bias_rate= 0;

y_i = [rec_pos; rec_clock_bias; rec_vel; rec_clock_bias_rate];
y_0 =y_i;


t_R =current_time; % minutes

R  =zeros(length(shift),length(shift));
%first compute error terms
for i=1:length(shift)
    elapsedtime=(JD_prop_to-ephemeris(i).jdsatepoch)*24*60;
    del_ADR(i) = acculumulated_delta_range_derivative(rec_pos, rec_clock_bias,elapsedtime, rec_vel, rec_clock_bias_rate, pos(i,:)',lambda,ephemeris(i,1), JD_prop_to);
    R(i,i) = (0.01)^2;
    error(i)= lambda*((shift(i)-del_ADR(i)));
end



% form the jacobian 
for i=1:length(shift)
A(i,:) = Jacobian_Psiaki_Row(rec_pos, rec_clock_bias,elapsedtime, rec_vel, rec_clock_bias_rate, pos(i),lambda,ephemeris(i), JD_prop_to);
end







iter = 0;
Error(iter+1,1) = abs(norm(error));


while(norm(error) >1e-3 && iter < 1000)
 

iter = iter+1
norm(error)
% tau =100;
% delta_y =pinv((A'*A))*A'*R*error';
% %delta_y = lschol(A'*R*A,A'*R*error');
% %delta_y = lsqr((A'*A),A'*error');
% %delta_y = (A'*R*A)\A'*R*error';
% 
% y_i = y_i +tau*delta_y;


rec_pos =y_i(1:3,1);
rec_clock_bias =y_i(4,1);
rec_vel =y_i(5:7,1);
rec_clock_bias_rate= y_i(8,1);


t_R =current_time; % minutes


R  =zeros(length(shift),length(shift));
%first compute error terms
 for i=1:length(shift)
    elapsedtime=(JD_prop_to-ephemeris(i).jdsatepoch)*24*60;
    del_ADR(i) = acculumulated_delta_range_derivative(rec_pos, rec_clock_bias,elapsedtime, rec_vel, rec_clock_bias_rate, pos(i,:)',lambda,ephemeris(i,1), JD_prop_to);
    R(i,i) = (0.01)^2;
    error(i)= lambda*((shift(i)-del_ADR(i)));
 end
 Error(iter+1,1) = norm(error);


for i=1:length(shift)
 A(i,:) = Jacobian_Psiaki_Row(rec_pos, rec_clock_bias,elapsedtime, rec_vel, rec_clock_bias_rate, pos(i),lambda,ephemeris(i), JD_prop_to);
end



 
end




rec_GDOP = sqrt(trace(inv(A'*A)))
state = y_i;
vpa(state)
vpa(state - y_0)

end