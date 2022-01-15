function [state rec_GDOP] = doppler_shift_positioning(shift, pos, lambda, current_time, ephemeris, JD_prop_to)

% inputs
% shifts are observed doppler shifts
% current time is the time in minutes since epoch (elapsed time) that all sats are in view.
% julian date is JD when sats are in view needed  to convert ECI to ECEF


%constants
c= 3e8;
y_i= zeros(8,1);
rec_pos =y_i(1:3,1);
rec_clock_bias =y_i(4,1);
rec_vel =y_i(5:7,1);
rec_clock_bias_rate= y_i(8,1);



t_R =current_time; % minutes

R  =zeros(length(shift),length(shift));
%first compute error terms
for i=1:length(shift)
    del_ADR(i) = acculumulated_delta_range_derivative(rec_pos, rec_clock_bias,t_R, rec_vel, rec_clock_bias_rate, pos(i),lambda,ephemeris(i), JD_prop_to);
    error(i)= -1*(del_ADR(i)+lambda*shift(i));
    R(i,i) = 0.01;
end

% form the jacobian 
for i=1:length(shift)
A(i,:) = Jacobian_Psiaki_Row(rec_pos, rec_clock_bias,t_R, rec_vel, rec_clock_bias_rate, pos(i),lambda,ephemeris(i), JD_prop_to);
end
vpa(A)

inv(A'*A)*A'



j=0;
while (abs(norm(error))>0.000000001)
 abs(norm(error))
  j=j+1
 delta_y =inv(A'*A)*A'*error';

 y_i=y_i+delta_y;
 rec_pos =y_i(1:3,1);
 rec_clock_bias =y_i(4,1);
 rec_vel =y_i(5:7,1);
 rec_clock_bias_rate= y_i(8,1);


t_R =current_time; % minutes

R  =zeros(length(shift),length(shift));
%first compute error terms
for i=1:length(shift)
    del_ADR(i) = acculumulated_delta_range_derivative(rec_pos, rec_clock_bias,t_R, rec_vel, rec_clock_bias_rate, pos(i),lambda,ephemeris(i), JD_prop_to);
    error(i)= -(del_ADR(i)+lambda*shift(i));
    R(i,i) = 0.01;
end

% form the jacobian 
for i=1:length(shift)
A(i,:) = Jacobian_Psiaki_Row(rec_pos, rec_clock_bias,t_R, rec_vel, rec_clock_bias_rate, pos(i),lambda,ephemeris(i), JD_prop_to);
end



end





rec_GDOP = 0;


% compute error terms

% compute update terms


%while loop

state = y_i;


end