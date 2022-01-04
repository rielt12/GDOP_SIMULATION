function [state rec_GDOP] = doppler_shift_positioning(shift, pos,vel, lambda, t_rec)
%constants
c= 3e8;

rec_pos = [0;0;0];
rec_clock_bias = 0;
rec_vel = [0;0;0];
rec_clock_bias_rate =0;
state = [rec_pos; rec_clock_bias; rec_vel; rec_clock_bias_rate];
t_R =t_rec;

R  =zeros(length(shift),length(shift));
%first compute derivatives
for i=1:length(shift)
    del_ADR(i) = acculumulated_delta_range_derivative(rec_pos, rec_clock_bias,t_R, rec_vel, rec_clock_bias_rate, pos(i),lambda);
    error(i)= ADR(i)+lambda*shift(i);
    R(i,i) = 0.01;
end
A 

% compute error terms

% compute update terms


%while loop




end