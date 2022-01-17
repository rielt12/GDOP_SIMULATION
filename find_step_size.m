function alpha = find_step_size(y_i,A,shift, pos, lambda, current_time, ephemeris, JD_prop_to, Error,iter, error)
alpha  = 2;
y_s = y_i;
delta_y =inv(A'*A)*A'*error';

size(Error)



while(Error(iter-1,1) > Error(iter-2,1))
display('entered')
 alpha = alpha/2;
 y_i=y_s+alpha*delta_y;
 rec_pos =y_i(1:3,1);
 rec_clock_bias =y_i(4,1);
 rec_vel =y_i(5:7,1);
 rec_clock_bias_rate= y_i(8,1);

 
 t_R =current_time; % minutes


 R  =zeros(length(shift),length(shift));
 %first compute error terms
 for i=1:length(shift)
    del_ADR(i) = acculumulated_delta_range_derivative(rec_pos, rec_clock_bias,t_R, rec_vel, rec_clock_bias_rate, pos(i),lambda,ephemeris(i), JD_prop_to);
    error(i)=lambda*(-del_ADR(i)/lambda-shift(i));
    R(i,i) = 0.01;
 end
Error(iter-1,1) = abs(norm(error));


if (Error(iter-1) < Error(iter-2)) 
display('succ')
end
if (Error(iter-1) > Error(iter-2))
display('fail')
end



end



end