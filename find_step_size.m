function alpha = find_step_size(y_i,A,shift, pos, lambda, current_time, ephemeris, JD_prop_to, Error,iter, error,R, delta_y)
alpha  = 1;
y_i =y_i + alpha*delta_y;
Error(iter+1) = find_error(y_i,shift, pos, lambda, current_time, ephemeris, JD_prop_to, error);


% alpha = 1;
% y_s = y_i;
% %tau =find_step_size(y_i,A,shift, pos, lambda, current_time, ephemeris, JD_prop_to, Error,iter, error,R, delta_y);
% Error(iter+1) = find_error(y_i,shift, pos, lambda, current_time, ephemeris, JD_prop_to, error);
% if (Error(iter+1) > Error(iter))
% alpha = alpha/2;
% y_s =y_s + alpha*delta_y;
% Error(iter+1) = find_error(y_s,shift, pos, lambda, current_time, ephemeris, JD_prop_to, error);
% end

while(Error(iter+1) > Error(iter))
% Error(iter)
% vpa(alpha) 
% Error(iter+1)
alpha = alpha/2;
y_i =y_i + alpha*delta_y;
Error(iter+1) = find_error(y_i,shift, pos, lambda, current_time, ephemeris, JD_prop_to, error);


end


