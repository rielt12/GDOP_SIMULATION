function alpha = find_step_size(y_i,A,shift, pos, lambda, current_time, ephemeris, JD_prop_to, Error,iter, error,R, delta_y)
alpha0  = 1;

g = @(alpha) find_error_1(y_i,shift, pos, lambda, current_time, ephemeris, JD_prop_to, error, alpha, delta_y);

[alpha,fval] = fminunc(g,alpha0);


% y_s = y_i;
% y_i =y_s + alpha*delta_y;
% Error(iter+1) = find_error(y_i,shift, pos, lambda, current_time, ephemeris, JD_prop_to, error);
% 
% 
% if(Error(iter+1) > Error(iter))
%     %disp('entered')
% end
% reps = 1;
% while(1)
% %disp('here')
% 
% % Error(iter)
% % vpa(alpha) 
% % Error(iter+1)
% alpha = alpha/2;
% y_i =y_s + alpha*delta_y;
% Error(iter+1) = find_error(y_i,shift, pos, lambda, current_time, ephemeris, JD_prop_to, error);
% if(Error(iter+1)< Error(iter))
%     break;
% end
% reps = reps+1;
% if (reps >50)
%     break;
% 
% end
% 
% if(vpa(Error(iter+1)) == vpa(Error(iter)))
%     %disp('equal')
% end
% 
% if(vpa(Error(iter+1)) < vpa(Error(iter)))
%     %disp('exited')
% end
end
