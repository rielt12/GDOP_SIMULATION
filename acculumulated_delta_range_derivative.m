function del_ADR = acculumulated_delta_range_derivative(r, del_R,t_R, v, del_R_rate, pos,lambda)
del_t_R = 0.1; 
ADR_1 = acculumulated_delta_range(r-(2*v*del_t_R)/(1+del_R_rate),del_R-(2*del_R_rate*del_t_R)/(1+del_R_rate),(t_R-2*del_t_R),pos, lambda);
ADR_2 = acculumulated_delta_range(r-(2*v*del_t_R)/(1+del_R_rate),del_R-(del_R_rate*del_t_R)/(1+del_R_rate),(t_R-del_t_R),pos, lambda);
ADR_3 = acculumulated_delta_range(r+(v*del_t_R)/(1+del_R_rate),del_R+(del_R_rate*del_t_R)/(1+del_R_rate),(t_R+del_t_R),pos, lambda);
ADR_4 = acculumulated_delta_range(r+(2*v*del_t_R)/(1+del_R_rate),del_R+(2*del_R_rate*del_t_R)/(1+del_R_rate),(t_R+2*del_t_R),pos, lambda);

del_ADR = (ADR_1-8*ADR_2+8*ADR_3-ADR_4)/(12*del_t_R);
end