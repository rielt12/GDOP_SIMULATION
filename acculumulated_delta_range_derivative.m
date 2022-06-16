function del_ADR = acculumulated_delta_range_derivative(r, del_R,t_R, v, del_R_rate, pos,lambda, ephemeris, JD_prop_to)
del_t_R = 0.01; 
ADR_1 = acculumulated_delta_range(r-(2*v*del_t_R)/(1+del_R_rate),del_R-(2*del_R_rate*del_t_R)/(1+del_R_rate),(t_R-2*del_t_R),pos, lambda,ephemeris, JD_prop_to);
ADR_2 = acculumulated_delta_range(r-(2*v*del_t_R)/(1+del_R_rate),del_R-(del_R_rate*del_t_R)/(1+del_R_rate),(t_R-del_t_R),pos, lambda,ephemeris, JD_prop_to);
ADR_3 = acculumulated_delta_range(r+(v*del_t_R)/(1+del_R_rate),del_R+(del_R_rate*del_t_R)/(1+del_R_rate),(t_R+del_t_R),pos, lambda, ephemeris, JD_prop_to);
ADR_4 = acculumulated_delta_range(r+(2*v*del_t_R)/(1+del_R_rate),del_R+(2*del_R_rate*del_t_R)/(1+del_R_rate),(t_R+2*del_t_R),pos, lambda,ephemeris,JD_prop_to);


del_ADR = (ADR_1-8*ADR_2+8*ADR_3-ADR_4)/(12*del_t_R);
end