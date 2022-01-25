function [Del_ADR_Del_r_dot] = del_ADR_del_r_dot(r, del_R,t_R, v, del_R_rate, pos,lambda, ephemeris, JD_prop_to)
del_t_R = 0.1; 

% derivative wrt to del_r_dot
del_R_rate_b = del_R_rate+1e-6;
ADR_1_a = acculumulated_delta_range(r-(2*v*del_t_R)/(1+del_R_rate),del_R-(2*del_R_rate*del_t_R)/(1+del_R_rate),(t_R-2*del_t_R),pos, lambda,ephemeris, JD_prop_to);
ADR_1_b = acculumulated_delta_range(r-(2*v*del_t_R)/(1+del_R_rate_b),del_R-(2*del_R_rate_b*del_t_R)/(1+del_R_rate_b),(t_R-2*del_t_R),pos, lambda,ephemeris, JD_prop_to);
D_ADR_1 = (ADR_1_b-ADR_1_a)/1e-6;
ADR_2_a = acculumulated_delta_range(r-(2*v*del_t_R)/(1+del_R_rate),del_R-(del_R_rate*del_t_R)/(1+del_R_rate),(t_R-del_t_R),pos, lambda,ephemeris, JD_prop_to);
ADR_2_b = acculumulated_delta_range(r-(2*v*del_t_R)/(1+del_R_rate_b),del_R-(del_R_rate_b*del_t_R)/(1+del_R_rate_b),(t_R-del_t_R),pos, lambda,ephemeris, JD_prop_to);
D_ADR_2 = (ADR_2_b-ADR_2_a)/1e-6;
ADR_3_a = acculumulated_delta_range(r+(v*del_t_R)/(1+del_R_rate),del_R+(del_R_rate*del_t_R)/(1+del_R_rate),(t_R+del_t_R),pos, lambda, ephemeris, JD_prop_to);
ADR_3_b = acculumulated_delta_range(r+(v*del_t_R)/(1+del_R_rate_b),del_R+(del_R_rate_b*del_t_R)/(1+del_R_rate_b),(t_R+del_t_R),pos, lambda, ephemeris, JD_prop_to);
D_ADR_3 = (ADR_3_b-ADR_3_a)/1e-6;
ADR_4_a = acculumulated_delta_range(r+(2*v*del_t_R)/(1+del_R_rate),del_R+(2*del_R_rate*del_t_R)/(1+del_R_rate),(t_R+2*del_t_R),pos, lambda,ephemeris,JD_prop_to);
ADR_4_b = acculumulated_delta_range(r+(2*v*del_t_R)/(1+del_R_rate_b),del_R+(2*del_R_rate_b*del_t_R)/(1+del_R_rate_b),(t_R+2*del_t_R),pos, lambda,ephemeris,JD_prop_to);
D_ADR_4 = (ADR_4_b-ADR_4_a)/1e-6;
Del_ADR_Del_r_dot = (D_ADR_1-8*D_ADR_2+8*D_ADR_3-D_ADR_4)/(12*del_t_R);


end