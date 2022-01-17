function [del_ADR_del_v_x, del_ADR_del_v_y, del_ADR_del_v_z] = del_ADR_del_v(r, del_R,t_R, v, del_R_rate, pos,lambda, ephemeris, JD_prop_to)
del_t_R = 0.1; 

% derivative wrt to v_x
v_b(1,1) = v(1,1) + 0.01;
v_b(2,1) = v(2,1);
v_b(3,1) = v(3,1);
ADR_1_a = acculumulated_delta_range(r-(2*v*del_t_R)/(1+del_R_rate),del_R-(2*del_R_rate*del_t_R)/(1+del_R_rate),(t_R-2*del_t_R),pos, lambda,ephemeris, JD_prop_to);
ADR_1_b = acculumulated_delta_range(r-(2*v_b*del_t_R)/(1+del_R_rate),del_R-(2*del_R_rate*del_t_R)/(1+del_R_rate),(t_R-2*del_t_R),pos, lambda,ephemeris, JD_prop_to);
D_ADR_1 = (ADR_1_b-ADR_1_a)/0.01;
ADR_2_a = acculumulated_delta_range(r-(2*v*del_t_R)/(1+del_R_rate),del_R-(del_R_rate*del_t_R)/(1+del_R_rate),(t_R-del_t_R),pos, lambda,ephemeris, JD_prop_to);
ADR_2_b = acculumulated_delta_range(r-(2*v_b*del_t_R)/(1+del_R_rate),del_R-(del_R_rate*del_t_R)/(1+del_R_rate),(t_R-del_t_R),pos, lambda,ephemeris, JD_prop_to);
D_ADR_2 = (ADR_2_b-ADR_2_a)/0.01;
ADR_3_a = acculumulated_delta_range(r+(v*del_t_R)/(1+del_R_rate),del_R+(del_R_rate*del_t_R)/(1+del_R_rate),(t_R+del_t_R),pos, lambda, ephemeris, JD_prop_to);
ADR_3_b = acculumulated_delta_range(r+(v_b*del_t_R)/(1+del_R_rate),del_R+(del_R_rate*del_t_R)/(1+del_R_rate),(t_R+del_t_R),pos, lambda, ephemeris, JD_prop_to);
D_ADR_3 = (ADR_3_b-ADR_3_a)/0.01;
ADR_4_a = acculumulated_delta_range(r+(2*v*del_t_R)/(1+del_R_rate),del_R+(2*del_R_rate*del_t_R)/(1+del_R_rate),(t_R+2*del_t_R),pos, lambda,ephemeris,JD_prop_to);
ADR_4_b = acculumulated_delta_range(r+(2*v_b*del_t_R)/(1+del_R_rate),del_R+(2*del_R_rate*del_t_R)/(1+del_R_rate),(t_R+2*del_t_R),pos, lambda,ephemeris,JD_prop_to);
D_ADR_4 = (ADR_4_b-ADR_4_a)/0.01;
del_ADR_del_v_x = (D_ADR_1-8*D_ADR_2+8*D_ADR_3-D_ADR_4)/(12*del_t_R);

% derivative wrt to v_y
v_b(1,1) = v(1,1);
v_b(2,1) = v(2,1)+0.01;
v_b(3,1) = v(3,1);
ADR_1_a = acculumulated_delta_range(r-(2*v*del_t_R)/(1+del_R_rate),del_R-(2*del_R_rate*del_t_R)/(1+del_R_rate),(t_R-2*del_t_R),pos, lambda,ephemeris, JD_prop_to);
ADR_1_b = acculumulated_delta_range(r-(2*v_b*del_t_R)/(1+del_R_rate),del_R-(2*del_R_rate*del_t_R)/(1+del_R_rate),(t_R-2*del_t_R),pos, lambda,ephemeris, JD_prop_to);
D_ADR_1 = (ADR_1_b-ADR_1_a)/0.01;
ADR_2_a = acculumulated_delta_range(r-(2*v*del_t_R)/(1+del_R_rate),del_R-(del_R_rate*del_t_R)/(1+del_R_rate),(t_R-del_t_R),pos, lambda,ephemeris, JD_prop_to);
ADR_2_b = acculumulated_delta_range(r-(2*v_b*del_t_R)/(1+del_R_rate),del_R-(del_R_rate*del_t_R)/(1+del_R_rate),(t_R-del_t_R),pos, lambda,ephemeris, JD_prop_to);
D_ADR_2 = (ADR_2_b-ADR_2_a)/0.01;
ADR_3_a = acculumulated_delta_range(r+(v*del_t_R)/(1+del_R_rate),del_R+(del_R_rate*del_t_R)/(1+del_R_rate),(t_R+del_t_R),pos, lambda, ephemeris, JD_prop_to);
ADR_3_b = acculumulated_delta_range(r+(v_b*del_t_R)/(1+del_R_rate),del_R+(del_R_rate*del_t_R)/(1+del_R_rate),(t_R+del_t_R),pos, lambda, ephemeris, JD_prop_to);
D_ADR_3 = (ADR_3_b-ADR_3_a)/0.01;
ADR_4_a = acculumulated_delta_range(r+(2*v*del_t_R)/(1+del_R_rate),del_R+(2*del_R_rate*del_t_R)/(1+del_R_rate),(t_R+2*del_t_R),pos, lambda,ephemeris,JD_prop_to);
ADR_4_b = acculumulated_delta_range(r+(2*v_b*del_t_R)/(1+del_R_rate),del_R+(2*del_R_rate*del_t_R)/(1+del_R_rate),(t_R+2*del_t_R),pos, lambda,ephemeris,JD_prop_to);
D_ADR_4 = (ADR_4_b-ADR_4_a)/0.01;
del_ADR_del_v_y = (D_ADR_1-8*D_ADR_2+8*D_ADR_3-D_ADR_4)/(12*del_t_R);

% derivative wrt to v_z
v_b(1,1) = v(1,1);
v_b(2,1) = v(2,1);
v_b(3,1) = v(3,1)+0.01;
ADR_1_a = acculumulated_delta_range(r-(2*v*del_t_R)/(1+del_R_rate),del_R-(2*del_R_rate*del_t_R)/(1+del_R_rate),(t_R-2*del_t_R),pos, lambda,ephemeris, JD_prop_to);
ADR_1_b = acculumulated_delta_range(r-(2*v_b*del_t_R)/(1+del_R_rate),del_R-(2*del_R_rate*del_t_R)/(1+del_R_rate),(t_R-2*del_t_R),pos, lambda,ephemeris, JD_prop_to);
D_ADR_1 = (ADR_1_b-ADR_1_a)/0.01;
ADR_2_a = acculumulated_delta_range(r-(2*v*del_t_R)/(1+del_R_rate),del_R-(del_R_rate*del_t_R)/(1+del_R_rate),(t_R-del_t_R),pos, lambda,ephemeris, JD_prop_to);
ADR_2_b = acculumulated_delta_range(r-(2*v_b*del_t_R)/(1+del_R_rate),del_R-(del_R_rate*del_t_R)/(1+del_R_rate),(t_R-del_t_R),pos, lambda,ephemeris, JD_prop_to);
D_ADR_2 = (ADR_2_b-ADR_2_a)/0.01;
ADR_3_a = acculumulated_delta_range(r+(v*del_t_R)/(1+del_R_rate),del_R+(del_R_rate*del_t_R)/(1+del_R_rate),(t_R+del_t_R),pos, lambda, ephemeris, JD_prop_to);
ADR_3_b = acculumulated_delta_range(r+(v_b*del_t_R)/(1+del_R_rate),del_R+(del_R_rate*del_t_R)/(1+del_R_rate),(t_R+del_t_R),pos, lambda, ephemeris, JD_prop_to);
D_ADR_3 = (ADR_3_b-ADR_3_a)/0.01;
ADR_4_a = acculumulated_delta_range(r+(2*v*del_t_R)/(1+del_R_rate),del_R+(2*del_R_rate*del_t_R)/(1+del_R_rate),(t_R+2*del_t_R),pos, lambda,ephemeris,JD_prop_to);
ADR_4_b = acculumulated_delta_range(r+(2*v_b*del_t_R)/(1+del_R_rate),del_R+(2*del_R_rate*del_t_R)/(1+del_R_rate),(t_R+2*del_t_R),pos, lambda,ephemeris,JD_prop_to);
D_ADR_4 = (ADR_4_b-ADR_4_a)/0.01;
del_ADR_del_v_z = (D_ADR_1-8*D_ADR_2+8*D_ADR_3-D_ADR_4)/(12*del_t_R);



end