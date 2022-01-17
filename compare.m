load('last_sat_sim.mat')
pos1 = r_ecef;
v1 = v_rec;
lambda1 = lambda;
elapsedtime1 = elapsedtime;
rec_pos1 = ecef_ref;
del_r_dot1 = del_r_dot;
del_r1 = del_r;
ephemeris1 = eph;
JD_prop_to_1= JD_prop_to;
del_ADR1 = acculumulated_delta_range_derivative(rec_pos1, del_r1,elapsedtime1, v1, del_r_dot1, pos1,lambda1, ephemeris1, JD_prop_to_1);

load('last_sat_NLS.mat')
pos2 = pos(8,:)';
v2 = rec_vel;
lambda2 = lambda;
elapsedtime2 = elapsedtime;
rec_pos2 = rec_pos;
del_r_dot2 = rec_clock_bias_rate;
del_r2 = rec_clock_bias;
ephemeris2 = ephemeris(8,:);
JD_prop_to_2= JD_prop_to;
del_ADR2 = acculumulated_delta_range_derivative(rec_pos2, del_r2,elapsedtime2, v2, del_r_dot2, pos2,lambda2, ephemeris2, JD_prop_to_2);

del_r2==del_r1 
del_r_dot2 == del_r_dot1
rec_pos2 == rec_pos1
elapsedtime2 == elapsedtime1
lambda2 == lambda1
v2 == v1
pos2' == pos1
tf = isequaln(ephemeris1,ephemeris2)
JD_prop_to_2== JD_prop_to_1;
del_ADR2==del_ADR1 




