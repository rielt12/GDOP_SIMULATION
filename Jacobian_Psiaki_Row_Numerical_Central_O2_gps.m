function A  = Jacobian_Psiaki_Row_Numerical_Central_O2_gps(r, c_del_R,v, c_del_R_rate, shift, pos, lambda, current_time, ephemeris, svids)
c=3e8;
A = zeros(1,8);
h =1;





r_x = r;
r_x(1,1) = r(1,1)-h;
y_i = [r_x; c_del_R; v; c_del_R_rate];
del_ADR_1 = compute_cost_doppler_gps(y_i , shift, pos, lambda, current_time, ephemeris, svids);
r_x = r;
r_x(1,1) = r(1,1)+h;
y_i = [r_x; c_del_R; v; c_del_R_rate];
del_ADR_2 = compute_cost_doppler_gps(y_i , shift, pos, lambda, current_time, ephemeris, svids);
r_x = r;
r_x(1,1) = r(1,1)+2*h;
y_i = [r_x; c_del_R; v; c_del_R_rate];
del_ADR_3 = compute_cost_doppler_gps(y_i , shift, pos, lambda, current_time, ephemeris, svids);
r_x = r;
r_x(1,1) = r(1,1)-2*h;
y_i = [r_x; c_del_R; v; c_del_R_rate];
del_ADR_4 = compute_cost_doppler_gps(y_i , shift, pos, lambda, current_time, ephemeris, svids);

A(1,1) =(-del_ADR_3+8*del_ADR_2-8*del_ADR_1+del_ADR_4)/(12*h);



r_y = r;
r_y(2,1) = r(2,1)-h;
y_i = [r_y; c_del_R; v; c_del_R_rate];
del_ADR_1 = compute_cost_doppler_gps(y_i , shift, pos, lambda, current_time, ephemeris, svids);
r_y = r;
r_y(2,1) = r(2,1)+h;
y_i = [r_y; c_del_R; v; c_del_R_rate];
del_ADR_2 = compute_cost_doppler_gps(y_i , shift, pos, lambda, current_time, ephemeris, svids);
r_y = r;
r_y(2,1) = r(2,1)+2*h;
y_i = [r_y; c_del_R; v; c_del_R_rate];
del_ADR_3 = compute_cost_doppler_gps(y_i , shift, pos, lambda, current_time, ephemeris, svids);
r_y = r;
r_y(2,1) = r(2,1)-2*h;
y_i = [r_y; c_del_R; v; c_del_R_rate];
del_ADR_4 = compute_cost_doppler_gps(y_i , shift, pos, lambda, current_time, ephemeris, svids);

A(1,2) =(-del_ADR_3+8*del_ADR_2-8*del_ADR_1+del_ADR_4)/(12*h);


r_z = r;
r_z(3,1) = r(3,1)-h;
y_i = [r_z; c_del_R; v; c_del_R_rate];
del_ADR_1 = compute_cost_doppler_gps(y_i , shift, pos, lambda, current_time, ephemeris, svids);
r_z = r;
r_z(3,1) = r(3,1)+h;
y_i = [r_z; c_del_R; v; c_del_R_rate];
del_ADR_2 = compute_cost_doppler_gps(y_i , shift, pos, lambda, current_time, ephemeris, svids);
r_z = r;
r_z(3,1) = r(3,1)+2*h;
y_i = [r_z; c_del_R; v; c_del_R_rate];
del_ADR_3 = compute_cost_doppler_gps(y_i , shift, pos, lambda, current_time, ephemeris, svids);
r_z = r;
r_z(3,1) = r(3,1)-2*h;
y_i = [r_z; c_del_R; v; c_del_R_rate];
del_ADR_4 = compute_cost_doppler_gps(y_i , shift, pos, lambda, current_time, ephemeris, svids);

A(1,3) =(-del_ADR_3+8*del_ADR_2-8*del_ADR_1+del_ADR_4)/(12*h);


del_R_pert = c_del_R - h;
y_i = [r; del_R_pert; v; c_del_R_rate];
del_ADR_1 = compute_cost_doppler_gps(y_i , shift, pos, lambda, current_time, ephemeris, svids);
del_R_pert = c_del_R + h;
y_i = [r; del_R_pert; v; c_del_R_rate];
del_ADR_2 = compute_cost_doppler_gps(y_i , shift, pos, lambda, current_time, ephemeris, svids);
del_R_pert = c_del_R + 2*h;
y_i = [r; del_R_pert; v; c_del_R_rate];
del_ADR_3 = compute_cost_doppler_gps(y_i , shift, pos, lambda, current_time, ephemeris, svids);
del_R_pert = c_del_R - 2*h;
y_i = [r; del_R_pert; v; c_del_R_rate];
del_ADR_4 = compute_cost_doppler_gps(y_i , shift, pos, lambda, current_time, ephemeris, svids);
A(1,4) =(-del_ADR_3+8*del_ADR_2-8*del_ADR_1+del_ADR_4)/(12*h);


v_x = v;
v_x(1,1) = v(1,1)-h;
y_i = [r; c_del_R; v_x; c_del_R_rate];
del_ADR_1 = compute_cost_doppler_gps(y_i , shift, pos, lambda, current_time, ephemeris, svids);
v_x = v;
v_x(1,1) = v(1,1)+h;
y_i = [r; c_del_R; v_x; c_del_R_rate];
del_ADR_2 = compute_cost_doppler_gps(y_i , shift, pos, lambda, current_time, ephemeris, svids);
v_x = v;
v_x(1,1) = v(1,1)+2*h;
y_i = [r; c_del_R; v_x; c_del_R_rate];
del_ADR_3 = compute_cost_doppler_gps(y_i , shift, pos, lambda, current_time, ephemeris, svids);
v_x = v;
v_x(1,1) = v(1,1)-2*h;
y_i = [r; c_del_R; v_x; c_del_R_rate];
del_ADR_4 = compute_cost_doppler_gps(y_i , shift, pos, lambda, current_time, ephemeris, svids);

A(1,5) =(-del_ADR_3+8*del_ADR_2-8*del_ADR_1+del_ADR_4)/(12*h);


v_y = v;
v_y(2,1) = v(2,1)-h;
y_i = [r; c_del_R; v_y; c_del_R_rate];
del_ADR_1 = compute_cost_doppler_gps(y_i , shift, pos, lambda, current_time, ephemeris, svids);
v_y = v;
v_y(2,1) = v(2,1)+h;
y_i = [r; c_del_R; v_y; c_del_R_rate];
del_ADR_2 = compute_cost_doppler_gps(y_i , shift, pos, lambda, current_time, ephemeris, svids);
v_y = v;
v_y(2,1) = v(2,1)+2*h;
y_i = [r; c_del_R; v_y; c_del_R_rate];
del_ADR_3 = compute_cost_doppler_gps(y_i , shift, pos, lambda, current_time, ephemeris, svids);
v_y = v;
v_y(2,1) = v(2,1)-2*h;
y_i = [r; c_del_R; v_y; c_del_R_rate];
del_ADR_4 = compute_cost_doppler_gps(y_i , shift, pos, lambda, current_time, ephemeris, svids);
A(1,6) =(-del_ADR_3+8*del_ADR_2-8*del_ADR_1+del_ADR_4)/(12*h);


v_z = v;
v_z(3,1) = v(3,1)-h;
y_i = [r; c_del_R; v_z; c_del_R_rate];
del_ADR_1 = compute_cost_doppler_gps(y_i , shift, pos, lambda, current_time, ephemeris, svids);
v_z = v;
v_z(3,1) = v(3,1)+h;
y_i = [r; c_del_R; v_z; c_del_R_rate];
del_ADR_2 = compute_cost_doppler_gps(y_i , shift, pos, lambda, current_time, ephemeris, svids);
v_z = v;
v_z(3,1) = v(3,1)+2*h;
y_i = [r; c_del_R; v_z; c_del_R_rate];
del_ADR_3 = compute_cost_doppler_gps(y_i , shift, pos, lambda, current_time, ephemeris, svids);
v_z = v;
v_z(3,1) = v(3,1)-2*h;
y_i = [r; c_del_R; v_z; c_del_R_rate];
del_ADR_4 = compute_cost_doppler_gps(y_i , shift, pos, lambda, current_time, ephemeris, svids);

A(1,7) =(-del_ADR_3+8*del_ADR_2-8*del_ADR_1+del_ADR_4)/(12*h);




del_R_rate_pert = c_del_R_rate-h;
y_i = [r; c_del_R; v; del_R_rate_pert];
del_ADR_1 = compute_cost_doppler_gps(y_i , shift, pos, lambda, current_time, ephemeris, svids);
del_R_rate_pert = c_del_R_rate+h;
y_i = [r; c_del_R; v; del_R_rate_pert];
del_ADR_2 = compute_cost_doppler_gps(y_i , shift, pos, lambda, current_time, ephemeris, svids);
del_R_rate_pert = c_del_R_rate+2*h;
y_i = [r; c_del_R; v; del_R_rate_pert];
del_ADR_3 = compute_cost_doppler_gps(y_i , shift, pos, lambda, current_time, ephemeris, svids);
del_R_rate_pert = c_del_R_rate-2*h;
y_i = [r; c_del_R; v; del_R_rate_pert];
del_ADR_4 = compute_cost_doppler_gps(y_i , shift, pos, lambda, current_time, ephemeris, svids);
A(1,8) =(-del_ADR_3+8*del_ADR_2-8*del_ADR_1+del_ADR_4)/(12*h);





end