function [AN1, AN2, AN3] = del_del_r_del_ADR_del_t(r, del_R,t_R, v, del_R_rate, pos,lambda, ephemeris, JD_prop_to)
c=3e8;
t_prop = norm(pos-r)/c; % approximation of propagation time
A = earth_rot_mat(t_prop);

[~, r_prop_eci, v_prop_eci] = sgp4(ephemeris,t_R-del_R/60-t_prop/60);
r_prop_eci = r_prop_eci.*1000;
v_prop_eci = v_prop_eci.*1000; % convert to meters
    
%convert to ecef
JD_prop_to_dt =datetime(JD_prop_to-t_R*(1/60)*(1/24)-del_R*(1/3600)*(1/24)-t_prop*(1/3600)*(1/24),'convertfrom','juliandate');  % convert minutes to days
[y,mo,d] = ymd(JD_prop_to_dt); 
[h,m,s] = hms(JD_prop_to_dt);
time  =[y,mo,d,h,m,s];
[r_prop_ecef,~]=eci2ecef(r_prop_eci',v_prop_eci',time(1),time(2),time(3),time(4),time(5),time(6));

vec = r-A*r_prop_ecef;

del_t_R = 0.1; 

ADR_1 = vec(1,1)/acculumulated_delta_range_T(r-(2*v*del_t_R)/(1+del_R_rate),del_R-(2*del_R_rate*del_t_R)/(1+del_R_rate),(t_R-2*del_t_R),pos, lambda,ephemeris, JD_prop_to);
ADR_2 = vec(1,1)/acculumulated_delta_range_T(r-(2*v*del_t_R)/(1+del_R_rate),del_R-(del_R_rate*del_t_R)/(1+del_R_rate),(t_R-del_t_R),pos, lambda,ephemeris, JD_prop_to);
ADR_3 = vec(1,1)/acculumulated_delta_range_T(r+(v*del_t_R)/(1+del_R_rate),del_R+(del_R_rate*del_t_R)/(1+del_R_rate),(t_R+del_t_R),pos, lambda, ephemeris, JD_prop_to);
ADR_4 = vec(1,1)/acculumulated_delta_range_T(r+(2*v*del_t_R)/(1+del_R_rate),del_R+(2*del_R_rate*del_t_R)/(1+del_R_rate),(t_R+2*del_t_R),pos, lambda,ephemeris, JD_prop_to);

AN1 = (ADR_1-8*ADR_2+8*ADR_3-ADR_4)/(12*del_t_R);

ADR_1 = vec(2,1)/acculumulated_delta_range_T(r-(2*v*del_t_R)/(1+del_R_rate),del_R-(2*del_R_rate*del_t_R)/(1+del_R_rate),(t_R-2*del_t_R),pos, lambda,ephemeris,JD_prop_to);
ADR_2 = vec(2,1)/acculumulated_delta_range_T(r-(2*v*del_t_R)/(1+del_R_rate),del_R-(del_R_rate*del_t_R)/(1+del_R_rate),(t_R-del_t_R),pos, lambda,ephemeris,JD_prop_to);
ADR_3 = vec(2,1)/acculumulated_delta_range_T(r+(v*del_t_R)/(1+del_R_rate),del_R+(del_R_rate*del_t_R)/(1+del_R_rate),(t_R+del_t_R),pos, lambda, ephemeris, JD_prop_to);
ADR_4 = vec(2,1)/acculumulated_delta_range_T(r+(2*v*del_t_R)/(1+del_R_rate),del_R+(2*del_R_rate*del_t_R)/(1+del_R_rate),(t_R+2*del_t_R),pos, lambda,ephemeris, JD_prop_to);

AN2 = (ADR_1-8*ADR_2+8*ADR_3-ADR_4)/(12*del_t_R);


ADR_1 = vec(3,1)/acculumulated_delta_range_T(r-(2*v*del_t_R)/(1+del_R_rate),del_R-(2*del_R_rate*del_t_R)/(1+del_R_rate),(t_R-2*del_t_R),pos, lambda,ephemeris,JD_prop_to);
ADR_2 = vec(3,1)/acculumulated_delta_range_T(r-(2*v*del_t_R)/(1+del_R_rate),del_R-(del_R_rate*del_t_R)/(1+del_R_rate),(t_R-del_t_R),pos, lambda,ephemeris,JD_prop_to);
ADR_3 = vec(3,1)/acculumulated_delta_range_T(r+(v*del_t_R)/(1+del_R_rate),del_R+(del_R_rate*del_t_R)/(1+del_R_rate),(t_R+del_t_R),pos, lambda, ephemeris,JD_prop_to);
ADR_4 = vec(3,1)/acculumulated_delta_range_T(r+(2*v*del_t_R)/(1+del_R_rate),del_R+(2*del_R_rate*del_t_R)/(1+del_R_rate),(t_R+2*del_t_R),pos, lambda,ephemeris,JD_prop_to);

AN3 = (ADR_1-8*ADR_2+8*ADR_3-ADR_4)/(12*del_t_R);

end