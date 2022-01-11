function [AN4] = del_del_del_r_del_ADR_del_t(r, del_R,t_R, v, del_R_rate, pos,lambda, ephemeris, JD_prop_to)
c=3e8;
t_prop = norm(pos-r)/c; % approximation of propagation time
A = earth_rot_mat(t_prop);

[~, r_prop_eci, v_prop_eci] = sgp4(ephemeris,t_R-del_R/60-t_prop/60);
r_prop_eci = r_prop_eci.*1000;
v_prop_eci = v_prop_eci.*1000; % convert to meters
    
%convert to ecef
JD_prop_to_dt =datetime(JD_prop_to-t_R*(1/60)*(1/24)-del_R*(1/3600)*(1/24)-t_prop*(1/3600)*(1/24),'convertfrom','juliandate');  % convert to days
[y,mo,d] = ymd(JD_prop_to_dt); 
[h,m,s] = hms(JD_prop_to_dt);
time  =[y,mo,d,h,m,s];
[r_prop_ecef,~]=eci2ecef(r_prop_eci',v_prop_eci',time(1),time(2),time(3),time(4),time(5),time(6));

vec = r-A*r_prop_ecef;

% numerical derivative : perturb del_r by 0.01 and propagate and take the difference divided by 0.01
[~, r_prop_eci_pert, v_prop_eci_pert] = sgp4(ephemeris,t_R-(del_R+0.01)/60-t_prop/60);
r_prop_eci_pert = r_prop_eci_pert.*1000;
v_prop_eci_pert = v_prop_eci_pert.*1000; % convert to meters

%convert to ecef
JD_prop_to_dt =datetime(JD_prop_to-t_R*(1/60)*(1/24)-(del_R+0.01)*(1/3600)*(1/24)-t_prop*(1/3600)*(1/24),'convertfrom','juliandate');  % convert to days
[y,mo,d] = ymd(JD_prop_to_dt); 
[h,m,s] = hms(JD_prop_to_dt);
time  =[y,mo,d,h,m,s];
[r_prop_ecef_pert,~]=eci2ecef(r_prop_eci_pert',v_prop_eci_pert',time(1),time(2),time(3),time(4),time(5),time(6));

derivative = (r_prop_ecef_pert-r_prop_ecef)./0.01;

dot_prod = dot(vec,derivative);


del_t_R = 0.1; 

% the c's cancel out in the finite difference

ADR_1 = dot_prod/acculumulated_delta_range_T(r-(2*v*del_t_R)/(1+del_R_rate),del_R-(2*del_R_rate*del_t_R)/(1+del_R_rate),(t_R-2*del_t_R),pos, lambda,ephemeris,JD_prop_to);
ADR_2 = dot_prod/acculumulated_delta_range_T(r-(2*v*del_t_R)/(1+del_R_rate),del_R-(del_R_rate*del_t_R)/(1+del_R_rate),(t_R-del_t_R),pos, lambda,ephemeris, JD_prop_to);
ADR_3 = dot_prod/acculumulated_delta_range_T(r+(v*del_t_R)/(1+del_R_rate),del_R+(del_R_rate*del_t_R)/(1+del_R_rate),(t_R+del_t_R),pos, lambda, ephemeris,JD_prop_to);
ADR_4 = dot_prod/acculumulated_delta_range_T(r+(2*v*del_t_R)/(1+del_R_rate),del_R+(2*del_R_rate*del_t_R)/(1+del_R_rate),(t_R+2*del_t_R),pos, lambda,ephemeris,JD_prop_to);

AN4 = (ADR_1-8*ADR_2+8*ADR_3-ADR_4)/(12*del_t_R);



end