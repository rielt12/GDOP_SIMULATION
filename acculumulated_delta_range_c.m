function ADR = acculumulated_delta_range_c(r, cdel_R,t_R,pos, lambda, ephemeris, JD_prop_to)
c=3e8;
del_ion=0;
del_trop = 0;
bias_beat = 0;
sat_clock_offset =0;

t_prop = norm(pos-r)/c; % approximation of propagation time
A = earth_rot_mat(t_prop);

del_R = cdel_R/c;



[~, r_prop_eci, v_prop_eci] = sgp4(ephemeris,t_R-(cdel_R/c)/60-t_prop/60);
r_prop_eci = r_prop_eci.*1000;
v_prop_eci = v_prop_eci.*1000; % convert to meters
    
%convert to ecef
JD_prop_to_dt =datetime(JD_prop_to-t_R*(1/60)*(1/24)-(cdel_R/c)*(1/3600)*(1/24)-t_prop*(1/3600)*(1/24),'convertfrom','juliandate');  % convert to days
[y,mo,d] = ymd(JD_prop_to_dt); 
[h,m,s] = hms(JD_prop_to_dt);
time  =[y,mo,d,h,m,s];
[r_prop_ecef,~]=eci2ecef(r_prop_eci',v_prop_eci',time(1),time(2),time(3),time(4),time(5),time(6));
  

ADR = sqrt(dot(r-A*r_prop_ecef,r-A*r_prop_ecef))+c*del_R-c*sat_clock_offset+c*(del_trop-del_ion)+lambda*bias_beat; 