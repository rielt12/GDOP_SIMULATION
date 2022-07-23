function  F= solve_t_prop(t_prop,t_R,del_R, ephemeris,JD_prop_to)
c= 3e8;
%t_prop initial guess 
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

r_prop_ecef(1,1) = r_prop_ecef(1,1);
r_j = norm(r_prop_ecef);

c*t_prop;
F = r_j - c*t_prop;


end