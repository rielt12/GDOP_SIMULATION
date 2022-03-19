function A  = Jacobian_Psiaki_Row_Approx(r, del_R,t_R, v, del_R_rate, pos, shift,lambda, ephemeris,JD_prop_to)

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
[r_prop_ecef,v_prop_ecef]=eci2ecef(r_prop_eci',v_prop_eci',time(1),time(2),time(3),time(4),time(5),time(6));

rho = r-A*r_prop_ecef./norm(r-A*r_prop_ecef);

% numerical derivative : perturb t_r by 1e-6 and propagate and take the difference divided by 0.01
[~, r_prop_eci_pert, v_prop_eci_pert] = sgp4(ephemeris,t_R+1e-6-(del_R)/60-t_prop/60);
r_prop_eci_pert = r_prop_eci_pert.*1000;
v_prop_eci_pert = v_prop_eci_pert.*1000; % convert to meters

%convert to ecef
JD_prop_to_dt =datetime(JD_prop_to-t_R+1e-6*(1/60)*(1/24)-(del_R)*(1/3600)*(1/24)-t_prop*(1/3600)*(1/24),'convertfrom','juliandate');  % convert to days
[y,mo,d] = ymd(JD_prop_to_dt); 
[h,m,s] = hms(JD_prop_to_dt);
time  =[y,mo,d,h,m,s];
[r_prop_ecef_pert,v_prop_ecef_pert]=eci2ecef(r_prop_eci_pert',v_prop_eci_pert',time(1),time(2),time(3),time(4),time(5),time(6));

rho1 = r-A*r_prop_ecef_pert./norm(r-A*r_prop_ecef_pert);

rho_dot = (rho1-rho)./1e-6;
v_dot = v_prop_ecef_pert-v_prop_ecef./1e-6;


clear A;
A(1,1)=rho_dot(1,1);
A(1,2)=rho_dot(2,1);
A(1,3)=rho_dot(3,1);


A(1,4)=dot(rho,v_dot)+dot(rho_dot,v);


A(1,5) = rho(1,1);
A(1,6) = rho(2,1);
A(1,7) = rho(3,1);

A(1,8) = c;


end

