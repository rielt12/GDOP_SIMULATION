function A = earth_rot_mat_deriv(del_t_prop)
w_e = 7.2921150e-5; 
A = [-w_e*sin(w_e*del_t_prop), w_e*cos(w_e*del_t_prop), 0;
     -w_e*cos(w_e*del_t_prop), -w_e*sin(w_e*del_t_prop), 0;
      0, 0,0];
end