function A = earth_rot_mat(del_t_prop)
w_e = 7.2921150e-5; 
A = [cos(w_e*del_t_prop), sin(w_e*del_t_prop), 0;
     -sin(w_e*del_t_prop), cos(w_e*del_t_prop), 0;
      0, 0,1];
end