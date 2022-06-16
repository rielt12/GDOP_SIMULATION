function ADR = acculumulated_delta_range_sat_rec(r, del_R,t_R,pos, lambda,eph,svid)
c=3e8;
del_ion=0;
del_trop = 0;
bias_beat = 0;
sat_clock_offset =0;

t_prop = norm(pos-r)/c; % approximation of propagation time
A = earth_rot_mat(t_prop);

%convert to ecef 
% currently does not propage backwards

sat_pos =get_satpos(t_R-del_R-t_prop,svid,eph,2);
r_prop_ecef = sat_pos(1:3,1);


    
ADR = sqrt(dot(r-A*r_prop_ecef,r-A*r_prop_ecef))+c*(del_R-sat_clock_offset)+c*(del_trop-del_ion)+lambda*bias_beat; 