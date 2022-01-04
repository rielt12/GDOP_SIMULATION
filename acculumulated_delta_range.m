function ADR = acculumulated_delta_range(r, del_R,t_R,pos, lambda)
c=3e8;
del_ion=0;
del_trop = 0;
bias_beat = 0;
sat_clock_offset =0;
t_R = 0;
t_prop = norm(pos-r)/c; % approximation of propagation time
A = earth_rot_mat(t_prop);
ADR = sqrt(dot(r-A*pos,r-A*pos))+c*(del_R-sat_clock_offset)+c*(del_trop-del_ion)+lambda*bias_beat; % approximation of when satellite time is taken: currently sat time=rec_time
end