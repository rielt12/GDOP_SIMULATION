function [gamma,nwd] = line_search_doppler(y_i, delta_y, N_g ,shift, pos, lambda,ephemeris,JD_prop_to)

Gamma  = 1;
prev_cost = compute_cost_doppler(y_i , shift, pos, lambda,ephemeris,JD_prop_to);
prev_cost = 0.5*norm(prev_cost)^2;
for j=[N_g:-1:1]
    
Gamma = j/N_g;


cost = compute_cost_doppler(y_i -Gamma*delta_y, shift, pos, lambda,ephemeris,JD_prop_to);
cost = 0.5*norm(cost)^2;

if(cost<prev_cost)
gamma = Gamma;
nwd =0;
break;
end
if(j == 1)
  display('no way down')
    nwd = 1;
  gamma = 1;
end


end


end