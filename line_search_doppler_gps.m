function [gamma,nwd] = line_search_doppler_gps(y_i, delta_y, N_g ,shift, pos, lambda, current_time, ephemeris, svids)

Gamma  = 1;
prev_cost = compute_cost_doppler_gps(y_i , shift, pos, lambda, current_time, ephemeris, svids);
prev_cost = 0.5*norm(prev_cost)^2;
for j=[N_g:-1:1]
 %for j =[1:N_g]

Gamma = j/N_g;
c = 3e8;

delta_y(4,1) = delta_y(4,1)/c; % helps with convergence
delta_y(8,1) = delta_y(8,1);
%y_i = y_i -Gamma*delta_y;
y_i(1) = y_i(1)-Gamma*delta_y(1,1);
y_i(2) = y_i(2)-Gamma*delta_y(2,1);
y_i(3) = y_i(3)-Gamma*delta_y(3,1);
 y_i(4) = y_i(4)-Gamma*delta_y(4,1);
 y_i(5) = y_i(5)-Gamma*delta_y(5,1);
 y_i(6) = y_i(6)-Gamma*delta_y(6,1);
 y_i(7) = y_i(7)-Gamma*delta_y(7,1);
 y_i(8) = y_i(8)-Gamma*delta_y(8,1);
y_i(4:8,1)  = [0;0;0;0;0]; % static condition % and no clockbias
%y_i(5:7) = [0;0;0]; % static condition

new = y_i;
cost = compute_cost_doppler_gps(new, shift, pos, lambda, current_time, ephemeris, svids);
cost = 0.5*norm(cost)^2;

if(cost<prev_cost)
gamma = Gamma;
nwd =0;
break;
end


if(j==1)
  display('no way down')
  gamma = 1;  
  nwd = 1;
end



if(cost>prev_cost)
continue;
end





end


end