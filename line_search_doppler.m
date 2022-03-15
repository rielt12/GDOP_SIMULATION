function [gamma,nwd] = line_search_doppler(y_i, delta_y, N_g ,shift, pos, lambda,ephemeris,JD_prop_to)

Gamma  = 1;
prev_cost = compute_cost_doppler(y_i , shift, pos, lambda,ephemeris,JD_prop_to);
for j=[1:N_g]
    
Gamma = j/N_g;

y_i(1,1) = y_i(1,1)+Gamma*delta_y(1,1);
y_i(2,1) = y_i(2,1)-Gamma*delta_y(2,1);
% y_i(3,1) = y_i(3,1)-Gamma*delta_y(3,1);
% y_i(4,1) = y_i(4,1)+Gamma*delta_y(4,1);
%  y_i(5,1) = y_i(5,1)+Gamma*delta_y(5,1);
% y_i(6,1) = y_i(6,1)-Gamma*delta_y(6,1);
% y_i(7,1) = y_i(7,1)-Gamma*delta_y(7,1);
  %y_i(8,1) = y_i(8,1)+Gamma*delta_y(8,1);
cost = compute_cost_doppler(y_i , shift, pos, lambda,ephemeris,JD_prop_to);


if(cost<prev_cost)
gamma = Gamma;
nwd =0;
break;
end
if(j == N_g)
  display('no way down')
    nwd = 1;
  gamma = 1;
end


end


end