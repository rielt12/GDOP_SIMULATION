function gamma = line_search(x_prev, delta_x, N_g ,y)

Gamma  = 0;
prev_cost = compute_cost(x_prev,y);
for j=[1:N_g]
    j;
Gamma = j/N_g;
cost = compute_cost(x_prev-Gamma*delta_x,y);


if(cost<prev_cost)
gamma = Gamma;
break;
end

end


end