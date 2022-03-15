function Kost = compute_cost_doppler(y_i , shift, pos, lambda,ephemeris,JD_prop_to)

rec_pos = y_i(1:3,1);
rec_clock_bias = y_i(4,1);
rec_vel = y_i(5:7,1);
rec_clock_bias_rate = y_i(8);


for i=1:length(shift)
    elapsedtime=(JD_prop_to-ephemeris(i).jdsatepoch)*24*60;
    del_ADR(i) = acculumulated_delta_range_derivative(rec_pos, rec_clock_bias,elapsedtime, rec_vel, rec_clock_bias_rate, pos(i,:)',lambda,ephemeris(i,1), JD_prop_to);
    R(i,i) = (0.01)^2;
    error(i)= lambda*(shift(i)-del_ADR(i));
end

Kost = 0.5*norm(error)^2;

end