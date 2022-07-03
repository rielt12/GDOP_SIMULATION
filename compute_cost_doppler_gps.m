function Kost = compute_cost_doppler_gps(y_i , shift, pos, lambda, current_time, ephemeris, svids)
c=3e8;

rec_pos = y_i(1:3,1);
rec_clock_bias = y_i(4,1)/c;
rec_vel = y_i(5:7,1);
rec_clock_bias_rate = y_i(8,1)/c;


bias = [ 
106.56481186797157079126918688416
  106.106679178774356842041015625
108.28082670271470533407409675419
108.43689345071717866630933713168
106.06789046277521038064151071012
107.72690038631384368272847495973
106.26553490012906877382192760706
 105.1043247431520626378187444061
 106.40465396145975773833924904466
];

m_b = mean(bias);

for i=1:length(shift)
    del_ADR(i) = acculumulated_delta_range_derivative_gps(rec_pos, rec_clock_bias,current_time, rec_vel, rec_clock_bias_rate, pos(i,:)',lambda,ephemeris, svids(i));
    R(i,i) = (0.01)^2;
    error(i)= lambda*(del_ADR(i)-bias(i)-shift(i));
end



Kost = error;

end