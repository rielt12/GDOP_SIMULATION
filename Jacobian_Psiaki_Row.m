function A  = Jacobian_Psiaki_Row(r, del_R,t_R, v, del_R_rate, pos,lambda, ephemeris,JD_prop_to)

A = zeros(1,8);

[AN1, AN2, AN3] = del_del_r_del_ADR_del_t(r, del_R,t_R, v, del_R_rate, pos,lambda, ephemeris,JD_prop_to);
A(1,1)=AN1;
A(1,2)=AN2;
A(1,3)=AN3;

[AN4] = del_del_del_r_del_ADR_del_t(r, del_R,t_R, v, del_R_rate, pos,lambda, ephemeris,JD_prop_to);
A(1,4)=AN4;


[AN5,AN6,AN7] = del_ADR_del_v(r, del_R,t_R, v, del_R_rate, pos,lambda, ephemeris, JD_prop_to);
A(1,5:7) = [AN5,AN6,AN7];

AN8 = del_ADR_del_r_dot(r, del_R,t_R, v, del_R_rate, pos,lambda, ephemeris, JD_prop_to);
A(1,8) = AN8;


end

