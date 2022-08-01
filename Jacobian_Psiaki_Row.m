function A  = Jacobian_Psiaki_Row(r, c_del_R,t_R, v, c_del_R_rate, pos,lambda, ephemeris,JD_prop_to)

c= 3e8;
del_t_R = 0.1;

del_R = c_del_R/c;
del_R_rate = c_del_R_rate/c;

A = zeros(1,8);

[ADR_T1,vec1,unit1,vel_dot_unit1] = acculumulated_delta_range_T(r-(2*v*del_t_R)/(1+del_R_rate),c*(del_R-(2*del_R_rate*del_t_R)/(1+del_R_rate)),(t_R-2*del_t_R),pos, lambda,ephemeris, JD_prop_to);
[ADR_T2,vec2,unit2,vel_dot_unit2] = acculumulated_delta_range_T(r-(v*del_t_R)/(1+del_R_rate),c*(del_R-(del_R_rate*del_t_R)/(1+del_R_rate)),(t_R-del_t_R),pos, lambda,ephemeris, JD_prop_to);
[ADR_T3,vec3,unit3,vel_dot_unit3] = acculumulated_delta_range_T(r+(v*del_t_R)/(1+del_R_rate),c*(del_R+(del_R_rate*del_t_R)/(1+del_R_rate)),(t_R+del_t_R),pos, lambda, ephemeris, JD_prop_to);
[ADR_T4,vec4,unit4,vel_dot_unit4] = acculumulated_delta_range_T(r+(2*v*del_t_R)/(1+del_R_rate),c*(del_R+(2*del_R_rate*del_t_R)/(1+del_R_rate)),(t_R+2*del_t_R),pos, lambda,ephemeris,JD_prop_to);

A(1,1)= -lambda*(unit1(1,1)-8*unit2(1,1)+8*unit3(1,1)+unit4(1,1))*(1/(12*del_t_R));
A(1,2)= -lambda*(unit1(1,2)-8*unit2(1,2)+8*unit3(1,2)+unit4(1,2))*(1/(12*del_t_R));
A(1,3)= -lambda*(unit1(1,3)-8*unit2(1,3)+8*unit3(1,3)+unit4(1,3))*(1/(12*del_t_R));

[AN4] =-lambda*(vel_dot_unit1-8*vel_dot_unit2+8*vel_dot_unit3+vel_dot_unit4)*(1/(12*c*del_t_R));
A(1,4)=AN4;


A(1,5)=-lambda*(-2*unit1(1,1)+8*unit2(1,1)+8*unit3(1,1)+2*unit4(1,1))*(1/(12));
A(1,6)=-lambda*(-2*unit1(1,2)+8*unit2(1,2)+8*unit3(1,2)+2*unit4(1,2))*(1/(12));
A(1,7)=-lambda*(-2*unit1(1,3)+8*unit2(1,3)+8*unit3(1,3)+2*unit4(1,3))*(1/(12));

A(1,8) = lambda*1;


end

