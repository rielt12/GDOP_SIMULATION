function [state rec_GDOP, Error, delta_y_0,error_0,A_0, del_ADR] = doppler_shift_positioning(shift, pos, lambda, current_time, ephemeris, JD_prop_to)

max_iter = 15000;

% inputs
% shifts are observed doppler shifts
% current time is the time in minutes since epoch (elapsed time) that all sats are in view.
% julian date is JD when sats are in view needed  to convert ECI to ECEF


%constants
c= 3e8;
longitude=-77.3790;
latitude=38.9649; 
H=0; 
%1085027.1960869177710264921188354
%-4845789.2537157507613301277160645
%3989288.0202742042019963264465332
[x,y,z]=lla2ecef_AB(latitude*(2*pi/360),longitude*(pi/180),H); % this function takes east longitude. change if need be.
rec_pos(1,1) =  x+20e3;
rec_pos(2,1) =  y;
rec_pos(3,1) =  z;
%rec_pos(1,1) = 0;
%rec_pos(2,1) = 0;
%rec_pos(3,1) = 0;
rec_clock_bias =0  ;
rec_vel = [0;0;0];
rec_clock_bias_rate= 0;

y_i = [rec_pos; rec_clock_bias; rec_vel; rec_clock_bias_rate];
y_0 =y_i;


t_R =current_time; % minutes

R  =zeros(length(shift),length(shift));
%first compute error terms
for i=1:length(shift)
    elapsedtime=(JD_prop_to-ephemeris(i).jdsatepoch)*24*60;
    del_ADR(i) = acculumulated_delta_range_derivative(rec_pos, rec_clock_bias,elapsedtime, rec_vel, rec_clock_bias_rate, pos(i,:)',lambda,ephemeris(i,1), JD_prop_to);
    R(i,i) = (0.01)^2;
    error(i)= lambda*(shift(i)-del_ADR(i));
end
error_0= error';

% form the jacobian 

% for i=1:length(shift)
% A(i,:) = Jacobian_Psiaki_Row(rec_pos, c_rec_clock_bias/c,elapsedtime, rec_vel, c_rec_clock_bias_rate/c, pos(i),lambda,ephemeris(i), JD_prop_to);
% end
% A_1 = A;

for i=1:length(shift)
A(i,:) = Jacobian_Psiaki_Row_Numerical(rec_pos, rec_clock_bias,elapsedtime, rec_vel, rec_clock_bias_rate, pos(i),lambda,ephemeris(i), JD_prop_to);
end
A_0 = A;


delta_y_0 =inv((A'*A))*A'*error';






Bazooka = 0;
iter = 1;
Error(iter,1) = norm(error);
while(norm(error) >1e-8 && iter < max_iter)

iter  
norm(error)

%delta_y =inv((A'*R*A))*A'*R*error';
%delta_y = lschol(A'*R*A,A'*R*error');
%delta_y = (A'*R*A)\A'*R*error';
%delta_y = lsqr((A'*A),A'*error');
%delta_y =  inv(A'*A)*A'*error';
delta_y = pinv(A)*error';

if(iter ==1)
tau = 1;
end

if(iter>1 && Error(iter)>Error(iter-1) && tau > 1e-3)
%y_i = state_record(iter-2,:)';
y_i = y_0;
tau=tau/2;
Bazooka = Bazooka+1;
end
% if(Bazooka > 15)
%     tau = 1;
%     Bazooka  = 0;
% end

tau
y_i = y_i + tau*delta_y;
state_record(iter,:) = y_i;

%y_i(1,1) = y_i(1,1)+tau*delta_y(1,1);
%y_i(2,1) = y_i(2,1)+tau*delta_y(2,1);
%  y_i(3,1) = y_i(3,1)+tau*delta_y(3,1);
%  y_i(4,1) = y_i(4,1)+tau*delta_y(4,1);
%  y_i(5,1) = y_i(5,1)+tau*delta_y(5,1);
%  y_i(6,1) = y_i(6,1)+tau*delta_y(6,1);
%  y_i(7,1) = y_i(7,1)+tau*delta_y(7,1);
%  y_i(8,1) = y_i(8,1)+tau*delta_y(8,1);



rec_pos =y_i(1:3,1);
rec_clock_bias =y_i(4,1);
rec_vel =y_i(5:7,1);
rec_clock_bias_rate= y_i(8,1);


t_R =current_time; % minutes


R  =zeros(length(shift),length(shift));
%first compute error terms
for i=1:length(shift)
    elapsedtime=(JD_prop_to-ephemeris(i).jdsatepoch)*24*60;
    del_ADR(i) = acculumulated_delta_range_derivative(rec_pos, rec_clock_bias,elapsedtime, rec_vel, rec_clock_bias_rate, pos(i,:)',lambda,ephemeris(i,1), JD_prop_to);
    R(i,i) = (0.01)^2;
    error(i)= lambda*(shift(i)-del_ADR(i));
end
 Error(iter+1,1) = norm(error);


for i=1:length(shift)
 A(i,:) = Jacobian_Psiaki_Row_Numerical(rec_pos, rec_clock_bias,elapsedtime, rec_vel, rec_clock_bias_rate, pos(i),lambda,ephemeris(i), JD_prop_to);
end


% plotting shifts
figure(400)
clf
num=1:9;
scatter(num,shift,'MarkerFaceColor','blue')
hold on
scatter(num,del_ADR,'MarkerFaceColor','red')
xlabel('satellite number')
ylabel('ADR Derivatives')
legend('measured','estimated')

iter = iter+1;

end




rec_GDOP = sqrt(trace(inv(A'*A)));
state = y_i;


end