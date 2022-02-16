
load('starlink_ephemeris_altitude.mat')


% We use heavens above to find all the satellites visible and we collect
% the doppler shift
% We collect the relevant shifts



shift(1,1) = satrec.STARLINK1339.doppler_shift(end_t,1);
shift(2,1) = satrec.STARLINK2170.doppler_shift(end_t,1);
shift(3,1) = satrec.STARLINK2712.doppler_shift(end_t,1);
shift(4,1) = satrec.STARLINK2549.doppler_shift(end_t,1);
shift(5,1) = satrec.STARLINK1488.doppler_shift(end_t,1);
shift(6,1) = satrec.STARLINK1578.doppler_shift(end_t,1);
shift(7,1) = satrec.STARLINK2722.doppler_shift(end_t,1);
shift(8,1) = satrec.STARLINK1554.doppler_shift(end_t,1);
%shift(9,1) = satrec.STARLINK2129.doppler_shift(end_t,1);

elevation(1,1) = satrec.STARLINK1339.elevation(end_t,1);
elevation(2,1) = satrec.STARLINK2170.elevation(end_t,1);
elevation(3,1) = satrec.STARLINK2712.elevation(end_t,1);
elevation(4,1) = satrec.STARLINK2549.elevation(end_t,1);
elevation(5,1) = satrec.STARLINK1488.elevation(end_t,1);
elevation(6,1) = satrec.STARLINK1578.elevation(end_t,1);
elevation(7,1) = satrec.STARLINK2722.elevation(end_t,1);
elevation(8,1) = satrec.STARLINK1554.elevation(end_t,1);
%elevation(9,1) = satrec.STARLINK2129.elevation(end_t,1);

pos(1,:) = satrec.STARLINK1339.pos(end_t,:);
pos(2,:) = satrec.STARLINK2170.pos(end_t,:);
pos(3,:) = satrec.STARLINK2712.pos(end_t,:);
pos(4,:) = satrec.STARLINK2549.pos(end_t,:);
pos(5,:) = satrec.STARLINK1488.pos(end_t,:);
pos(6,:) = satrec.STARLINK1578.pos(end_t,:);
pos(7,:) = satrec.STARLINK2722.pos(end_t,:);
pos(8,:) = satrec.STARLINK1554.pos(end_t,:);
%pos(9,:) = satrec.STARLINK2129.pos(end_t,:);


vel(1,:) = satrec.STARLINK1339.vel(end_t,:);
vel(2,:) = satrec.STARLINK2170.vel(end_t,:);
vel(3,:) = satrec.STARLINK2712.vel(end_t,:);
vel(4,:) = satrec.STARLINK2549.vel(end_t,:);
vel(5,:) = satrec.STARLINK1488.vel(end_t,:);
vel(6,:) = satrec.STARLINK1578.vel(end_t,:);
vel(7,:) = satrec.STARLINK2722.vel(end_t,:);
vel(8,:) = satrec.STARLINK1554.vel(end_t,:);
%vel(9,:) = satrec.STARLINK2129.vel(end_t,:);

ephemeris(1,1) = satrec.STARLINK1339;
ephemeris(2,1) = satrec.STARLINK2170;
ephemeris(3,1) = satrec.STARLINK2712;
ephemeris(4,1) = satrec.STARLINK2549;
ephemeris(5,1) = satrec.STARLINK1488;
ephemeris(6,1) = satrec.STARLINK1578;
ephemeris(7,1) = satrec.STARLINK2722;
ephemeris(8,1) = satrec.STARLINK1554;
%ephemeris(9,1) = satrec.STARLINK2129;



% symbolic check that elevation is positive
for i=1:length(shift)
if(elevation(i,1) <0 )
    display("there is a problem")
end
end


% Goal to calculate GDOP and position
current_time = elapsedtime;

[state ~, Error, delta_y_0,error_0,A_0] = doppler_shift_positioning(shift, pos, lambda, current_time, ephemeris, JD_prop_to);
vpa(state)

figure(1000)
plot(Error)
xlabel('iterations')
ylabel('error')
title('error vs. iterations')

