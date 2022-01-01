clear all
close all
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
shift(9,1) = satrec.STARLINK2129.doppler_shift(end_t,1);

elevation(1,1) = satrec.STARLINK1339.elevation(end_t,1);
elevation(2,1) = satrec.STARLINK2170.elevation(end_t,1);
elevation(3,1) = satrec.STARLINK2712.elevation(end_t,1);
elevation(4,1) = satrec.STARLINK2549.elevation(end_t,1);
elevation(5,1) = satrec.STARLINK1488.elevation(end_t,1);
elevation(6,1) = satrec.STARLINK1578.elevation(end_t,1);
elevation(7,1) = satrec.STARLINK2722.elevation(end_t,1);
elevation(8,1) = satrec.STARLINK1554.elevation(end_t,1);
elevation(9,1) = satrec.STARLINK2129.elevation(end_t,1);

% symbolic check that elevation is positive
for i=1:9
    i
if(elevation(i,1) <0 )
display("there is a problem")
end
end


% Goal to calculate GDOP and position



