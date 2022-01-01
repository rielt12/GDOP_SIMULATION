function gst=theta_g(y,m,d,h,min,s)

[~,jday0,UT]=Julian(y,m,d,h,min,s);

T0=(jday0-2451545)/36525;

theta_go=100.4606184+36000.77004*T0+0.000387933*T0^2-2.583*10^-8*T0^3;

if theta_go>360
    q=floor(theta_go/360);
    theta_go=theta_go-q*360;
end

gst=theta_go+360.98564724*UT/24;

if gst>360
    gst=360-gst;
end


end