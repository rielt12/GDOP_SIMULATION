function [r_ecef,v_ecef]=eci2ecef(r,v,y,m,d,h,min,s)

gst=theta_g(y,m,d,h,min,s);

if gst>360
    gst=gst-360;
end
    

rotation_matrix(1,1)=cosd(gst);
rotation_matrix(2,1)=sind(gst);
rotation_matrix(3,1)=0;
rotation_matrix(1,2)=-sind(gst);
rotation_matrix(2,2)=cosd(gst);
rotation_matrix(3,2)=0;
rotation_matrix(1,3)=0;
rotation_matrix(2,3)=0;
rotation_matrix(3,3)=1;

r_ecef=transpose(rotation_matrix)*r;
v_ecef=transpose(rotation_matrix)*v;

end
