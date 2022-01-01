function [range,az,el]=eci2rangeazel(eci,site ,lst, latitude)
range=eci-site;



Rx(1,1)=1;
Rx(1,2)=0;
Rx(1,3)=0;
Rx(2,1)=0;
Rx(2,2)=cosd(90-latitude);
Rx(2,3)=sind(90-latitude);
Rx(3,1)=0;
Rx(3,2)=-sind(90-latitude);
Rx(3,3)=cosd(90-latitude);

Rz(1,1)=cosd(90+lst);
Rz(1,2)=sind(90+lst);
Rz(1,3)=0;
Rz(2,1)=-sind(90+lst);
Rz(2,2)=cosd(90+lst);
Rz(2,3)=0;
Rz(3,1)=0;
Rz(3,2)=0;
Rz(3,3)=1;

range_enz=Rx*Rz*transpose(range);
az=atan2d(range_enz(1,1),range_enz(2,1));
el=atan2d(range_enz(3,1),sqrt(range_enz(1,1)^2+range_enz(2,1)^2));


range=norm(range_enz);
end
