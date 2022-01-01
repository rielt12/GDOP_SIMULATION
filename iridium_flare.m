function [angular_difference]=iridium_flare(EQ_X,EQ_Y,EQ_Z,V_X, V_Y, V_Z,TOPO_X,TOPO_Y,TOPO_Z,SUN,ROT_1, ROT_2)
% algorithm based off of: http://home.comcast.net/~skysat/algo.txt
% inputs:
% geocentric location of sat: EQ_X,EQ_Y,EQ_Z
% geocentric velocity of the sat: V_X, V_Y, V_Z
% topocentric position of the sat: TOPO_X, TOPO_Y, TOPO_Z,
% geocentric position of the SUN. 
% two mirror rotation angles: ROT_1, ROT_2

% basic premise of algorithm is that if a ray from an observer on earth to
% the sat is reflected onto the sun then a flare will happen.

% We first make a local frame attached to the Sat.
% the x axis is in the direction of the velocity
% the y axis is in the direction of angular momentum
% the z axis completes the right handed frame

pos=[EQ_X,EQ_Y,EQ_Z];
vel=[V_X,V_Y,V_Z];

XX(1,1)=(1/norm(vel))*(vel(1,1));
XX(2,1)=(1/norm(vel))*(vel(1,2));
XX(3,1)=(1/norm(vel))*(vel(1,3));
temp=cross(pos,vel);

YY(1,1)=(1/norm(temp))*temp(1,1);
YY(2,1)=(1/norm(temp))*temp(1,2);
YY(3,1)=(1/norm(temp))*temp(1,3);

temp=cross(XX,YY);
ZZ(1,1)=(1/norm(temp))*temp(1,1);
ZZ(2,1)=(1/norm(temp))*temp(2,1);
ZZ(3,1)=(1/norm(temp))*temp(3,1);

%% 
% we create a transformation table to go from the sat frame to the
% equatorial frame.
%%
Sat2eq(1,1)=  XX(1,1);
Sat2eq(1,2) = YY(1,1);
Sat2eq(1,3) = ZZ(1,1);
Sat2eq(2,1) = XX(2,1);
Sat2eq(2,2) = YY(2,1);
Sat2eq(2,3) = ZZ(2,1);
Sat2eq(3,1) = XX(3,1);
Sat2eq(3,2) = YY(3,1);
Sat2eq(3,3) = ZZ(3,1);


%%
% We now go from that sat frame to the mirror frame. we rotate first about
% the YY axis by ROT_1 and then subsequently about the new 3 axis by ROT_2.
% We do this for each of the sat frame's coordinate axes. 
%%

TT(1,1)= 1;
TT(2,1)=  0;
TT(3,1)=  0;
NTT(1,1) = TT(1,1) * cosd(ROT_1) -  TT(3,1)  * sind(ROT_1);
NTT(2,1) = TT(2,1);
NTT(3,1) = TT(1,1) * sind(ROT_1) +  TT(3,1)  * cosd(ROT_1);
TT(1,1) = NTT(1,1) * cosd(ROT_2) +  NTT(2,1) * sind(ROT_2);
TT(2,1) =-NTT(1,1) * sind(ROT_2) +  NTT(2,1) * cosd(ROT_2);
TT(3,1)=  NTT(3,1);
XX(1,1) = TT(1,1) * Sat2eq(1, 1)+ TT(2,1) * Sat2eq(1,2) + TT(3,1) * Sat2eq(1,3);
XX(2,1) = TT(1,1) * Sat2eq(2, 1)+ TT(2,1) * Sat2eq(2,2) + TT(3,1) * Sat2eq(2,3);
XX(3,1) = TT(1,1) * Sat2eq(3, 1)+ TT(2,1) * Sat2eq(3,2) + TT(3,1) * Sat2eq(3,3);

TT(1,1) = 0;
TT(2,1) = 1;
TT(3,1)=0;
NTT(1,1) = TT(1,1)* cosd(ROT_1) -  TT(3,1)*sind(ROT_1);
NTT(2,1) = TT(2,1);
NTT(3,1) = TT(1,1)* sind(ROT_1) +  TT(3,1)*cosd(ROT_1);
TT(1,1) = NTT(1,1)* cosd(ROT_2) +  NTT(2,1)*sind(ROT_2);
TT(2,1) =-NTT(1,1)* sind(ROT_2) +  NTT(2,1)*cosd(ROT_2);
TT(3,1)=  NTT(3,1);
YY(1,1) = TT(1,1) * Sat2eq(1,1) + TT(2,1) * Sat2eq(1,2) + TT(3,1) * Sat2eq(1,3);
YY(2,1) = TT(1,1) * Sat2eq(2,1) + TT(2,1) * Sat2eq(2,2) + TT(3,1) * Sat2eq(2,3);
YY(3,1) = TT(1,1) * Sat2eq(3,1) + TT(2,1) * Sat2eq(3,2) + TT(3,1) * Sat2eq(3,3);

TT(1,1) = 0;
TT(2,1) = 0;
TT(3,1)=  1;
NTT(1,1) = TT(1,1) * cosd(ROT_1) -  TT(3,1)  * sind(ROT_1);
NTT(2,1) = TT(2,1);
NTT(3,1) = TT(1,1) * sind(ROT_1) +  TT(3,1)  * cosd(ROT_1);
TT(1,1) = NTT(1,1) * cosd(ROT_2) +  NTT(2,1) * sind(ROT_2);
TT(2,1) =-NTT(1,1) * sind(ROT_2) +  NTT(2,1) * cosd(ROT_2);
TT(3,1)=  NTT(3,1);
ZZ(1,1) = TT(1,1) * Sat2eq(1,1) + TT(2,1) * Sat2eq(1,2) + TT(3,1) * Sat2eq(1,3);
ZZ(2,1) = TT(1,1) * Sat2eq(2,1) + TT(2,1) * Sat2eq(2,2) + TT(3,1) * Sat2eq(2,3);
ZZ(3,1) = TT(1,1) * Sat2eq(3,1) + TT(2,1) * Sat2eq(3,2) + TT(3,1) * Sat2eq(3,3);

Mirror2eq(1, 1) = XX(1,1);
Mirror2eq(1, 2) = YY(1,1);
Mirror2eq(1, 3) = ZZ(1,1);
Mirror2eq(2, 1) = XX(2,1);
Mirror2eq(2, 2) = YY(2,1);
Mirror2eq(2, 3) = ZZ(2,1);
Mirror2eq(3, 1) = XX(3,1);
Mirror2eq(3, 2) = YY(3,1);
Mirror2eq(3, 3) = ZZ(3,1);     


%%
% We now convert the observer-sat vector into the mirror coordinate system.
% We create the reflection by negating the x coordinate. (The x coordinate
% in the mirror system is the normal to the mirror.) We then convert back
% to the equatorial system and find the angle between this vector (the
% reflection) and the direction of the sun. If the angle is very small then
% there is an irridium flare.
%%

TT(1,1) = TOPO_X * Mirror2eq(1,1) + TOPO_Y * Mirror2eq(2,1) + TOPO_Z * Mirror2eq(3, 1);
TT(2,1) = TOPO_X * Mirror2eq(1, 2) +TOPO_Y * Mirror2eq(2, 2) + TOPO_Z * Mirror2eq(3, 2);
TT(3,1) = TOPO_X * Mirror2eq(1, 3) + TOPO_Y * Mirror2eq(2, 3) + TOPO_Z*Mirror2eq(3, 3);

 TT(1,1) = -TT(1,1);

if  TT(1,1)>=0
 SUN_REF_X = TT(1,1) * Mirror2eq(1, 1) + TT(2,1) * Mirror2eq(1, 2) + TT(3,1) * Mirror2eq(1, 3);
 SUN_REF_Y = TT(1,1) * Mirror2eq(2, 1) + TT(2,1) * Mirror2eq(2, 2) + TT(3,1) * Mirror2eq(2, 3);
 SUN_REF_Z = TT(1,1) * Mirror2eq(3, 1) + TT(2,1) * Mirror2eq(3, 2) + TT(3,1) * Mirror2eq(3, 3);
 
 SUN_REF(1,1)=SUN_REF_X;
 SUN_REF(2,1)=SUN_REF_Y;
 SUN_REF(3,1)=SUN_REF_Z;
 
 
 SUN=transpose(SUN);
 angular_difference=Spherdist(SUN_REF, SUN);
 

 
 
else
    angular_difference=180;
end


end