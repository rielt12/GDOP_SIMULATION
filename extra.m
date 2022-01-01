% code to compute azimuth and elevation

[x,y,z]=lla2ecef_AB(latitude*(2*pi/360),-(longitude*(pi/180)),H); % this function takes east longitude. change if need be.
ecef_ref=[x/1000;y/1000;z/1000];
eci_ref=ECEFtoECI(JD2,ecef_ref); % reference (ground station) location in ECI
eci_ref=transpose(eci_ref);
GMST = JD2GMST(JD2);
LST = GMST-longitude;
[range,az,el]=eci2rangeazel(r.*1000,eci_ref.*1000 ,LST, latitude);

S.a = 'test';
S.b = 1;
S.c = rand(2);
fNames = fieldnames(S);
for n = 1:length(fNames)
    disp(S.(fNames{n}))
end


Rel_Times = {};
for n = 1:length(combos(1e6,:))
 Rel_Times(n,:) = {satrec.(fNames{combos(1,n)}).relevant_time};  
end

%mintersect(Rel_Times{:}')
mintersect(satrec.(fNames{1}).relevant_time,satrec.(fNames{2}).relevant_time, satrec.(fNames{3}).relevant_time,satrec.(fNames{4}).relevant_time, satrec.(fNames{5}).relevant_time, satrec.(fNames{6}).relevant_time, satrec.(fNames{7}).relevant_time, satrec.(fNames{8}).relevant_time...
 satrec.(fNames{9}).relevant_time, satrec.(fNames{10}).relevant_time, satrec.(fNames{11}).relevant_time, satrec.(fNames{12}).relevant_time, satrec.(fNames{13}).relevant_time)
