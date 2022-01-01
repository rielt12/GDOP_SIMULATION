function [jday,jday0,UT]=Julian(y,m,d,h,min,s)

jday0=367*y-floor(7*(y+floor((m+9)/12))/4)+floor(275*m/9)+d+1721013.5;

UT=h+min/60+s/3600;

jday=jday0+UT/24;
end