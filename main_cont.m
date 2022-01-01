clear all
close all
load('ephemeris_altitude.mat')

% find time when at lest 8 elevations are greater than zero.

el_mask = 0;
for n = 1:length(fNames)
    i = 0;
    satrec.(fNames{n}).relevant_time=[];
    for t = 1:end_t
    if(satrec.(fNames{n}).elevation(t,1)  >el_mask)
        i=i+1;
      satrec.(fNames{n}).relevant_time(i,1) = t;
    end
    end
end

combos = nchoosek(1:length(fNames),8);

Rel_Times = {};
for n = 1:length(combos(1e6,:))
 Rel_Times(n,:) = {satrec.(fNames{combos(1,n)}).relevant_time};  
end

mintersect(Rel_Times{:})
%mintersect(satrec.(fNames{1}).relevant_time,satrec.(fNames{2}).relevant_time, satrec.(fNames{3}).relevant_time,satrec.(fNames{4}).relevant_time, satrec.(fNames{5}).relevant_time, satrec.(fNames{6}).relevant_time, satrec.(fNames{7}).relevant_time, satrec.(fNames{8}).relevant_time)
