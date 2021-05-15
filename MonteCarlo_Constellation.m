clc;clear all;


angdist = 24; %maximum angular distance between ground tracks at equator
numplane = 12;
plane_sep = angdist./(numplane);
plane_val = plane_sep:plane_sep:angdist;
plane_val = [0,plane_val];
% plane_val = plane_val(1:end-1);

numSamp = 100000;
max_trans_per_plane = zeros(numSamp,1);
for j = 1:numSamp
    
    x = rand([20,1])*angdist; % generate random values from 1 to max angular distance

    plane_idx = zeros(length(x),1);
    for i = 1:length(x)
        %for each randomly selected target location, check to see the closest plane 
        c = abs(x(i)-plane_val);
        [~, plane_idx(i)] = min(c);
        if plane_idx(i) == 1
            plane_idx(i) = length(plane_val);
        end
    end
    
    %Count the number of times each plane is used
    for k = 1:length(plane_val)
        npp(k) = sum(plane_idx(:) == k);
    end
    
        

    max_trans_per_plane(j) = max(npp);
    
end


% plane_idx = plane_idx; %don't include the start/end twice
histogram(max_trans_per_plane,'Normalization','probability')
hold on