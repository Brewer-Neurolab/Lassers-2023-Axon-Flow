function [area,conduction_time]= Spike_Hist_Area_210406(hist_vec,z,current_elec_pair)
% Sam Lassers Brewer Lab 1/8/21
% Starting out, this function has an arbitrary bounding method
% Make minimum bin size (get distribution of spikes) 
%% TODO
% Create a running list of last created png images created. 
% Compare area to the rest of the histogram. If it does not capture at
% least 80% of the area, count as indeterminate (-1). Check with the lab
% appropriate area to capture.

%%
% Flip dataset if feedback
if z==0
    hist_vec=flip(hist_vec);
end

% Set first four bins=0 because impossible conduction times
% hist_vec(1:4)=0; %Handled elsewhere

% Find location of max datapoint 
max_pk_loc=find(hist_vec==max(hist_vec));

% User defined width of summation
hist_width_to_sum=2;

% Find bounds to histogram indecies to be summed

% Left most
if max_pk_loc-hist_width_to_sum>0
    left_bound=max_pk_loc-hist_width_to_sum;
else
    left_bound=1;
end

% Right most
size_vec=size(hist_vec);
if max_pk_loc+hist_width_to_sum<size_vec(2)
    right_bound=max_pk_loc+hist_width_to_sum;
else
    right_bound=size_vec(2);
end

% Get area
area=sum(hist_vec(left_bound:right_bound));

% Bar graph
figure(2)
%For 40 khz
%x=(0.05:0.05:1);
%For 25khz
x=(0:1/12:1);
conduction_time=x(hist_vec==max(hist_vec));
b=bar(x,hist_vec);
%Sets summed bars to red
b.FaceColor='flat';
for i=left_bound:right_bound
    b.CData(i,:)=[1 0 0];
end


if z==0
    xlabel('Conduction time (ms) Feedback');
    name_graph=strcat(string(current_elec_pair),char('_Feedback'),char('_area a= '),string(area));
elseif z==1
    xlabel('Conduction time (ms) Feedforward');
    name_graph=strcat(string(current_elec_pair),char('_Feedforward'),char('_area a= '),string(area));
else
    error('Invalid Direction')
end
ylabel('Spike Counts')
hold on
title(name_graph)
saveas(gcf,strcat(name_graph,'.png'))

end