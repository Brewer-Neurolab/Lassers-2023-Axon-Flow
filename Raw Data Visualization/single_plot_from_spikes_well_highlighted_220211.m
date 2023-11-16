%% Wave Clus Output Plotting
% load both the raw data and times files manually
%this is for wells
%clc
%close all
%% loading
current_program=dbstack();
current_program=current_program.name;

% load spike dyn from folder above current in dir
spike_dyn_path=fullfile('..\','well_spike_dynamics_table*.mat');
dir_spike_dyn=dir(spike_dyn_path);
spike_dyn=dir_spike_dyn.name;
spike_dyn=fullfile("..\",spike_dyn);
load(spike_dyn)

%load matching table from one dir up
load('..\matching_table_wells_CW.mat')

%get current system FID
my_folder=string(pwd);
only_folder=strsplit(my_folder,{'\'});
only_folder=only_folder(end);
up_dir=dir('..\');
up_folders=string({up_dir([up_dir.isdir]).name});
up_folders=up_folders(3:end);
fid=find(contains(up_folders,only_folder));

%get last loaded electrode
historypath = com.mathworks.mlservices.MLCommandHistoryServices.getSessionHistory;
lastexpression = string(historypath(end));
if any(contains(lastexpression,current_program))
    lastexpression=string(historypath(end-1));
end

splt_last=strsplit(lastexpression,{'_','.','load','(',')',char(39)});
vb = cellfun(@isempty, splt_last);
splt_last=splt_last(~vb);
last_electrode=splt_last(1);

burst_bounds=well_spike_dynamics_table.BurstBounds(well_spike_dynamics_table.fi==fid & well_spike_dynamics_table.channel_name==last_electrode);

if isempty(burst_bounds)
    error("Burst_bounds variable empty. Please load spikes and raw data manually.")
end

burst_bounds=burst_bounds{1}./25000;

%% plotting
format long
figure('units','normalized','OuterPosition',[0 0 1 1])
hold on
data=highpass(data,300,25000);
x=[1:length(data)];
x=x./25000;
plot(x,data,'Color','k')

%cluster range inputted by user to for single/multi cluster display

%colors=distinguishable_colors(num_clus,{'w','b'});
xlim([0,x(end)])
y_low=min(data);
y_low=y_low(1)-10;
y_hi=max(data);
y_hi=y_hi(1)+10;
ylim([y_low,y_hi])
ylabel("uV")
xlabel("seconds")
yline(-threshold,'Color','r')
%title(strcat(channels_to_plot(k),' FID',string(j)))
pts=unique((index)./1e3);
[~,ia,~]=intersect(round(x,8),round(pts,8));
plot(pts,data(ia),'.','Color',[0 0 1],'MarkerSize',10)
ax=gca;
ax.FontSize=20;
%yline(threshold)
title(strcat(only_folder,'-',last_electrode))
% patch_y=repmat([y_low,y_hi],size(burst_bounds,1),1);
%xticks(0:0.1:300)

%target_bd=0.124; %CA1 no stim
%target_bd=0.103; %CA1 5hfs
target_bd=0.087; %CA1 40hfs

bd_bounds=[target_bd-(0.1*target_bd),target_bd+(0.1*target_bd)];

for i=1:size(burst_bounds,1)
    f=[1 2 3 4];
    v=[burst_bounds(i,1) y_low;burst_bounds(i,2) y_low;burst_bounds(i,2) y_hi;burst_bounds(i,1) y_hi];
    %v=[burst_bounds(i,1) -Inf;burst_bounds(i,2) -Inf;burst_bounds(i,2) Inf;burst_bounds(i,1) Inf];
    bd=burst_bounds(i,2)-burst_bounds(i,1);
    if bd>=bd_bounds(1) && bd<=bd_bounds(2)
        patch('Faces',f,'Vertices',v,'FaceColor','red','FaceAlpha',0.3)
    else
        patch('Faces',f,'Vertices',v,'FaceColor','yellow','FaceAlpha',0.3)
    end
    
end
hold off

%% Cycle through average spike rate windows

%loading table
load(strcat('D:\Brewer lab data\HFS\No Stim\Wells_5SD_500maxSD\ECDGCA3CA1 19908 150729 150823 d25 5minspont0001.h5\freq_bar','\',last_electrode,'_table.mat'))

for i=1:300/10%length(t.nostim{5})+5
    
    xlims=[t.nostim{5}(i)-5,t.nostim{5}(i)+5];
%     xlims=[t.("5hfs"){5}(i)-5,t.("5hfs"){5}(i)+5];
%     xlims=[t.("40hfs"){5}(i)-5,t.("40hfs"){5}(i)+5];
    
    xlim(xlims);
    
%     xlim([t.nostim{5}(i)-5,t.nostim{5}(i)+5])
    %xlim([t.("5hfs"){5}(i)-5,t.("5hfs"){5}(i)+5])
%     xlim([t.("40hfs"){5}(i)-5,t.("40hfs"){5}(i)+5])
    
    xlims=xlims*sr;
    data_segment=data([xlims(1):xlims])
    ylim('auto')
    % only 10s windows
    %xlim([10*(i-1),10*i])
    
    disp("Press enter for next window.")
    pause
end