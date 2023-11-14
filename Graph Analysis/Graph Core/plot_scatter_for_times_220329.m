%% Scatter plots for AFR 300s Linear model

%% load prerequisites
clear
clc
close all
%Loading prerequisites
load("D:\Brewer lab data\HFS\Temporal Analysis\No Stim\dataInfo.mat")
% load("D:\Brewer lab data\HFS\No Stim\23-Nov-2021_A\allregion_unit_matched_stim.mat")
load("D:\Brewer lab data\HFS\Temporal Analysis\No Stim\full_idx_allregion_unit_matched_stim.mat")

% load("D:\Brewer lab data\HFS\Temporal Analysis\40 HFS\dataInfo.mat")
% load("D:\Brewer lab data\HFS\Theta Stim\03-Mar-2022_B\allregion_unit_matched_stim.mat")
% load("D:\Brewer lab data\HFS\Temporal Analysis\40 HFS\full_idx_allregion_unit_matched_stim.mat")

% load("D:\Brewer lab data\HFS\Temporal Analysis\5 HFS\dataInfo.mat")
% load("D:\Brewer lab data\HFS\HFS Stim\24-Nov-2021_A\allregion_unit_matched_stim.mat")
% load("D:\Brewer lab data\HFS\Temporal Analysis\5 HFS\full_idx_allregion_unit_matched_stim.mat")

allregion_unit_matched=convert_allregion_unit_matched_220222(allregion_unit_matched_stim);

% partial index
% data_folder_addr = "D:\Brewer lab data\HFS\No Stim\xcorr_pseudo_times";
% data_folder_addr = "D:\Brewer lab data\HFS\Theta Stim\xcorr_pseudo_times";
% data_folder_addr = "D:\Brewer lab data\HFS\HFS Stim\xcorr_pseudo_times";

% full spike index
data_folder_addr="D:\Brewer lab data\HFS\No Stim\full_index_pseudo_times";
% data_folder_addr="D:\Brewer lab data\HFS\Theta Stim\full_index_pseudo_times";
% data_folder_addr="D:\Brewer lab data\HFS\HFS Stim\full_index_pseudo_times";

data_folder_dir=dir(data_folder_addr);
data_folder_isdir=[data_folder_dir.isdir];
data_folder_names=string({data_folder_dir(data_folder_isdir).name});
data_folder_names=data_folder_names(3:end);
%data_folder_names=erase(data_folder_names,"_mat_files")';
all_region_order=[];
for i=1:length(data_folder_names)
    all_region_order(i)=find(contains(data_folder_names,dataInfo.meaName(i)));
end
allregion_unit_matched=allregion_unit_matched(all_region_order);

%no stim
source=load('D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output 1s bin AFR\concatenated_source_table.mat');
target=load('D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output 1s bin AFR\concatenated_target_table.mat');

% 40 HFS
% source = load('D:\Brewer lab data\HFS\Temporal Analysis\40 HFS\AFR Output\concatenated_source_table.mat');
% target = load('D:\Brewer lab data\HFS\Temporal Analysis\40 HFS\AFR Output\concatenated_target_table.mat');

% source=load("D:\Brewer lab data\HFS\Temporal Analysis\40 HFS\AFR Output Full Index\concatenated_source_table.mat");
% target=load("D:\Brewer lab data\HFS\Temporal Analysis\40 HFS\AFR Output Full Index\concatenated_target_table.mat");

% 5 HFS
% source = load('D:\Brewer lab data\HFS\Temporal Analysis\5 HFS\AFR Output\concatenated_source_table.mat');
% target = load('D:\Brewer lab data\HFS\Temporal Analysis\5 HFS\AFR Output\concatenated_target_table.mat');

% source=load("D:\Brewer lab data\HFS\Temporal Analysis\5 HFS\AFR Output Full Index\concatenated_source_table.mat");
% target=load("D:\Brewer lab data\HFS\Temporal Analysis\5 HFS\AFR Output Full Index\concatenated_target_table.mat");

NodeTable = create_node_table_200820;
regList = ["EC","DG","CA3","CA1"];
source_regList_labels = ["EC->tunnels","DG->tunnels","CA3->tunnels","CA1->tunnels"];
target_regList_labels = ["tunnels->EC","tunnels->DG","tunnels->CA3","tunnels->CA1"];

%% Source J7-u1, L9 no stim fb
find_in_AFR_table = @(table, tunnel_name, well_name, fi) find(contains(table.tunnel_names, tunnel_name) & contains(table.well_names,well_name) & table.fi==fi)
ids = find_in_AFR_table(source.AFR_table, 'J7-u1', 'L9',2);
rowi = ids;
%bini=300;
plot_AFR_scatter_times_source(source.AFR_table,rowi)

rsq=source.AFR_table.rsq{rowi}
slope=source.AFR_table.slope{rowi}
pval=source.AFR_table.p_val{rowi}

saveas(gcf,"D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output full idx\time_depn_graph_plots\2\LM Scatter\J7u1L9sourcefb.png")
%% source G11-u1 K11 nostim FF

find_in_AFR_table = @(table, tunnel_name, well_name, fi) find(contains(table.tunnel_names, tunnel_name) & contains(table.well_names,well_name) & table.fi==fi)
ids = find_in_AFR_table(source.AFR_table, 'G11-u1', 'K11',2);
rowi = ids;
%bini=300;
plot_AFR_scatter_times_source(source.AFR_table,rowi)

saveas(gcf,"D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output full idx\time_depn_graph_plots\2\LM Scatter\G11u1K11sourceff.png")

%% Source B10, A7
find_in_AFR_table = @(table, tunnel_name, well_name, fi) find(contains(table.tunnel_names, tunnel_name) & contains(table.well_names,well_name) & table.fi==fi)
ids = find_in_AFR_table(source.AFR_table, 'A7-u1', 'B10',2);
rowi = ids;
%bini=300;
plot_AFR_scatter_times_source(source.AFR_table,rowi)
rsq=source.AFR_table.rsq{rowi}
slope=source.AFR_table.slope{rowi}
pval=source.AFR_table.p_val{rowi}
saveas(gcf,"D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output full idx\time_depn_graph_plots\2\LM Scatter\A7u1B10sourceff.png")

%% target nostim
find_in_AFR_table = @(table, tunnel_name, well_name, fi) find(contains(table.tunnel_names, tunnel_name) & contains(table.well_names,well_name) & table.fi==fi)
ids = find_in_AFR_table(target.AFR_table, 'H7-u1', 'J9',2);
rowi = ids;
%bini=300;
plot_AFR_scatter_times_target(target.AFR_table,rowi)

saveas(gcf,"D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output full idx\time_depn_graph_plots\2\LM Scatter\H7u1J9targetff.png")

%% target nostim
find_in_AFR_table = @(table, tunnel_name, well_name, fi) find(contains(table.tunnel_names, tunnel_name) & contains(table.well_names,well_name) & table.fi==fi)
ids = find_in_AFR_table(target.AFR_table, 'H7-u1', 'K9',2);
rowi = ids;
%bini=300;
plot_AFR_scatter_times_target(target.AFR_table,rowi)

saveas(gcf,"D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output full idx\time_depn_graph_plots\2\LM Scatter\H7u1K9targetff.png")

%% target nostim
find_in_AFR_table = @(table, tunnel_name, well_name, fi) find(contains(table.tunnel_names, tunnel_name) & contains(table.well_names,well_name) & table.fi==fi)
ids = find_in_AFR_table(target.AFR_table, 'G3-u1', 'H2',2);
rowi = ids;
%bini=300;
plot_AFR_scatter_times_target(target.AFR_table,rowi)

saveas(gcf,"D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output full idx\time_depn_graph_plots\2\LM Scatter\G3u1H2targetff.png")

%% target nostim
find_in_AFR_table = @(table, tunnel_name, well_name, fi) find(contains(table.tunnel_names, tunnel_name) & contains(table.well_names,well_name) & table.fi==fi)
ids = find_in_AFR_table(target.AFR_table, 'C7-u1', 'C8',2);
rowi = ids;
%bini=300;
plot_AFR_scatter_times_target(target.AFR_table,rowi)

saveas(gcf,"D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output full idx\time_depn_graph_plots\2\LM Scatter\C7u1C8targetff.png")

%% target no stim J6 L3
find_in_AFR_table = @(table, tunnel_name, well_name, fi) find(contains(table.tunnel_names, tunnel_name) & contains(table.well_names,well_name) & table.fi==fi)
ids = find_in_AFR_table(target.AFR_table, 'J6-u1', 'L3',2);
rowi = ids;
%bini=300;
plot_AFR_scatter_times_target(target.AFR_table,rowi)
rsq=target.AFR_table.rsq{rowi}
slope=target.AFR_table.slope{rowi}
pval=target.AFR_table.p_val{rowi}
saveas(gcf,"D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output full idx\time_depn_graph_plots\2\LM Scatter\J6u1L3targetff.png")

%% target no stim J6 L3 100ms AFR
% source=load('D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output full idx\concatenated_source_table.mat');
target=load('D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output full idx\concatenated_target_table.mat');
find_in_AFR_table = @(table, tunnel_name, well_name, fi) find(contains(table.tunnel_names, tunnel_name) & contains(table.well_names,well_name) & table.fi==fi)
ids = find_in_AFR_table(target.AFR_table, 'J6-u1', 'L3',2);
rowi = ids;

for i=1:length(target.AFR_table.slope)
    plot_partial_AFR_scatter_times_target(target.AFR_table, rowi, i)
    saveas(gcf,"D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output full idx\time_depn_graph_plots\2\LM Scatter\J6u1L3targetff100ms\bin"+string(i)+".png")
    
end

%% target no stim J6 L3 100ms 300s AFR
% source=load('D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output full idx\concatenated_source_table.mat');
target=load('D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output 100ms bin AFR\concatenated_target_table.mat');
find_in_AFR_table = @(table, tunnel_name, well_name, fi) find(contains(table.tunnel_names, tunnel_name) & contains(table.well_names,well_name) & table.fi==fi)
ids = find_in_AFR_table(target.AFR_table, 'J6-u1', 'L3',2);
rowi = ids;

plot_AFR_scatter_times_target(target.AFR_table, rowi)
% saveas(gcf,"D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output full idx\time_depn_graph_plots\2\LM Scatter\J6u1L3targetff100ms\bin"+string(i)+".png")
    
%% target no stim J6 H2
find_in_AFR_table = @(table, tunnel_name, well_name, fi) find(contains(table.tunnel_names, tunnel_name) & contains(table.well_names,well_name) & table.fi==fi)
ids = find_in_AFR_table(target.AFR_table, 'J6-u1', 'H2',2);
rowi = ids;
%bini=300;
plot_AFR_scatter_times_target(target.AFR_table,rowi)

saveas(gcf,"D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output full idx\time_depn_graph_plots\2\LM Scatter\J6u1H2targetff.png")

%% target no stim A6 B3

find_in_AFR_table = @(table, tunnel_name, well_name, fi) find(contains(table.tunnel_names, tunnel_name) & contains(table.well_names,well_name) & table.fi==fi)
ids = find_in_AFR_table(target.AFR_table, 'A6-u1', 'B3',2);
rowi = ids;
%bini=300;
plot_AFR_scatter_times_target(target.AFR_table,rowi)

saveas(gcf,"D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output full idx\time_depn_graph_plots\2\LM Scatter\A6u1B3targetff.png")
%% target no stim A6 B4

find_in_AFR_table = @(table, tunnel_name, well_name, fi) find(contains(table.tunnel_names, tunnel_name) & contains(table.well_names,well_name) & table.fi==fi)
ids = find_in_AFR_table(target.AFR_table, 'A6-u1', 'B4',2);
rowi = ids;
%bini=300;
plot_AFR_scatter_times_target(target.AFR_table,rowi)
rsq=target.AFR_table.rsq{rowi}
slope=target.AFR_table.slope{rowi}
pval=target.AFR_table.p_val{rowi}
saveas(gcf,"D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output full idx\time_depn_graph_plots\2\LM Scatter\A6u1B4targetff.png")
%% target no stim A6 D4

find_in_AFR_table = @(table, tunnel_name, well_name, fi) find(contains(table.tunnel_names, tunnel_name) & contains(table.well_names,well_name) & table.fi==fi)
ids = find_in_AFR_table(target.AFR_table, 'A6-u1', 'B4',2);
rowi = ids;
%bini=300;
plot_AFR_scatter_times_target(target.AFR_table,rowi)

saveas(gcf,"D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output full idx\time_depn_graph_plots\2\LM Scatter\A6u1D4targetff.png")
%% target no stim A6 D4

find_in_AFR_table = @(table, tunnel_name, well_name, fi) find(contains(table.tunnel_names, tunnel_name) & contains(table.well_names,well_name) & table.fi==fi)
ids = find_in_AFR_table(target.AFR_table, 'G3-u1', 'L3',2);
rowi = ids;
%bini=300;
plot_AFR_scatter_times_target(target.AFR_table,rowi)

% saveas(gcf,"D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output full idx\time_depn_graph_plots\2\LM Scatter\A6u1D4targetff.png")

%% 40 HFS quick load prereq
close all
%Loading prerequisites

load("D:\Brewer lab data\HFS\Temporal Analysis\40 HFS\dataInfo.mat")
% load("D:\Brewer lab data\HFS\Theta Stim\03-Mar-2022_B\allregion_unit_matched_stim.mat")
load("D:\Brewer lab data\HFS\Temporal Analysis\40 HFS\full_idx_allregion_unit_matched_stim.mat")

allregion_unit_matched=convert_allregion_unit_matched_220222(allregion_unit_matched_stim);

% partial index
% data_folder_addr = "D:\Brewer lab data\HFS\No Stim\xcorr_pseudo_times";
% data_folder_addr = "D:\Brewer lab data\HFS\Theta Stim\xcorr_pseudo_times";
% data_folder_addr = "D:\Brewer lab data\HFS\HFS Stim\xcorr_pseudo_times";

% full spike index
% data_folder_addr="D:\Brewer lab data\HFS\No Stim\full_index_pseudo_times";
data_folder_addr="D:\Brewer lab data\HFS\Theta Stim\full_index_pseudo_times";
% data_folder_addr="D:\Brewer lab data\HFS\HFS Stim\full_index_pseudo_times";

data_folder_dir=dir(data_folder_addr);
data_folder_isdir=[data_folder_dir.isdir];
data_folder_names=string({data_folder_dir(data_folder_isdir).name});
data_folder_names=data_folder_names(3:end);
%data_folder_names=erase(data_folder_names,"_mat_files")';
all_region_order=[];
for i=1:length(data_folder_names)
    all_region_order(i)=find(contains(data_folder_names,dataInfo.meaName(i)));
end
allregion_unit_matched=allregion_unit_matched(all_region_order);

% 40 HFS
source = load('D:\Brewer lab data\HFS\Temporal Analysis\40 HFS\AFR Output 1bin 300s LM\concatenated_source_table.mat');
target = load('D:\Brewer lab data\HFS\Temporal Analysis\40 HFS\AFR Output 1bin 300s LM\concatenated_target_table.mat');

NodeTable = create_node_table_200820;
regList = ["EC","DG","CA3","CA1"];
source_regList_labels = ["EC->tunnels","DG->tunnels","CA3->tunnels","CA1->tunnels"];
target_regList_labels = ["tunnels->EC","tunnels->DG","tunnels->CA3","tunnels->CA1"];

%% soruce 40 HFS
find_in_AFR_table = @(table, tunnel_name, well_name, fi) find(contains(table.tunnel_names, tunnel_name) & contains(table.well_names,well_name) & table.fi==fi)
ids = find_in_AFR_table(source.AFR_table, 'D7-u1', 'C10',2);
rowi = ids;
%bini=300;
plot_AFR_scatter_times_source(source.AFR_table,rowi)

% saveas(gcf,"D:\Brewer lab data\HFS\Temporal Analysis\40 HFS\AFR Output full idx\time_depn_graph_plots\2\LM Scatter\G11u1K11sourceff.png")

%% target 40 HFS G5 L3
find_in_AFR_table = @(table, tunnel_name, well_name, fi) find(contains(table.tunnel_names, tunnel_name) & contains(table.well_names,well_name) & table.fi==fi)
ids = find_in_AFR_table(target.AFR_table, 'G3-u1', 'K4',2);
rowi = ids;
%bini=300;
plot_AFR_scatter_times_target(target.AFR_table,rowi)

% saveas(gcf,"D:\Brewer lab data\HFS\Temporal Analysis\40 HFS\AFR Output full idx\time_depn_graph_plots\2\LM Scatter\G11u1K11sourceff.png")


%% Related functions
function plot_AFR_scatter_times_source(AFR_table, rowi)
timi = 1:3000;
tunnel_AFR = AFR_table.tunnel_AFR{rowi}(timi) ;
well_AFR =  AFR_table.well_AFR{rowi}(timi);
% Plot scatter
figure( 'Position', [100 100 700 600])
z_idx = tunnel_AFR == 0 | well_AFR == 0;
%scatter(tunnel_AFR(~(z_idx)),well_AFR(~(z_idx)),'filled')
scatter(well_AFR(~(z_idx)),tunnel_AFR(~(z_idx)),'filled')
% xlabel(char(AFR_table.tunnel_reg(rowi) + " Tunnel Spike Rate (Hz)"));
% ylabel(char(AFR_table.well_reg(rowi) + " Well Spike Rate (Hz)"))

ylabel(char(AFR_table.tunnel_reg(rowi) + " Tunnel Spike Rate (Hz)"));
xlabel(char(AFR_table.well_reg(rowi) + " Well Spike Rate (Hz)"))

if AFR_table.if_ff(rowi)
title({char("Tunnel unit - "+AFR_table.tunnel_names(rowi)+ ' (ff)'),...
    char("Well el - "+AFR_table.well_names(rowi))})
else
    title({char("Tunnel unit - "+AFR_table.tunnel_names(rowi)+ ' (fb)'),...
    char("Well el - "+AFR_table.well_names(rowi))})
end

% Squaring the axis
x_end = get(gca,'Xlim'); x_end = x_end(2);
y_end = get(gca,'Ylim'); y_end = y_end(2);
ax_end = max(x_end, y_end);

ylim([-inf ax_end]), xlim([-inf, ax_end])
hold on,line([0 ax_end], [0 ax_end], 'lineStyle','--',...
        'Color','k')
axis square

% plotting linear line
b1= AFR_table.slope{rowi}(1);
b0 = AFR_table.intercept{rowi}(1);
mdl_y = b0 + b1*ax_end;
line([0 ax_end], [b0 +mdl_y], 'color','r')
hold off
end

function plot_AFR_scatter_times_target(AFR_table, rowi)
timi = 1:3000;
tunnel_AFR = AFR_table.tunnel_AFR{rowi}(timi) ;
well_AFR =  AFR_table.well_AFR{rowi}(timi);
% Plot scatter
figure( 'Position', [100 100 700 600])
z_idx = tunnel_AFR == 0 | well_AFR == 0;
scatter(tunnel_AFR(~(z_idx)),well_AFR(~(z_idx)),'filled')
% scatter(well_AFR(~(z_idx)),tunnel_AFR(~(z_idx)),'filled')
xlabel(char(AFR_table.tunnel_reg(rowi) + " Tunnel Spike Rate (Hz)"));
ylabel(char(AFR_table.well_reg(rowi) + " Well Spike Rate (Hz)"))

% ylabel(char(AFR_table.tunnel_reg(rowi) + " Tunnel Spike Rate (Hz)"));
% xlabel(char(AFR_table.well_reg(rowi) + " Well Spike Rate (Hz)"))

if AFR_table.if_ff(rowi)
title({char("Tunnel unit - "+AFR_table.tunnel_names(rowi)+ ' (ff)'),...
    char("Well el - "+AFR_table.well_names(rowi))})
else
    title({char("Tunnel unit - "+AFR_table.tunnel_names(rowi)+ ' (fb)'),...
    char("Well el - "+AFR_table.well_names(rowi))})
end

% Squaring the axis
x_end = get(gca,'Xlim'); x_end = x_end(2);
y_end = get(gca,'Ylim'); y_end = y_end(2);
ax_end = max(x_end, y_end);

ylim([-inf ax_end]), xlim([-inf, ax_end])
hold on,line([0 ax_end], [0 ax_end], 'lineStyle','--',...
        'Color','k')
axis square

% plotting linear line
b1= AFR_table.slope{rowi}(1);
b0 = AFR_table.intercept{rowi}(1);
mdl_y = b0 + b1*ax_end;
line([0 ax_end], [b0 +mdl_y], 'color','r')
hold off
end

function plot_partial_AFR_scatter_times_target(AFR_table, rowi, bini)
timi = ((bini-1)*10+1):((bini-1)*10+10);
tunnel_AFR = AFR_table.tunnel_AFR{rowi}(timi) ;
well_AFR =  AFR_table.well_AFR{rowi}(timi);
% Plot scatter
figure( 'Position', [100 100 700 600])
z_idx = tunnel_AFR == 0 | well_AFR == 0;
scatter(tunnel_AFR(~(z_idx)),well_AFR(~(z_idx)),'filled')
xlabel(char(AFR_table.tunnel_reg(rowi) + " Tunnel Spike Rate (Hz)"));
ylabel(char(AFR_table.well_reg(rowi) + " Well Spike Rate (Hz)"))
title({char("Tunnel unit - "+AFR_table.tunnel_names(rowi)),...
    char("Well el - "+AFR_table.well_names(rowi)),...
    char("Time idx - "+ bini)})

% Squaring the axis
x_end = get(gca,'Xlim'); x_end = x_end(2);
y_end = get(gca,'Ylim'); y_end = y_end(2);
ax_end = max(x_end, y_end);

ylim([-inf ax_end]), xlim([-inf, ax_end])
hold on,line([0 ax_end], [0 ax_end], 'lineStyle','--',...
        'Color','k')
axis square

% plotting linear line
b1= AFR_table.slope{rowi}(bini);
b0 = AFR_table.intercept{rowi}(bini);
mdl_y = b0 + b1*ax_end;
line([0 ax_end], [b0 +mdl_y], 'color','r')
hold off
end

function plot_partial_AFR_scatter_times_source(AFR_table, rowi, bini)
timi = ((bini-1)*10+1):((bini-1)*10+10);
tunnel_AFR = AFR_table.tunnel_AFR{rowi}(timi) ;
well_AFR =  AFR_table.well_AFR{rowi}(timi);
% Plot scatter
figure( 'Position', [100 100 700 600])
z_idx = tunnel_AFR == 0 | well_AFR == 0;
scatter(well_AFR(~(z_idx)),tunnel_AFR(~(z_idx)),'filled')
ylabel(char(AFR_table.tunnel_reg(rowi) + " Tunnel Spike Rate (Hz)"));
xlabel(char(AFR_table.well_reg(rowi) + " Well Spike Rate (Hz)"))
title({char("Tunnel unit - "+AFR_table.tunnel_names(rowi)),...
    char("Well el - "+AFR_table.well_names(rowi)),...
    char("Time idx - "+ bini)})

% Squaring the axis
x_end = get(gca,'Xlim'); x_end = x_end(2);
y_end = get(gca,'Ylim'); y_end = y_end(2);
ax_end = max(x_end, y_end);

ylim([-inf ax_end]), xlim([-inf, ax_end])
hold on,line([0 ax_end], [0 ax_end], 'lineStyle','--',...
        'Color','k')
axis square

% plotting linear line
b1= AFR_table.slope{rowi}(bini);
b0 = AFR_table.intercept{rowi}(bini);
mdl_y = b0 + b1*ax_end;
line([0 ax_end], [b0 +mdl_y], 'color','r')
hold off
end