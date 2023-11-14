%% Plots a scatter plot from only the points that go into creating the weights

%% Load Prerequisite data
clear 
clc
% no stimulation
load("D:\Brewer lab data\HFS\Temporal Analysis\No Stim\dataInfo.mat")
load("D:\Brewer lab data\HFS\Temporal Analysis\No Stim\full_idx_allregion_unit_matched_stim.mat")
data_folder_addr="D:\Brewer lab data\HFS\No Stim\full_index_pseudo_times";
source_1=load('D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output full idx\concatenated_source_table.mat');
target_1=load('D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output full idx\concatenated_target_table.mat');
source_2=load('D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output 100ms bin AFR\concatenated_source_table.mat');
target_2=load('D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output 100ms bin AFR\concatenated_target_table.mat');


allregion_unit_matched=convert_allregion_unit_matched_220222(allregion_unit_matched_stim);

data_folder_dir=dir(data_folder_addr);
data_folder_isdir=[data_folder_dir.isdir];
data_folder_names=string({data_folder_dir(data_folder_isdir).name});
data_folder_names=data_folder_names(3:end);
all_region_order=[];
for i=1:length(data_folder_names)
    all_region_order(i)=find(contains(data_folder_names,dataInfo.meaName(i)));
end
allregion_unit_matched=allregion_unit_matched(all_region_order);
NodeTable = create_node_table_200820;
regList = ["EC","DG","CA3","CA1"];
source_regList_labels = ["EC->tunnels","DG->tunnels","CA3->tunnels","CA1->tunnels"];
target_regList_labels = ["tunnels->EC","tunnels->DG","tunnels->CA3","tunnels->CA1"];

%% create new afr table
%source
AFR_table=select_edge_points(source_1.AFR_table,source_2.AFR_table,0.1,0.2);
AFR_table=compute_allregion_time_dependent_linear_model(AFR_table,300,20,"source");
source.AFR_table=AFR_table;

%target
AFR_table=select_edge_points(target_1.AFR_table,target_2.AFR_table,0.1,0.2);
AFR_table=compute_allregion_time_dependent_linear_model(AFR_table,300,20,"target");
target.AFR_table=AFR_table;

%% No Stimulation plots
%% target no stim J6 L3
find_in_AFR_table = @(table, tunnel_name, well_name, fi) find(contains(table.tunnel_names, tunnel_name) & contains(table.well_names,well_name) & table.fi==fi)
ids = find_in_AFR_table(target.AFR_table, 'J6-u1', 'L3',2);
rowi = ids;
%bini=300;
plot_AFR_scatter_times_target(target.AFR_table,rowi)
rsq=target.AFR_table.rsq{rowi}
slope=target.AFR_table.slope{rowi}
pval=target.AFR_table.p_val{rowi}
% saveas(gcf,"D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output full idx\time_depn_graph_plots\2\LM Scatter\J6u1L3targetff.png")

%% target no stim H7 L10
find_in_AFR_table = @(table, tunnel_name, well_name, fi) find(contains(table.tunnel_names, tunnel_name) & contains(table.well_names,well_name) & table.fi==fi)
ids = find_in_AFR_table(target.AFR_table, 'H7-u1', 'L10',2);
rowi = ids;
%bini=300;
plot_AFR_scatter_times_target(target.AFR_table,rowi)
rsq=target.AFR_table.rsq{rowi}
slope=target.AFR_table.slope{rowi}
pval=target.AFR_table.p_val{rowi}

%% Related Functions
function new_afr=select_edge_points(AFR1, AFR2, slope_thresh, rsq_thresh)
new_afr=AFR2;
% afr 1 should have 100ms bins
% afr 1 and 2 should have the same height
new_afr=removevars(new_afr,{'slope','rsq','intercept','p_val'});

for rows=1:height(AFR1)
    find_bins=find(AFR1.slope{rows}>slope_thresh & AFR1.rsq{rows}>rsq_thresh);
    include_afr=zeros(1,length(AFR1.tunnel_AFR{rows}));

    for i=1:length(find_bins)
        range=(find_bins(i)*10)-9;
        include_afr(range:find_bins(i)*10)=1;
    end
    include_afr=logical(include_afr);
    new_afr.tunnel_AFR{rows}(~include_afr)=0;
    
end
end

function time_depn_AFR_table = compute_allregion_time_dependent_linear_model(AFR_table, bin_lm_s, no_points_thresh, well_type)
% Computes time dependent linear mdl for AFR_table with time bins = time_s
time_vec = AFR_table.Properties.CustomProperties.time_vec;
time_diff = mean(unique(diff(time_vec)));
total_time_s = time_vec(end);
time_depn_AFR_table = AFR_table;
time_depn_AFR_table.slope = cell(height(time_depn_AFR_table),1);
time_depn_AFR_table.rsq = cell(height(time_depn_AFR_table),1);
time_depn_AFR_table.intercept= cell(height(time_depn_AFR_table),1);
time_depn_AFR_table.p_val = cell(height(time_depn_AFR_table),1);

for t_start = 0.5*time_diff:bin_lm_s:total_time_s - 0.5*time_diff
    for rowi = 1:height(AFR_table)
        % removing zeros
        idx = time_vec >= t_start & time_vec <=t_start+bin_lm_s-time_diff;
        tunnel_AFR = AFR_table.tunnel_AFR{rowi}(idx);
        well_AFR = AFR_table.well_AFR{rowi}(idx);
        z_idx = tunnel_AFR == 0 |  well_AFR == 0;
        tunnel_AFR(z_idx) = []; well_AFR(z_idx) = [];
        if length(tunnel_AFR) > no_points_thresh
            % fitting linear regression model
            if well_type == "target"
                mdl = fitlm(tunnel_AFR, well_AFR);
            else 
                mdl = fitlm(well_AFR, tunnel_AFR);
            end
            coeff = mdl.Coefficients.Estimate;                              %[intercept, slope]
            pval = mdl.Coefficients{2,4};
            rsq = mdl.Rsquared.Ordinary;

            % adding to the table
            time_depn_AFR_table.slope{rowi} = [time_depn_AFR_table.slope{rowi}, coeff(2)];
            time_depn_AFR_table.rsq{rowi} = [time_depn_AFR_table.rsq{rowi}, rsq];
            time_depn_AFR_table.intercept{rowi} = [time_depn_AFR_table.intercept{rowi}, coeff(1)];
            time_depn_AFR_table.p_val{rowi} = [time_depn_AFR_table.p_val{rowi}, pval];
        else 
            
            time_depn_AFR_table.slope{rowi} = [time_depn_AFR_table.slope{rowi}, NaN];
            time_depn_AFR_table.rsq{rowi} = [time_depn_AFR_table.rsq{rowi}, NaN];
            time_depn_AFR_table.intercept{rowi} = [time_depn_AFR_table.intercept{rowi}, NaN];
            time_depn_AFR_table.p_val{rowi} = [time_depn_AFR_table.p_val{rowi}, NaN];
        end

    end
   
end
end

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



