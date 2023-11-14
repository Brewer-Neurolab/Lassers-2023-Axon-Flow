%% Cumulative well-tunnel/subregion-tunnel scatter plots
% Sam Lassers
% This program uses data from AFR tables to create a scatter plot of
% frequency responses over an entire recording either between a specified
% electrode-tunnel pair, or between an entire subergion of electrodes and
% the entirety of a set of tunnels. 

%% set up 

load("D:\Brewer lab data\HFS\Temporal Analysis\No Stim\dataInfo.mat")
load("D:\Brewer lab data\HFS\No Stim\23-Nov-2021_A\allregion_unit_matched_stim.mat")
% load("D:\Brewer lab data\HFS\Temporal Analysis\40 HFS\dataInfo.mat")
% load("D:\Brewer lab data\HFS\Theta Stim\03-Mar-2022_B\allregion_unit_matched_stim.mat")
% load("D:\Brewer lab data\HFS\Temporal Analysis\5 HFS\dataInfo.mat")
% load("D:\Brewer lab data\HFS\HFS Stim\24-Nov-2021_A\allregion_unit_matched_stim.mat")
allregion_unit_matched=convert_allregion_unit_matched_220222(allregion_unit_matched_stim);

data_folder_addr = "D:\Brewer lab data\HFS\No Stim\xcorr_pseudo_times";
% data_folder_addr = "D:\Brewer lab data\HFS\Theta Stim\xcorr_pseudo_times";
% data_folder_addr = "D:\Brewer lab data\HFS\HFS Stim\xcorr_pseudo_times";
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

source = load('D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output 1s bin AFR\concatenated_source_table.mat');
target = load('D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output 1s bin AFR\concatenated_target_table.mat');
% source = load('D:\Brewer lab data\HFS\Temporal Analysis\40 HFS\AFR Output\concatenated_source_table.mat');
% target = load('D:\Brewer lab data\HFS\Temporal Analysis\40 HFS\AFR Output\concatenated_target_table.mat');
% source = load('D:\Brewer lab data\HFS\Temporal Analysis\5 HFS\AFR Output\concatenated_source_table.mat');
% target = load('D:\Brewer lab data\HFS\Temporal Analysis\5 HFS\AFR Output\concatenated_target_table.mat');

NodeTable = create_node_table_200820;
regList = ["EC","DG","CA3","CA1"];
source_regList_labels = ["EC->tunnels","DG->tunnels","CA3->tunnels","CA1->tunnels"];
target_regList_labels = ["tunnels->EC","tunnels->DG","tunnels->CA3","tunnels->CA1"];

%% well-tunnel

%% Plotting fi=1 axon c6, b3, timebin 113
find_in_AFR_table = @(table, tunnel_name, well_name, fi) find(contains(table.tunnel_names, tunnel_name) & contains(table.well_names,well_name) & table.fi==fi);
ids = find_in_AFR_table(target.AFR_table, 'G5-u1', 'L3', 2);
rowi = ids;
plot_AFR_scatter_times(source.AFR_table,rowi, "tunnel", 5)

%% subregion-tunnel

%%
find_in_AFR_table = @(table, well_reg, tunnel_reg, fi) find(contains(table.well_reg, well_reg) & contains(table.tunnel_reg, tunnel_reg) & (table.fi==fi));
ids = find_in_AFR_table(source.AFR_table, 'EC', 'EC-DG',2);
rowi = ids;
plot_AFR_scatter_times(source.AFR_table,rowi,"source",5)

%%
find_in_AFR_table = @(table, well_reg, tunnel_reg, fi) find(contains(table.well_reg, well_reg) & contains(table.tunnel_reg, tunnel_reg) & (table.fi==fi));
ids = find_in_AFR_table(source.AFR_table, 'CA3', 'CA3-CA1',2);
rowi = ids;
plot_AFR_scatter_times(source.AFR_table,rowi,"source",5)
%% Related functions
function plot_AFR_scatter_times(AFR_table, rowi, well_type, no_points_thresh)
tunnel_AFR=[];
well_AFR=[];

for i=1:length(rowi)
    tunnel_AFR = [tunnel_AFR, AFR_table.tunnel_AFR{rowi}];
    well_AFR =  [well_AFR, AFR_table.well_AFR{rowi}];
end

partial_AFR=AFR_table(rowi,:);

%time_dep_AFR_table=compute_allregion_time_dependent_linear_model(AFR_table,300,5,well_type);

% Plot scatter
figure( 'Position', [100 100 700 600])
z_idx = tunnel_AFR == 0 | well_AFR == 0;
if well_type=="tunnel"
    scatter(tunnel_AFR(~(z_idx)),well_AFR(~(z_idx)),'filled')
    xlabel(char(AFR_table.tunnel_reg(rowi(1)) + " Tunnel Spike Rate (Hz)"));
    ylabel(char(AFR_table.well_reg(rowi(1)) + " Well Spike Rate (Hz)"))
elseif well_type=="source"
    scatter(well_AFR(~(z_idx)),tunnel_AFR(~(z_idx)),'filled')
    ylabel(char(AFR_table.tunnel_reg(rowi(1)) + " Tunnel Spike Rate (Hz)"));
    xlabel(char(AFR_table.well_reg(rowi(1)) + " Well Spike Rate (Hz)"))
end
title({char("Tunnel unit - "+AFR_table.tunnel_names(rowi(1))+" thru "+AFR_table.tunnel_names(rowi(end))),...
    char("Well el - "+AFR_table.well_names(rowi(1)))})

% Squaring the axis
x_end = get(gca,'Xlim'); x_end = x_end(2);
y_end = get(gca,'Ylim'); y_end = y_end(2);
ax_end = max(x_end, y_end);

ylim([-inf ax_end]), xlim([-inf, ax_end])
hold on,line([0 ax_end], [0 ax_end], 'lineStyle','--',...
        'Color','k')
axis square

% calculate new LM

if length(tunnel_AFR(~z_idx)) > no_points_thresh
    % fitting linear regression model
    if well_type == "tunnel"
        mdl = fitlm(tunnel_AFR(~z_idx), well_AFR(~z_idx));
    else 
        mdl = fitlm(well_AFR(~z_idx), tunnel_AFR(~z_idx));
    end
    coeff = mdl.Coefficients.Estimate;                              %[intercept, slope]
    pval = mdl.Coefficients{2,4};
    rsq = mdl.Rsquared.Ordinary;

    % adding to the table
    slope =  coeff(2);
    rsq =  rsq;
    intercept =  coeff(1);
    p_val =  pval;
else 

    slope{rowi} =  NaN;
    rsq{rowi} =  NaN;
    intercept{rowi} =  NaN;
    p_val{rowi} =  NaN;
end

% plotting linear line
b1=[];
b0=[];

%for i=1:length(rowi)
    b1= slope;
    b0 = intercept;
%end

mdl_y = b0 + b1*ax_end;
line([0 ax_end], [b0 +mdl_y], 'color','r')
hold off
end

% function time_depn_AFR_table = compute_allregion_time_dependent_linear_model(AFR_table, bin_lm_s, no_points_thresh, well_type)
% % Computes time dependent linear mdl for AFR_table with time bins = time_s
% time_vec = AFR_table.Properties.CustomProperties.time_vec;
% time_diff = mean(unique(diff(time_vec)));
% total_time_s = time_vec(end);
% time_depn_AFR_table = AFR_table;
% time_depn_AFR_table.slope = cell(height(time_depn_AFR_table),1);
% time_depn_AFR_table.rsq = cell(height(time_depn_AFR_table),1);
% time_depn_AFR_table.intercept= cell(height(time_depn_AFR_table),1);
% time_depn_AFR_table.p_val = cell(height(time_depn_AFR_table),1);
% 
% for t_start = 0.5*time_diff:bin_lm_s:total_time_s - 0.5*time_diff
%     for rowi = 1:height(AFR_table)
%         % removing zeros
%         idx = time_vec >= t_start & time_vec <=t_start+bin_lm_s-time_diff;
%         tunnel_AFR = AFR_table.tunnel_AFR{rowi}(idx);
%         well_AFR = AFR_table.well_AFR{rowi}(idx);
%         z_idx = tunnel_AFR == 0 |  well_AFR == 0;
%         tunnel_AFR(z_idx) = []; well_AFR(z_idx) = [];
%         if length(tunnel_AFR) > no_points_thresh
%             % fitting linear regression model
%             if well_type == "target"
%                 mdl = fitlm(tunnel_AFR, well_AFR);
%             else 
%                 mdl = fitlm(well_AFR, tunnel_AFR);
%             end
%             coeff = mdl.Coefficients.Estimate;                              %[intercept, slope]
%             pval = mdl.Coefficients{2,4};
%             rsq = mdl.Rsquared.Ordinary;
% 
%             % adding to the table
%             time_depn_AFR_table.slope{rowi} = [time_depn_AFR_table.slope{rowi}, coeff(2)];
%             time_depn_AFR_table.rsq{rowi} = [time_depn_AFR_table.rsq{rowi}, rsq];
%             time_depn_AFR_table.intercept{rowi} = [time_depn_AFR_table.intercept{rowi}, coeff(1)];
%             time_depn_AFR_table.p_val{rowi} = [time_depn_AFR_table.p_val{rowi}, pval];
%         else 
%             
%             time_depn_AFR_table.slope{rowi} = [time_depn_AFR_table.slope{rowi}, NaN];
%             time_depn_AFR_table.rsq{rowi} = [time_depn_AFR_table.rsq{rowi}, NaN];
%             time_depn_AFR_table.intercept{rowi} = [time_depn_AFR_table.intercept{rowi}, NaN];
%             time_depn_AFR_table.p_val{rowi} = [time_depn_AFR_table.p_val{rowi}, NaN];
%         end
% 
%     end
%    
% end
% end