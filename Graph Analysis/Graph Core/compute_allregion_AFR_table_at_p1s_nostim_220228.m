clear
clc
%outfolder = "D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output";
% outfolder = "D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output 1s bin AFR";
outfolder = "D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output full idx";
% outfolder = "D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output 100ms bin AFR";
load("D:\Brewer lab data\HFS\Temporal Analysis\No Stim\dataInfo.mat")
names=dataInfo.meaName;
% MFR_threshold = 0.0333;
MFR_threshold = 0.2;
bin_AFR_s = 0.1; %original 0.1
bin_lm_s=1;  %original 1
no_points_thresh = 5;

%for i=1:length(names)
    
%     if ~exist(outfolder+'\'+names(i),'dir')
%         mkdir(outfolder+'\'+names(i))
%     end
    
    for tunnel_type = ["FF","FB"]
        for well_type = ["source","target"]
            AFR_table = compute_allregion_AFR_table(tunnel_type, well_type, MFR_threshold, bin_AFR_s, bin_lm_s, no_points_thresh, i);
            if ~exist(outfolder,'dir'), mkdir(outfolder); end
            save(outfolder+'\'+tunnel_type+"_"+well_type+"_AFR_table.mat", 'AFR_table')
        end
    end

    rearrange_data_220228(outfolder+'\')
%     rearrange_data_201130

%end

%% Funtions involved 
function AFR_table_w_lmfit = compute_allregion_AFR_table(tunnel_type, well_type, MFR_threshold, bin_AFR_s, bin_lm_s, no_points_thresh, fi)
%% Computes Avergae Firing Rate (AFR) for tunnels and wells (source/ target)
% and computes linear model parameters between (slope, p-value, r-squared) 

% I/P parameters -
% tunnel_type - FF/FB tunnels (string)
% well_type - source/target wells(string)
% MFR_threshold - Mean Firing rate (No. of spikes/ Total time) threshold of 
% number of spikes per channel above which the channel with be used
% (double)
% bin_lm_s - the time resolution of the time depn AFT computation (double)
% no_points_thresh - Minimum AFR values required for linear modelling 

% O/P - AFR_table - Table containing AFR vecs for all tunnel-well pairs


%% Loading prerequisites
load('D:\Brewer lab data\HFS\Temporal Analysis\No Stim\tunnelChannelList.mat')                                  %List of tunnel channels whether near target or source
load('D:\Brewer lab data\HFS\Temporal Analysis\No Stim\dataInfo.mat')                                           %Names of MEA's and their orientations
%load("D:\Brewer lab data\HFS\No Stim\23-Nov-2021_A\allregion_unit_matched_stim.mat")                 %allregion structure with showing matching of sorted unit between tunnel channel sets and their directions
load('D:\Brewer lab data\HFS\Temporal Analysis\No Stim\full_idx_allregion_unit_matched_stim.mat')

%convert allregion_unit_matching from automation data structures- SBL
%22/23/22
allregion_unit_matched=convert_allregion_unit_matched_220222(allregion_unit_matched_stim);

%all region unit matched stim comes in system name order, has to be
%reordered based on data info
% data_folder_addr = "D:\Brewer lab data\HFS\No Stim\xcorr_pseudo_times";

data_folder_addr = "D:\Brewer lab data\HFS\No Stim\full_index_pseudo_times";

data_folder_dir=dir(data_folder_addr);
data_folder_isdir=[data_folder_dir.isdir];
data_folder_names=string({data_folder_dir(data_folder_isdir).name});
data_folder_names=data_folder_names(3:end);
data_folder_names=erase(data_folder_names,"_mat_files");
all_region_order=[];
for i=1:length(data_folder_names)
    all_region_order(i)=find(contains(data_folder_names,dataInfo.meaName(i)));
end
allregion_unit_matched=allregion_unit_matched(all_region_order);
%% Initializing parameters

well_regList = ["EC","DG","CA3","CA1"];
ff_tunnel_regList = ["EC-DG","DG-CA3","CA3-CA1","CA1-EC"];
fb_tunnel_regList = ["DG-EC", "CA3-DG", "CA1-CA3", "EC-CA1"];

fs=25e3;

fi_vec=[]; w_reg={}; t_reg={};
tunnel_names_array = {}; well_names_array = {}; 
well_AFR = {}; tunnel_AFR = {}; 
rowi=1;

if strcmpi(tunnel_type, "FF") && strcmpi(well_type, "target")
    well_regi = [2,3,4,1];
elseif strcmpi(tunnel_type, "FF") && strcmpi(well_type, "source")
    well_regi = 1:4;
elseif strcmpi(tunnel_type, "FB") && strcmpi(well_type, "target")
    well_regi = 1:4;
elseif strcmpi(tunnel_type, "FB") && strcmpi(well_type, "source")
    well_regi = [2,3,4,1];
end

for fi = 1:9
    tunnel_unit_mea_ff = load(data_folder_addr+'\'+dataInfo.meaName{fi}+"_mat_files/mea_tunnels_ff");
    tunnel_unit_mea_fb = load(data_folder_addr+'\'+dataInfo.meaName{fi}+"_mat_files/mea_tunnels_fb");
    %tunnel_unit_mea=[];
    %% obtaining tunnel-channel indices
    ori = dataInfo.orientation{fi};
    if strcmpi(tunnel_type,"FF") && strcmpi(well_type,"target")
        if strcmpi(ori,"cw" )
            tunnel_chan_idx_cell = tunnelChannelList.near_target.cw;
        else
            tunnel_chan_idx_cell = tunnelChannelList.near_target.ccw;
        end
    elseif strcmpi(tunnel_type,"FF") && strcmpi(well_type,"source")
        if strcmpi(ori,"cw" )
            tunnel_chan_idx_cell = tunnelChannelList.near_source.cw;
        else
            tunnel_chan_idx_cell = tunnelChannelList.near_source.ccw;
        end
    elseif strcmpi(tunnel_type,"FB") && strcmpi(well_type,"target")
        if strcmpi(ori,"cw" )
            tunnel_chan_idx_cell = tunnelChannelList.near_source.cw;
        else
            tunnel_chan_idx_cell = tunnelChannelList.near_source.ccw;
        end
    elseif strcmpi(tunnel_type,"FB") && strcmpi(well_type,"source")
        if strcmpi(ori,"cw" )
            tunnel_chan_idx_cell = tunnelChannelList.near_target.cw;
        else
            tunnel_chan_idx_cell = tunnelChannelList.near_target.ccw;
        end
    end
    
    for regi = 1:4
        well_mea = load(data_folder_addr+'\'+dataInfo.meaName{fi}+"_mat_files/mea_"+well_regList(well_regi(regi)));
        
        %% collecting all the well channels
        %Computing MFR and select MFR >= 0.2
        MFR_vec = [];
        for chani = 1:length(well_mea.mea.spike_data)
            MFR_vec(chani) = length(find(well_mea.mea.spike_data{chani}))./(5*60);
        end
        well_count_table = table(well_mea.mea.channel_names, MFR_vec',...
            'VariableNames', {'Channel_name','MFR'});
        well_channel_list = string(well_count_table.Channel_name(well_count_table.MFR >= MFR_threshold))';
        
        if isempty(well_channel_list)
            disp("fi= "+fi+" regi="+regi+" skipped for AFR computation")
            continue
        end
        well_channel_list_idx = find_in_channel_list(well_mea.mea.channel_names, well_channel_list);
        no_wells = length(well_channel_list_idx);
        
        % Running loop to gather all channels
        well_spike_mat = [];
        coli=1;
        for chani=well_channel_list_idx
            well_spike_mat(:,coli) = well_mea.mea.spike_data{chani};
            coli=coli+1;
        end
        
        %% collecting all the tunnel units
        % selecting tunnel channel indices with appropriate dir
        tunnel_unit_spike_mat = []; tunnel_name_list = {};
        for chani=1:length(tunnel_chan_idx_cell{regi})
            tunnel_unit_dir = allregion_unit_matched{fi}{regi}{chani,3};
            tunnel_unit_idx_pair = allregion_unit_matched{fi}{regi}{chani,1};
            tunnel_unit_order = allregion_unit_matched{fi}{regi}{chani,2};
            tunnel_channel_idx = tunnel_chan_idx_cell{regi}(chani);
            tunnel_channel_name = tunnel_unit_mea_ff.mea.channel_names{tunnel_channel_idx}; %have the same channel names in same order, should not matter
            
            % Finding name of the tunnel channel
            if strfind(tunnel_unit_order,tunnel_channel_name) == 1
                u_pos = 1;
            elseif strfind(tunnel_unit_order,tunnel_channel_name) == 4 || ...
                strfind(tunnel_unit_order,tunnel_channel_name) == 5 
                u_pos = 2;
            else
                error('Incompatible channel names in unit match mat ');
            end
                
            if tunnel_type == "FF" || tunnel_type == "ff"
                % selecting FF units
                for ui=1:length(tunnel_unit_dir)
                    if tunnel_unit_dir(ui)
                        u_no = tunnel_unit_idx_pair(ui, u_pos);
                        tunnel_unit_spike_mat = [tunnel_unit_spike_mat, ...
                            tunnel_unit_mea_ff.mea.sorted_spike_data{tunnel_channel_idx}{u_no}];
                        tunnel_name_list = [tunnel_name_list, tunnel_channel_name+"-u"+u_no];
                    end
                end
            elseif tunnel_type == "FB" || tunnel_type == "fb"
                % selecting FB units
                for ui=1:length(tunnel_unit_dir)
                    if ~tunnel_unit_dir(ui) 
                        u_no = tunnel_unit_idx_pair(ui, u_pos);
                        tunnel_unit_spike_mat = [tunnel_unit_spike_mat, ...
                            tunnel_unit_mea_fb.mea.sorted_spike_data{tunnel_channel_idx}{u_no}];
                        tunnel_name_list = [tunnel_name_list, tunnel_channel_name+"-u"+u_no];
                    end
                end
            end
        end
        
        % skipping tunnels no appropriate data
        if isempty(tunnel_unit_spike_mat)
            disp("fi= "+fi+" regi="+regi+" skipped for AFR computation")
            continue
        end
        %% Analysis - Computing AFR for every channel pair
        for w_chani = 1:size(well_spike_mat,2)
            for t_chani = 1:size(tunnel_unit_spike_mat,2)
                well_AFR{rowi} = compute_AFR_200722(well_spike_mat(:,w_chani), bin_AFR_s, fs);
                tunnel_AFR{rowi} = compute_AFR_200722(tunnel_unit_spike_mat(:,t_chani), bin_AFR_s, fs);
                %% storing in appropriate data-structures
                fi_vec(rowi)=fi; w_reg{rowi}=well_regList(well_regi(regi));
                if strcmpi(tunnel_type,"FF")
                    t_reg{rowi}=ff_tunnel_regList(regi);
                else
                    t_reg{rowi}=fb_tunnel_regList(regi);
                end
                tunnel_names_array{rowi}=tunnel_name_list(t_chani);
                well_names_array{rowi}=well_channel_list(w_chani);
                rowi=rowi+1; 
            end
        end
        disp("fi= "+fi+" regi="+regi+" processed for AFR computation")
    end
end
AFR_table = table(fi_vec', string(t_reg'), string(tunnel_names_array'),...
    tunnel_AFR', string(w_reg'), string(well_names_array'), well_AFR',...
    'VariableNames', {'fi', 'tunnel_reg','tunnel_names', 'tunnel_AFR', ...
    'well_reg','well_names','well_AFR'});

%% Storing custom properties of 
AFR_table = addprop(AFR_table, {'tunnel_type','well_type','MFR_threshold','time_vec','if_ff',},...
    {'table','table','table','table','table'});
AFR_table.Properties.CustomProperties.tunnel_type = char(tunnel_type);
AFR_table.Properties.CustomProperties.well_type = char(well_type);
AFR_table.Properties.CustomProperties.MFR_threshold = MFR_threshold;
AFR_table.Properties.CustomProperties.time_vec =  (bin_AFR_s:bin_AFR_s:300) - bin_AFR_s/2;

%% Compute linear model parameters
AFR_table_w_lmfit = compute_allregion_time_dependent_linear_model(AFR_table, bin_lm_s, no_points_thresh, well_type);
clear AFR_table;
disp("fi= "+fi+" regi="+regi+" processed for lm computation")
%%
time_vec = AFR_table_w_lmfit.Properties.CustomProperties.time_vec;
time_diff = mean(unique(diff(time_vec)));
total_time_s = time_vec(end);
AFR_table_w_lmfit = addprop(AFR_table_w_lmfit, {'mdl_par_time_vec','no_points_threshold'},{'table','table'});
AFR_table_w_lmfit.Properties.CustomProperties.mdl_par_time_vec = (time_diff:bin_lm_s:total_time_s) + bin_lm_s/2 - time_diff;
AFR_table_w_lmfit.Properties.CustomProperties.no_points_threshold = no_points_thresh;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
