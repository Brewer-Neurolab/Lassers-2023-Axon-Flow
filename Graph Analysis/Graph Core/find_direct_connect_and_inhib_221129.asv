%% No Stim
clear
clc
close all
%Loading prerequisites
load("D:\Brewer lab data\HFS\Temporal Analysis\No Stim\dataInfo.mat")
load("D:\Brewer lab data\HFS\Temporal Analysis\No Stim\full_idx_allregion_unit_matched_stim.mat")

allregion_unit_matched=convert_allregion_unit_matched_220222(allregion_unit_matched_stim);

% full spike index
data_folder_addr="D:\Brewer lab data\HFS\No Stim\full_index_pseudo_times";

data_folder_dir=dir(data_folder_addr);
data_folder_isdir=[data_folder_dir.isdir];
data_folder_names=string({data_folder_dir(data_folder_isdir).name});
data_folder_names=data_folder_names(3:end)';
% data_folder_names=erase(data_folder_names,"_mat_files")';
all_region_order=[];
for i=1:length(data_folder_names)
    all_region_order(i)=find(contains(data_folder_names,dataInfo.meaName(i)));
end
allregion_unit_matched=allregion_unit_matched(all_region_order);

%no stim
source = load('D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output Full Index slide 10sLM\concatenated_source_table.mat');
target = load('D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output Full Index slide 10sLM\concatenated_target_table.mat');

NodeTable = create_node_table_200820;
regList = ["EC","DG","CA3","CA1"];
source_regList_labels = ["EC->tunnels","DG->tunnels","CA3->tunnels","CA1->tunnels"];
target_regList_labels = ["tunnels->EC","tunnels->DG","tunnels->CA3","tunnels->CA1"];
% powerlawfitfun = @(b,x) 10.^(b(2)*log10(x) + b(1)); 
powerlawfitfunc= @(b,x) b(1) + b(2)*log(x);
%% Create averaged edges
time_vec = [0.05:0.1:299.95];
no_ti = length(time_vec);
G_Slopes=[];
G_Rsq=[];
for graph_type = ["slope","rsq"]
    for fi=1:length(allregion_unit_matched)
        output_folder = "./time_depn_graph_plots/"+fi+"/"+graph_type+"/";
        if ~exist(output_folder,"dir"), mkdir(output_folder); end
        %         hold on
        G_well_cat=[];
        %G_tunnel_cat=[];
        for ti = 1:no_ti
            source_table_full = source.AFR_table(source.AFR_table.fi == fi,:);
            target_table_full = target.AFR_table(target.AFR_table.fi == fi,:);
            
            % Compute time subset
            source_subset = time_subset(source_table_full, ti, 0.1, 0.2, 0);
            target_subset = time_subset(target_table_full, ti, 0.1, 0.2, 0);
            
            if isempty(source_subset) && isempty(target_subset)
                continue;
            end
            
            % Computing Edge Table with only positive slopes
            [EdgeTable] = create_edge_table_200820(NodeTable, target_subset, source_subset, graph_type);
            
            G = digraph(EdgeTable,NodeTable);
            
            if isempty(G_well_cat)
                G_well_cat=G;
            else
                G_well_cat=addedge(G_well_cat,G.Edges);
            end
            
            % Plotting graphs
            
            %[gp,ax] = plot_graph(G,dataInfo.orientation{fi}, graph_type);
            %             if strcmp(graph_type,"slope")
            %                 plot_tunnels_in_graph_220309( allregion_unit_matched{fi},G,1)
            %             else
            %                 plot_tunnels_in_graph_220309( allregion_unit_matched{fi},G,0)
            %             end
            axis tight
            %title("Time = "+time_vec(ti)+"s", "FontSize",20)
        end
        %        close all
        %saveas(gcf,char(output_folder+graph_type+"_graph_"+strrep(string(time_vec(ti)),".","p")),'png')
        %        hold off
        
        % remove duplicate edges
        if ~isempty(G_well_cat)
            unique_edges=unique(string(G_well_cat.Edges.EndNodes),'row','stable');
            unique_edges=strcat(unique_edges(:,1),unique_edges(:,2));
            G_keep=[];
            G_wt_avg=[];
            endNodes=strcat(string(G_well_cat.Edges.EndNodes(:,1)),string(G_well_cat.Edges.EndNodes(:,2)));
            for uNodes=1:length(unique_edges)
                [max_W,~]=max(G_well_cat.Edges.Weight(strcmp(unique_edges(uNodes),endNodes)));
                %G_keep=addedge(G_keep,G_well_cat.Edges(to_keep,:));
                to_keep=find(strcmp(unique_edges(uNodes),endNodes)&G_well_cat.Edges.Weight==max_W(1));
                avg_wt=mean(G_well_cat.Edges.Weight(strcmp(unique_edges(uNodes),endNodes)));
                G_wt_avg=[G_wt_avg;avg_wt];
                G_keep=[G_keep;to_keep(1)];
            end
            to_delete=ones(size(string(G_well_cat.Edges.EndNodes),1),1);
            to_delete(G_keep)=0;
            to_delete=find(to_delete);
            G_well_cat=rmedge(G_well_cat,to_delete);
            G_well_cat.Edges.Weight=G_wt_avg;

            if graph_type=="slope"
                G_Slopes{fi}=G_well_cat;
            else
                G_Rsq{fi}=G_well_cat;
            end

        end
        disp(string(fi))
    end
    close all
end

save("D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output Full Index slide 10sLM\avg_G_Slope.mat","G_Slopes")
save("D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output Full Index slide 10sLM\avg_G_Rsq.mat","G_Rsq")
%% Load
load("D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output Full Index slide 10sLM\avg_G_Slope.mat")
load("D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output Full Index slide 10sLM\avg_G_Rsq.mat")

%% Find Percent Direct Connect Through Slopes
num_direct=0;
for i=1:length(G_Slopes)
    num_direct=num_direct+length(G_Slopes{i}.Edges.Weight(G_Slopes{i}.Edges.Weight<=1.05 ...
        & G_Slopes{i}.Edges.Weight>=0.95));
end

%% Find Strictly Inhibitory Axons, Find if consistently inhibitory

G_inhib=[];
for i=1:length(G_Slopes)
    G_inhib{i}=G_Slopes{i}.Edges(G_Slopes{i}.Edges.Weight<1,:);
end

num_inhib=sum(cellfun(@height,G_inhib));%% Finding inhibitory Axons Based on Slope
clear
clc
close all
%% HFS5

load("D:\Brewer lab data\HFS\Temporal Analysis\5 HFS\dataInfo.mat")
load("D:\Brewer lab data\HFS\Temporal Analysis\5 HFS\full_idx_allregion_unit_matched_stim.mat")

allregion_unit_matched=convert_allregion_unit_matched_220222(allregion_unit_matched_stim);

% full spike index
data_folder_addr="D:\Brewer lab data\HFS\Theta Stim\full_index_pseudo_times";

data_folder_dir=dir(data_folder_addr);
data_folder_isdir=[data_folder_dir.isdir];
data_folder_names=string({data_folder_dir(data_folder_isdir).name});
data_folder_names=data_folder_names(3:end)';
% data_folder_names=erase(data_folder_names,"_mat_files")';
all_region_order=[];
for i=1:length(data_folder_names)
    all_region_order(i)=find(contains(data_folder_names,dataInfo.meaName(i)));
end
allregion_unit_matched=allregion_unit_matched(all_region_order);

% 5 HFS
source=load("D:\Brewer lab data\HFS\Temporal Analysis\5 HFS\AFR Output Full Index slide 10sLM\concatenated_source_table.mat");
target=load("D:\Brewer lab data\HFS\Temporal Analysis\5 HFS\AFR Output Full Index slide 10sLM\concatenated_target_table.mat");

NodeTable = create_node_table_200820;
regList = ["EC","DG","CA3","CA1"];
source_regList_labels = ["EC->tunnels","DG->tunnels","CA3->tunnels","CA1->tunnels"];
target_regList_labels = ["tunnels->EC","tunnels->DG","tunnels->CA3","tunnels->CA1"];
% powerlawfitfun = @(b,x) 10.^(b(2)*log10(x) + b(1)); 
powerlawfitfunc= @(b,x) b(1) + b(2)*log(x);
%% Create averaged edges
time_vec = [0.05:0.1:299.95];
no_ti = length(time_vec);
G_Slopes=[];
G_Rsq=[];
for graph_type = ["slope","rsq"]
    for fi=1:length(allregion_unit_matched)
        output_folder = "./time_depn_graph_plots/"+fi+"/"+graph_type+"/";
        if ~exist(output_folder,"dir"), mkdir(output_folder); end
        %         hold on
        G_well_cat=[];
        %G_tunnel_cat=[];
        for ti = 1:no_ti
            source_table_full = source.AFR_table(source.AFR_table.fi == fi,:);
            target_table_full = target.AFR_table(target.AFR_table.fi == fi,:);
            
            % Compute time subset
            source_subset = time_subset(source_table_full, ti, 0.1, 0.2, 0);
            target_subset = time_subset(target_table_full, ti, 0.1, 0.2, 0);
            
            if isempty(source_subset) && isempty(target_subset)
                continue;
            end
            
            % Computing Edge Table with only positive slopes
            [EdgeTable] = create_edge_table_200820(NodeTable, target_subset, source_subset, graph_type);
            
            G = digraph(EdgeTable,NodeTable);
            
            if isempty(G_well_cat)
                G_well_cat=G;
            else
                G_well_cat=addedge(G_well_cat,G.Edges);
            end
            
            % Plotting graphs
            
            %[gp,ax] = plot_graph(G,dataInfo.orientation{fi}, graph_type);
            %             if strcmp(graph_type,"slope")
            %                 plot_tunnels_in_graph_220309( allregion_unit_matched{fi},G,1)
            %             else
            %                 plot_tunnels_in_graph_220309( allregion_unit_matched{fi},G,0)
            %             end
            axis tight
            %title("Time = "+time_vec(ti)+"s", "FontSize",20)
        end
        %        close all
        %saveas(gcf,char(output_folder+graph_type+"_graph_"+strrep(string(time_vec(ti)),".","p")),'png')
        %        hold off
        
        % remove duplicate edges
        if ~isempty(G_well_cat)
            unique_edges=unique(string(G_well_cat.Edges.EndNodes),'row','stable');
            unique_edges=strcat(unique_edges(:,1),unique_edges(:,2));
            G_keep=[];
            G_wt_avg=[];
            endNodes=strcat(string(G_well_cat.Edges.EndNodes(:,1)),string(G_well_cat.Edges.EndNodes(:,2)));
            for uNodes=1:length(unique_edges)
                [max_W,~]=max(G_well_cat.Edges.Weight(strcmp(unique_edges(uNodes),endNodes)));
                %G_keep=addedge(G_keep,G_well_cat.Edges(to_keep,:));
                to_keep=find(strcmp(unique_edges(uNodes),endNodes)&G_well_cat.Edges.Weight==max_W(1));
                avg_wt=mean(G_well_cat.Edges.Weight(strcmp(unique_edges(uNodes),endNodes)));
                G_wt_avg=[G_wt_avg;avg_wt];
                G_keep=[G_keep;to_keep(1)];
            end
            to_delete=ones(size(string(G_well_cat.Edges.EndNodes),1),1);
            to_delete(G_keep)=0;
            to_delete=find(to_delete);
            G_well_cat=rmedge(G_well_cat,to_delete);
            G_well_cat.Edges.Weight=G_wt_avg;

            if graph_type=="slope"
                G_Slopes{fi}=G_well_cat;
            else
                G_Rsq{fi}=G_well_cat;
            end

        end
        disp(string(fi))
    end
    close all
end

save("D:\Brewer lab data\HFS\Temporal Analysis\5 HFS\AFR Output Full Index slide 10sLM\avg_G_Slope.mat","G_Slopes")
save("D:\Brewer lab data\HFS\Temporal Analysis\5 HFS\AFR Output Full Index slide 10sLM\avg_G_Rsq.mat","G_Rsq")
%% Load
load("D:\Brewer lab data\HFS\Temporal Analysis\5 HFS\AFR Output Full Index slide 10sLM\avg_G_Slope.mat")
load("D:\Brewer lab data\HFS\Temporal Analysis\5 HFS\AFR Output Full Index slide 10sLM\avg_G_Rsq.mat")
%% Find Percent Direct Connect Through Slopes
num_direct=0;
for i=1:length(G_Slopes)
    num_direct=num_direct+length(G_Slopes{i}.Edges.Weight(G_Slopes{i}.Edges.Weight<=1.05 ...
        & G_Slopes{i}.Edges.Weight>=0.95));
end

%% Find Strictly Inhibitory Axons, Find if consistently inhibitory

G_inhib=[];
for i=1:length(G_Slopes)
    G_inhib{i}=G_Slopes{i}.Edges(G_Slopes{i}.Edges.Weight<1,:);
end

num_inhib=sum(cellfun(@height,G_inhib));

%% HFS40
clear
clc
close all
%Loading prerequisites
load("D:\Brewer lab data\HFS\Temporal Analysis\40 HFS\dataInfo.mat")
load("D:\Brewer lab data\HFS\Temporal Analysis\40 HFS\full_idx_allregion_unit_matched_stim.mat")

allregion_unit_matched=convert_allregion_unit_matched_220222(allregion_unit_matched_stim);

% full spike index
data_folder_addr="D:\Brewer lab data\HFS\HFS Stim\full_index_pseudo_times";

data_folder_dir=dir(data_folder_addr);
data_folder_isdir=[data_folder_dir.isdir];
data_folder_names=string({data_folder_dir(data_folder_isdir).name});
data_folder_names=data_folder_names(3:end)';
% data_folder_names=erase(data_folder_names,"_mat_files")';
all_region_order=[];
for i=1:length(data_folder_names)
    all_region_order(i)=find(contains(data_folder_names,dataInfo.meaName(i)));
end
allregion_unit_matched=allregion_unit_matched(all_region_order);

% 40 HFS
source=load("D:\Brewer lab data\HFS\Temporal Analysis\40 HFS\AFR Output Full Index slide 10sLM\concatenated_source_table.mat");
target=load("D:\Brewer lab data\HFS\Temporal Analysis\40 HFS\AFR Output Full Index slide 10sLM\concatenated_target_table.mat");

NodeTable = create_node_table_200820;
regList = ["EC","DG","CA3","CA1"];
source_regList_labels = ["EC->tunnels","DG->tunnels","CA3->tunnels","CA1->tunnels"];
target_regList_labels = ["tunnels->EC","tunnels->DG","tunnels->CA3","tunnels->CA1"];
% powerlawfitfun = @(b,x) 10.^(b(2)*log10(x) + b(1)); 
powerlawfitfunc= @(b,x) b(1) + b(2)*log(x);
%% Create averaged edges
time_vec = [0.05:0.1:299.95];
no_ti = length(time_vec);
G_Slopes=[];
G_Rsq=[];
for graph_type = ["slope","rsq"]
    for fi=1:length(allregion_unit_matched)
        output_folder = "./time_depn_graph_plots/"+fi+"/"+graph_type+"/";
        if ~exist(output_folder,"dir"), mkdir(output_folder); end
        %         hold on
        G_well_cat=[];
        %G_tunnel_cat=[];
        for ti = 1:no_ti
            source_table_full = source.AFR_table(source.AFR_table.fi == fi,:);
            target_table_full = target.AFR_table(target.AFR_table.fi == fi,:);
            
            % Compute time subset
            source_subset = time_subset(source_table_full, ti, 0.1, 0.2, 0);
            target_subset = time_subset(target_table_full, ti, 0.1, 0.2, 0);
            
            if isempty(source_subset) && isempty(target_subset)
                continue;
            end
            
            % Computing Edge Table with only positive slopes
            [EdgeTable] = create_edge_table_200820(NodeTable, target_subset, source_subset, graph_type);
            
            G = digraph(EdgeTable,NodeTable);
            
            if isempty(G_well_cat)
                G_well_cat=G;
            else
                G_well_cat=addedge(G_well_cat,G.Edges);
            end
            
            % Plotting graphs
            
            %[gp,ax] = plot_graph(G,dataInfo.orientation{fi}, graph_type);
            %             if strcmp(graph_type,"slope")
            %                 plot_tunnels_in_graph_220309( allregion_unit_matched{fi},G,1)
            %             else
            %                 plot_tunnels_in_graph_220309( allregion_unit_matched{fi},G,0)
            %             end
            axis tight
            %title("Time = "+time_vec(ti)+"s", "FontSize",20)
        end
        %        close all
        %saveas(gcf,char(output_folder+graph_type+"_graph_"+strrep(string(time_vec(ti)),".","p")),'png')
        %        hold off
        
        % remove duplicate edges
        if ~isempty(G_well_cat)
            unique_edges=unique(string(G_well_cat.Edges.EndNodes),'row','stable');
            unique_edges=strcat(unique_edges(:,1),unique_edges(:,2));
            G_keep=[];
            G_wt_avg=[];
            endNodes=strcat(string(G_well_cat.Edges.EndNodes(:,1)),string(G_well_cat.Edges.EndNodes(:,2)));
            for uNodes=1:length(unique_edges)
                [max_W,~]=max(G_well_cat.Edges.Weight(strcmp(unique_edges(uNodes),endNodes)));
                %G_keep=addedge(G_keep,G_well_cat.Edges(to_keep,:));
                to_keep=find(strcmp(unique_edges(uNodes),endNodes)&G_well_cat.Edges.Weight==max_W(1));
                avg_wt=mean(G_well_cat.Edges.Weight(strcmp(unique_edges(uNodes),endNodes)));
                G_wt_avg=[G_wt_avg;avg_wt];
                G_keep=[G_keep;to_keep(1)];
            end
            to_delete=ones(size(string(G_well_cat.Edges.EndNodes),1),1);
            to_delete(G_keep)=0;
            to_delete=find(to_delete);
            G_well_cat=rmedge(G_well_cat,to_delete);
            G_well_cat.Edges.Weight=G_wt_avg;

            if graph_type=="slope"
                G_Slopes{fi}=G_well_cat;
            else
                G_Rsq{fi}=G_well_cat;
            end

        end
        disp(string(fi))
    end
    close all
end

save("D:\Brewer lab data\HFS\Temporal Analysis\40 HFS\AFR Output Full Index slide 10sLM\avg_G_Slope.mat","G_Slopes")
save("D:\Brewer lab data\HFS\Temporal Analysis\40 HFS\AFR Output Full Index slide 10sLM\avg_G_Rsq.mat","G_Rsq")
%% Load
load("D:\Brewer lab data\HFS\Temporal Analysis\40 HFS\AFR Output Full Index slide 10sLM\avg_G_Slope.mat")
load("D:\Brewer lab data\HFS\Temporal Analysis\40 HFS\AFR Output Full Index slide 10sLM\avg_G_Rsq.mat")

%% Find Percent Direct Connect Through Slopes
num_direct=0;
for i=1:length(G_Slopes)
    num_direct=num_direct+length(G_Slopes{i}.Edges.Weight(G_Slopes{i}.Edges.Weight<=1.05 ...
        & G_Slopes{i}.Edges.Weight>=0.95));
end

%% Find Strictly Inhibitory Axons, Find if consistently inhibitory

G_inhib=[];
for i=1:length(G_Slopes)
    G_inhib{i}=G_Slopes{i}.Edges(G_Slopes{i}.Edges.Weight<1,:);
end

num_inhib=sum(cellfun(@height,G_inhib));
%% Functions
function subset_table = time_subset( table, ti, slope_thresh, rsq_thresh, if_pval_thresh)

subset_table=[];

if isempty(table)
    return
end

slope_mat = cell2mat(table.slope);
pval_mat = cell2mat(table.p_val);
rsq_mat = cell2mat(table.rsq);

slope_vec = slope_mat(:,ti);
pval_vec = pval_mat(:,ti);
rsq_vec = rsq_mat(:,ti);

if if_pval_thresh
    z_idx = isnan(slope_vec) | isnan(pval_vec) | isnan(rsq_vec) | ...
        slope_vec<rsq_thresh | rsq_vec<slope_thresh ;
else z_idx = isnan(slope_vec) | isnan(pval_vec) | isnan(rsq_vec) | ...
        slope_vec<rsq_thresh | rsq_vec<slope_thresh | pval_vec < 0.05;
end
table.slope = slope_vec;
table.p_val = pval_vec;
table.rsq = rsq_vec;

subset_table = table(~z_idx,:);


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
