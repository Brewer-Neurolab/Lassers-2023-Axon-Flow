%% Create Graphs with strength of connections varying with time representing slopes of AFR between nodes
% Reduction of degrees to generate summary figure

clear
clc
close all
%Loading prerequisites
nostim_datatable=(load("D:\Brewer lab data\HFS\Temporal Analysis\No Stim\dataInfo.mat"));
nostim_datatable=nostim_datatable.dataInfo;
nostim_allreg=(load("D:\Brewer lab data\HFS\Temporal Analysis\No Stim\full_idx_allregion_unit_matched_stim.mat"));
nostim_allreg=nostim_allreg.allregion_unit_matched_stim;

hfs5_datatable=(load("D:\Brewer lab data\HFS\Temporal Analysis\5 HFS\dataInfo.mat"));
hfs5_datatable=hfs5_datatable.dataInfo;
hfs5_allreg=(load("D:\Brewer lab data\HFS\Temporal Analysis\5 HFS\full_idx_allregion_unit_matched_stim.mat"));
hfs5_allreg=hfs5_allreg.allregion_unit_matched_stim;

hfs40_datatable=(load("D:\Brewer lab data\HFS\Temporal Analysis\40 HFS\dataInfo.mat"));
hfs40_datatable=hfs40_datatable.dataInfo;
hfs40_allreg=(load("D:\Brewer lab data\HFS\Temporal Analysis\40 HFS\full_idx_allregion_unit_matched_stim.mat"));
hfs40_allreg=hfs40_allreg.allregion_unit_matched_stim;

nostim_allregion_unit_matched=convert_allregion_unit_matched_220222(nostim_allreg);
hfs5_allregion_unit_matched=convert_allregion_unit_matched_220222(hfs5_allreg);
hfs40_allregion_unit_matched=convert_allregion_unit_matched_220222(hfs40_allreg);

% full spike index
data_folder_addr_nostim="D:\Brewer lab data\HFS\No Stim\full_index_pseudo_times";
data_folder_addr_5hfs="D:\Brewer lab data\HFS\Theta Stim\full_index_pseudo_times";
data_folder_addr_40hfs="D:\Brewer lab data\HFS\HFS Stim\full_index_pseudo_times";

nostim_allregion_unit_matched=order_allregion(nostim_datatable,nostim_allregion_unit_matched,data_folder_addr_nostim);
hfs5_allregion_unit_matched=order_allregion(hfs5_datatable, hfs5_allregion_unit_matched,data_folder_addr_5hfs);
hfs40_allregion_unit_matched=order_allregion(hfs40_datatable, hfs40_allregion_unit_matched, data_folder_addr_40hfs);

%no stim
nostim_source = (load('D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output Full Index slide 10sLM\concatenated_source_table.mat'));
nostim_target = (load('D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output Full Index slide 10sLM\concatenated_target_table.mat'));

% 5 HFS
hfs5_source=(load("D:\Brewer lab data\HFS\Temporal Analysis\5 HFS\AFR Output Full Index slide 10sLM\concatenated_source_table.mat"));
hfs5_target=(load("D:\Brewer lab data\HFS\Temporal Analysis\5 HFS\AFR Output Full Index slide 10sLM\concatenated_target_table.mat"));

% 40 HFS
hfs40_source=(load("D:\Brewer lab data\HFS\Temporal Analysis\40 HFS\AFR Output Full Index slide 10sLM\concatenated_source_table.mat"));
hfs40_target=(load("D:\Brewer lab data\HFS\Temporal Analysis\40 HFS\AFR Output Full Index slide 10sLM\concatenated_target_table.mat"));

NodeTable = create_node_table_200820;
regList = ["EC","DG","CA3","CA1"];
source_regList_labels = ["EC->tunnels","DG->tunnels","CA3->tunnels","CA1->tunnels"];
target_regList_labels = ["tunnels->EC","tunnels->DG","tunnels->CA3","tunnels->CA1"];
% powerlawfitfun = @(b,x) 10.^(b(2)*log10(x) + b(1));
powerlawfitfunc= @(b,x) b(1) + b(2)*log(x);

colors=[0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;...
    0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.3010 0.7450 0.9330;0.6350 0.0780 0.1840];

colormap(colors)
close all

% ff_colors=["#218832","#CA1A00","#EDB120"];
% ff_colors=["#3869EE","#BA0000","#EE9B00"];
% fb_colors=["#4fb7ff","#FF2303","#FFDA49"];

%old colors
% ff_colors=["#3869EE","#FF2303","#EE9B00"];
% fb_colors=["#4fb7ff","#FF5151","#FFDA49"];

%new colors
ff_colors=["#3869EE","#FF2303","#EE9B00"];
fb_colors=["#4fb7ff","#FF9B9B","#FFDA49"];

%10_2
% ff_colors=["#56A0F1","#CA1A00","#EDB120"];
% fb_colors=["#4FB7FF","#FF2303","#FFC843"];

%isolates one array
my_id=3;

% FID 3 only
% nostim_allregion_unit_matched=nostim_allregion_unit_matched(my_id);
% nostim_source.AFR_table=nostim_source.AFR_table(nostim_source.AFR_table.fi==my_id,:);
% nostim_target.AFR_table=nostim_target.AFR_table(nostim_target.AFR_table.fi==my_id,:);
% nostim_datatable=nostim_datatable(my_id,:);
% 
% hfs5_allregion_unit_matched=hfs5_allregion_unit_matched(my_id);
% hfs5_source.AFR_table=hfs5_source.AFR_table(hfs5_source.AFR_table.fi==my_id,:);
% hfs5_target.AFR_table=hfs5_target.AFR_table(hfs5_target.AFR_table.fi==my_id,:);
% hfs5_datatable=hfs5_datatable(my_id,:);
% 
% hfs40_allregion_unit_matched=hfs40_allregion_unit_matched(my_id);
% hfs40_source.AFR_table=hfs40_source.AFR_table(hfs40_source.AFR_table.fi==my_id,:);
% hfs40_target.AFR_table=hfs40_target.AFR_table(hfs40_target.AFR_table.fi==my_id,:);
% hfs40_datatable=hfs40_datatable(my_id,:);
%% Graph Analysis Calc
% No need to run this section if these data structures already exist

[nostim_slope_G,nostim_rsq_G]=cat_edges(nostim_allregion_unit_matched,nostim_source,nostim_target,NodeTable,nostim_datatable,my_id);
disp("No stim computed.")
[hfs5_slope_G,hfs5_rsq_G]=cat_edges(hfs5_allregion_unit_matched,hfs5_source,hfs5_target,NodeTable,hfs5_datatable,my_id);
disp("5 HFS computed.")
[hfs40_slope_G,hfs40_rsq_G]=cat_edges(hfs5_allregion_unit_matched,hfs40_source,hfs40_target,NodeTable,hfs40_datatable,my_id);
disp("40 HFS computed.")
% 
% save("D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output Full Index slide 10sLM\nostim_slope_G.mat","nostim_slope_G")
% save("D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output Full Index slide 10sLM\nostim_rsq_G.mat","nostim_rsq_G")
% 
% save("D:\Brewer lab data\HFS\Temporal Analysis\5 HFS\AFR Output Full Index slide 10sLM\hfs5_slope_G.mat","hfs5_slope_G")
% save("D:\Brewer lab data\HFS\Temporal Analysis\5 HFS\AFR Output Full Index slide 10sLM\hfs5_rsq_G.mat","hfs5_rsq_G")
% 
% save("D:\Brewer lab data\HFS\Temporal Analysis\40 HFS\AFR Output Full Index slide 10sLM\hfs40_slope_G.mat","hfs40_slope_G")
% save("D:\Brewer lab data\HFS\Temporal Analysis\40 HFS\AFR Output Full Index slide 10sLM\hfs40_rsq_G.mat","hfs40_rsq_G")

%% Related functions

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

function allregion_unit_matched=order_allregion(dataInfo,allreg,data_folder_addr)

data_folder_dir=dir(data_folder_addr);
data_folder_isdir=[data_folder_dir.isdir];
data_folder_names=string({data_folder_dir(data_folder_isdir).name});
data_folder_names=data_folder_names(3:end)';
%data_folder_names=erase(data_folder_names,"_mat_files")';
all_region_order=[];
for i=1:length(data_folder_names)
    all_region_order(i)=find(contains(data_folder_names,dataInfo.meaName(i)));
end
allregion_unit_matched=allreg(all_region_order);

end

function [G_cat_slopes,G_cat_rsq]=cat_edges(allregion_unit_matched,source,target,NodeTable,dataInfo,my_id)
% does not save tunnel-tunnel interactions in G for slopes or rsq
G_cat_slopes=[];
G_cat_rsq=[];

time_vec = [0.05:0.1:299.95];
no_ti = length(time_vec);
for graph_type = ["slope","rsq"]
    for fi=my_id
        %         output_folder = "./time_depn_graph_plots/"+fi+"/"+graph_type+"/";
        %         if ~exist(output_folder,"dir"), mkdir(output_folder); end
        figure('Renderer','painters','Position',[10 10 800 800])
        %         hold on
        G_well_cat=[];
        %G_tunnel_cat=[];
        for ti = 1:no_ti
            source_table_full = source.AFR_table(source.AFR_table.fi == fi,:);
            target_table_full = target.AFR_table(target.AFR_table.fi == fi,:);

            % Compute time subset
            source_subset = time_subset(source_table_full, ti, 0.1, 0.35, 0); %increased threshold to 0.35 from 0.2
            target_subset = time_subset(target_table_full, ti, 0.1, 0.35, 0); %increased threshold to 0.35 from 0.2

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

            [gp,ax] = plot_graph(G_well_cat,dataInfo.orientation{fi}, graph_type);
            if strcmp(graph_type,"slope")
                plot_tunnels_in_graph_220309( allregion_unit_matched{fi},G_well_cat,1)
                G_cat_slopes{fi}=G_well_cat;
            else
                plot_tunnels_in_graph_220309( allregion_unit_matched{fi},G_well_cat,0)
                G_cat_rsq{fi}=G_well_cat;
            end
            %title(graph_type+" FID "+string(fi))
            set(ax,'position',[0 0 1 1])
            %             saveas(gcf,char(output_folder+graph_type+"_graph_full_time"),'png')
        end
        disp(string(fi))
    end
%     close all
end

end

function [well_out_ff,well_out_fb,tunnel_out_ff,tunnel_out_fb]...
    =get_degrees(G,datatable,data_fldr,temporal_fldr)
well_out_ff=[];
tunnel_out_ff=[];
well_out_fb=[];
tunnel_out_fb=[];
regions=["EC","DG","CA3","CA1"];
% tunnel_reg=["EC","DG","CA3","CA1"];
for fi=1:length(G)
    Edges_ff=G{fi}.Edges(G{fi}.Edges.Code==0.1,:);
    Edges_fb=G{fi}.Edges(G{fi}.Edges.Code==0.6,:);

    G_ff=digraph(Edges_ff,G{fi}.Nodes);
    G_fb=digraph(Edges_fb,G{fi}.Nodes);

    out_degree_ff=outdegree(G_ff);
    out_degree_fb=outdegree(G_fb);

    for regi=1:4
        mea=load(data_fldr+'\mea_'+regions(regi)+'.mat');
        if datatable.orientation(fi)=="CW"
            well_out_ff{fi,regi}=out_degree_ff(ismember(G_ff.Nodes.Name,...
                mea.cw.channel_names));
            well_out_fb{fi,regi}=out_degree_fb(ismember(G_fb.Nodes.Name,...
                mea.cw.channel_names));

            well_elec_ff{fi,regi}=G_ff.Nodes.Name(ismember(G_ff.Nodes.Name,...
                mea.cw.channel_names));
            well_elec_fb{fi,regi}=G_fb.Nodes.Name(ismember(G_fb.Nodes.Name,...
                mea.cw.channel_names));

            tunnels=load(data_fldr+'\cw_vars_tunnel.mat');
            load(temporal_fldr+'\tunnelChannelList.mat');
            ff_tunnel_List=circshift(tunnelChannelList.near_target.cw,1);
            tun_elec_ff=tunnels.channel_names(ff_tunnel_List{regi});
            fb_tunnel_List=circshift(tunnelChannelList.near_source.cw,0);
            tun_elec_fb=tunnels.channel_names(fb_tunnel_List{regi});

            tunnel_out_ff{fi,regi}=out_degree_ff(ismember(G_ff.Nodes.Name,...
                tun_elec_ff));
            tunnel_ff_elec{fi,regi}=tun_elec_ff;
            tunnel_ff_orientation{fi,regi}="CW";

            tunnel_out_fb{fi,regi}=out_degree_fb(ismember(G_fb.Nodes.Name,...
                tun_elec_fb));
            tunnel_fb_elec{fi,regi}=tun_elec_fb;
            tunnel_fb_orientation{fi,regi}="CW";
        else
            well_out_ff{fi,regi}=out_degree_ff(ismember(G_ff.Nodes.Name,...
                mea.ccw.channel_names));
            well_out_fb{fi,regi}=out_degree_fb(ismember(G_fb.Nodes.Name,...
                mea.ccw.channel_names));

            well_elec_ff{fi,regi}=G_ff.Nodes.Name(ismember(G_ff.Nodes.Name,...
                mea.ccw.channel_names));
            well_elec_fb{fi,regi}=G_fb.Nodes.Name(ismember(G_fb.Nodes.Name,...
                mea.ccw.channel_names));

            load(temporal_fldr+'\tunnelChannelList.mat');
            tunnels=load(data_fldr+'\cw_vars_tunnel.mat');
            load(temporal_fldr+'\tunnelChannelList.mat');
            ff_tunnel_List=circshift(tunnelChannelList.near_target.ccw,1);
            tun_elec_ff=tunnels.channel_names(ff_tunnel_List{regi});
            fb_tunnel_List=circshift(tunnelChannelList.near_source.ccw,0);
            tun_elec_fb=tunnels.channel_names(fb_tunnel_List{regi});

            tunnel_out_ff{fi,regi}=out_degree_ff(ismember(G_ff.Nodes.Name,...
                tun_elec_ff));
            tunnel_ff_elec{fi,regi}=tun_elec_ff;
            tunnel_ff_orientation{fi,regi}="CCW";

            tunnel_out_fb{fi,regi}=out_degree_fb(ismember(G_fb.Nodes.Name,...
                tun_elec_fb));
            tunnel_fb_elec{fi,regi}=tun_elec_fb;
            tunnel_fb_orientation{fi,regi}="CCW";
        end
    end
end

well_out_ff=cell2table([well_out_ff,well_elec_ff,tunnel_ff_orientation(:,1)],...
    "VariableNames",["EC","DG","CA3","CA1","EC Elec","DG Elec","CA3 Elec","CA1 Elec","Orientation"]);
well_out_fb=cell2table([well_out_fb,well_elec_fb,tunnel_fb_orientation(:,1)],...
    "VariableNames",["EC","DG","CA3","CA1","EC Elec","DG Elec","CA3 Elec","CA1 Elec","Orientation"]);
tunnel_out_ff=cell2table([tunnel_out_ff,tunnel_ff_elec,tunnel_ff_orientation(:,1)],...
    "VariableNames",["EC-DG","DG-CA3","CA3-CA1","CA1-EC",...
    "EC-DG Elec","DG-CA3 Elec","CA3-CA1 Elec","CA1-EC Elec","Orientation"]);
tunnel_out_fb=cell2table([tunnel_out_fb,tunnel_fb_elec,tunnel_fb_orientation(:,1)],...
    "VariableNames",["EC-DG","DG-CA3","CA3-CA1","CA1-EC",...
    "EC-DG Elec","DG-CA3 Elec","CA3-CA1 Elec","CA1-EC Elec","Orientation"]);
end

function [well_out_ff,well_out_fb,tunnel_out_ff,tunnel_out_fb]...
    =get_out_degrees(G,datatable,data_fldr,temporal_fldr)
well_out_ff=[];
tunnel_out_ff=[];
well_out_fb=[];
tunnel_out_fb=[];
regions=["EC","DG","CA3","CA1"];
% tunnel_reg=["EC","DG","CA3","CA1"];
for fi=1:length(G)
    Edges_ff=G{fi}.Edges(G{fi}.Edges.Code==0.1,:);
    Edges_fb=G{fi}.Edges(G{fi}.Edges.Code==0.6,:);

    G_ff=digraph(Edges_ff,G{fi}.Nodes);
    G_fb=digraph(Edges_fb,G{fi}.Nodes);

    out_degree_ff=outdegree(G_ff);
    out_degree_fb=outdegree(G_fb);

    for regi=1:4
        mea=load(data_fldr+'\mea_'+regions(regi)+'.mat');
        if datatable.orientation(fi)=="CW"
            well_out_ff{fi,regi}=out_degree_ff(ismember(G_ff.Nodes.Name,...
                mea.cw.channel_names));
            well_out_fb{fi,regi}=out_degree_fb(ismember(G_fb.Nodes.Name,...
                mea.cw.channel_names));

            well_elec_ff{fi,regi}=G_ff.Nodes.Name(ismember(G_ff.Nodes.Name,...
                mea.cw.channel_names));
            well_elec_fb{fi,regi}=G_fb.Nodes.Name(ismember(G_fb.Nodes.Name,...
                mea.cw.channel_names));

            tunnels=load(data_fldr+'\cw_vars_tunnel.mat');
            load(temporal_fldr+'\tunnelChannelList.mat');
            ff_tunnel_List=circshift(tunnelChannelList.near_target.cw,1);
            tun_elec_ff=tunnels.channel_names(ff_tunnel_List{regi});
            fb_tunnel_List=circshift(tunnelChannelList.near_source.cw,0);
            tun_elec_fb=tunnels.channel_names(fb_tunnel_List{regi});

            tunnel_out_ff{fi,regi}=out_degree_ff(ismember(G_ff.Nodes.Name,...
                tun_elec_ff));
            tunnel_ff_elec{fi,regi}=tun_elec_ff;
            tunnel_ff_orientation{fi,regi}="CW";

            tunnel_out_fb{fi,regi}=out_degree_fb(ismember(G_fb.Nodes.Name,...
                tun_elec_fb));
            tunnel_fb_elec{fi,regi}=tun_elec_fb;
            tunnel_fb_orientation{fi,regi}="CW";
        else
            well_out_ff{fi,regi}=out_degree_ff(ismember(G_ff.Nodes.Name,...
                mea.ccw.channel_names));
            well_out_fb{fi,regi}=out_degree_fb(ismember(G_fb.Nodes.Name,...
                mea.ccw.channel_names));

            well_elec_ff{fi,regi}=G_ff.Nodes.Name(ismember(G_ff.Nodes.Name,...
                mea.ccw.channel_names));
            well_elec_fb{fi,regi}=G_fb.Nodes.Name(ismember(G_fb.Nodes.Name,...
                mea.ccw.channel_names));

            load(temporal_fldr+'\tunnelChannelList.mat');
            tunnels=load(data_fldr+'\cw_vars_tunnel.mat');
            load(temporal_fldr+'\tunnelChannelList.mat');
            ff_tunnel_List=circshift(tunnelChannelList.near_target.ccw,1);
            tun_elec_ff=tunnels.channel_names(ff_tunnel_List{regi});
            fb_tunnel_List=circshift(tunnelChannelList.near_source.ccw,0);
            tun_elec_fb=tunnels.channel_names(fb_tunnel_List{regi});

            tunnel_out_ff{fi,regi}=out_degree_ff(ismember(G_ff.Nodes.Name,...
                tun_elec_ff));
            tunnel_ff_elec{fi,regi}=tun_elec_ff;
            tunnel_ff_orientation{fi,regi}="CCW";

            tunnel_out_fb{fi,regi}=out_degree_fb(ismember(G_fb.Nodes.Name,...
                tun_elec_fb));
            tunnel_fb_elec{fi,regi}=tun_elec_fb;
            tunnel_fb_orientation{fi,regi}="CCW";
        end
    end
end

well_out_ff=cell2table([well_out_ff,well_elec_ff,tunnel_ff_orientation(:,1)],...
    "VariableNames",["EC","DG","CA3","CA1","EC Elec","DG Elec","CA3 Elec","CA1 Elec","Orientation"]);
well_out_fb=cell2table([well_out_fb,well_elec_fb,tunnel_fb_orientation(:,1)],...
    "VariableNames",["EC","DG","CA3","CA1","EC Elec","DG Elec","CA3 Elec","CA1 Elec","Orientation"]);
tunnel_out_ff=cell2table([tunnel_out_ff,tunnel_ff_elec,tunnel_ff_orientation(:,1)],...
    "VariableNames",["EC-DG","DG-CA3","CA3-CA1","CA1-EC",...
    "EC-DG Elec","DG-CA3 Elec","CA3-CA1 Elec","CA1-EC Elec","Orientation"]);
tunnel_out_fb=cell2table([tunnel_out_fb,tunnel_fb_elec,tunnel_fb_orientation(:,1)],...
    "VariableNames",["EC-DG","DG-CA3","CA3-CA1","CA1-EC",...
    "EC-DG Elec","DG-CA3 Elec","CA3-CA1 Elec","CA1-EC Elec","Orientation"]);
end

function [well_out_ff,well_out_fb,tunnel_out_ff,tunnel_out_fb]...
    =get_in_degrees(G,datatable,data_fldr,temporal_fldr)
well_out_ff=[];
tunnel_out_ff=[];
well_out_fb=[];
tunnel_out_fb=[];
regions=["EC","DG","CA3","CA1"];
% tunnel_reg=["EC","DG","CA3","CA1"];
for fi=1:length(G)
    Edges_ff=G{fi}.Edges(G{fi}.Edges.Code==0.1,:);
    Edges_fb=G{fi}.Edges(G{fi}.Edges.Code==0.6,:);

    G_ff=digraph(Edges_ff,G{fi}.Nodes);
    G_fb=digraph(Edges_fb,G{fi}.Nodes);

    out_degree_ff=indegree(G_ff);
    out_degree_fb=indegree(G_fb);

    for regi=1:4
        mea=load(data_fldr+'\mea_'+regions(regi)+'.mat');
        if datatable.orientation(fi)=="CW"
            well_out_ff{fi,regi}=out_degree_ff(ismember(G_ff.Nodes.Name,...
                mea.cw.channel_names));
            well_out_fb{fi,regi}=out_degree_fb(ismember(G_fb.Nodes.Name,...
                mea.cw.channel_names));

            well_elec_ff{fi,regi}=G_ff.Nodes.Name(ismember(G_ff.Nodes.Name,...
                mea.cw.channel_names));
            well_elec_fb{fi,regi}=G_fb.Nodes.Name(ismember(G_fb.Nodes.Name,...
                mea.cw.channel_names));

            tunnels=load(data_fldr+'\cw_vars_tunnel.mat');
            load(temporal_fldr+'\tunnelChannelList.mat');
            ff_tunnel_List=circshift(tunnelChannelList.near_source.cw,1);
            tun_elec_ff=tunnels.channel_names(ff_tunnel_List{regi});
            fb_tunnel_List=circshift(tunnelChannelList.near_target.cw,0);
            tun_elec_fb=tunnels.channel_names(fb_tunnel_List{regi});

            tunnel_out_ff{fi,regi}=out_degree_ff(ismember(G_ff.Nodes.Name,...
                tun_elec_ff));
            tunnel_ff_elec{fi,regi}=tun_elec_ff;
            tunnel_ff_orientation{fi,regi}="CW";

            tunnel_out_fb{fi,regi}=out_degree_fb(ismember(G_fb.Nodes.Name,...
                tun_elec_fb));
            tunnel_fb_elec{fi,regi}=tun_elec_fb;
            tunnel_fb_orientation{fi,regi}="CW";
        else
            well_out_ff{fi,regi}=out_degree_ff(ismember(G_ff.Nodes.Name,...
                mea.ccw.channel_names));
            well_out_fb{fi,regi}=out_degree_fb(ismember(G_fb.Nodes.Name,...
                mea.ccw.channel_names));

            well_elec_ff{fi,regi}=G_ff.Nodes.Name(ismember(G_ff.Nodes.Name,...
                mea.ccw.channel_names));
            well_elec_fb{fi,regi}=G_fb.Nodes.Name(ismember(G_fb.Nodes.Name,...
                mea.ccw.channel_names));

            load(temporal_fldr+'\tunnelChannelList.mat');
            tunnels=load(data_fldr+'\cw_vars_tunnel.mat');
            load(temporal_fldr+'\tunnelChannelList.mat');
            ff_tunnel_List=circshift(tunnelChannelList.near_source.ccw,1);
            tun_elec_ff=tunnels.channel_names(ff_tunnel_List{regi});
            fb_tunnel_List=circshift(tunnelChannelList.near_target.ccw,0);
            tun_elec_fb=tunnels.channel_names(fb_tunnel_List{regi});

            tunnel_out_ff{fi,regi}=out_degree_ff(ismember(G_ff.Nodes.Name,...
                tun_elec_ff));
            tunnel_ff_elec{fi,regi}=tun_elec_ff;
            tunnel_ff_orientation{fi,regi}="CCW";

            tunnel_out_fb{fi,regi}=out_degree_fb(ismember(G_fb.Nodes.Name,...
                tun_elec_fb));
            tunnel_fb_elec{fi,regi}=tun_elec_fb;
            tunnel_fb_orientation{fi,regi}="CCW";
        end
    end
end

well_out_ff=cell2table([well_out_ff,well_elec_ff,tunnel_ff_orientation(:,1)],...
    "VariableNames",["EC","DG","CA3","CA1","EC Elec","DG Elec","CA3 Elec","CA1 Elec","Orientation"]);
well_out_fb=cell2table([well_out_fb,well_elec_fb,tunnel_fb_orientation(:,1)],...
    "VariableNames",["EC","DG","CA3","CA1","EC Elec","DG Elec","CA3 Elec","CA1 Elec","Orientation"]);
tunnel_out_ff=cell2table([tunnel_out_ff,tunnel_ff_elec,tunnel_ff_orientation(:,1)],...
    "VariableNames",["EC-DG","DG-CA3","CA3-CA1","CA1-EC",...
    "EC-DG Elec","DG-CA3 Elec","CA3-CA1 Elec","CA1-EC Elec","Orientation"]);
tunnel_out_fb=cell2table([tunnel_out_fb,tunnel_fb_elec,tunnel_fb_orientation(:,1)],...
    "VariableNames",["EC-DG","DG-CA3","CA3-CA1","CA1-EC",...
    "EC-DG Elec","DG-CA3 Elec","CA3-CA1 Elec","CA1-EC Elec","Orientation"]);
end

function [well_out_ff,well_out_fb,tunnel_out_ff,tunnel_out_fb]...
    =get_edgecount_out(G,datatable,data_fldr,temporal_fldr)
well_out_ff=[];
tunnel_out_ff=[];
well_out_fb=[];
tunnel_out_fb=[];
regions=["EC","DG","CA3","CA1"];
% tunnel_reg=["EC","DG","CA3","CA1"];
for fi=1:length(G)
    Edges_ff=G{fi}.Edges(G{fi}.Edges.Code==0.1,:);
    Edges_fb=G{fi}.Edges(G{fi}.Edges.Code==0.6,:);

    G_ff=digraph(Edges_ff,G{fi}.Nodes);
    G_fb=digraph(Edges_fb,G{fi}.Nodes);

    %     s_ff=string(G_ff.Edges.EndNodes(:,1));
    %     t_ff=string(G_ff.Edges.EndNodes(:,2));
    %     s_fb=string(G_fb.Edges.EndNodes(:,1));
    %     t_fb=string(G_fb.Edges.EndNodes(:,2));

    %     [s_ff,t_ff]=findedge(G_ff);
    %     [s_fb,t_fb]=findedge(G_fb);

    %     out_degree_ff=edgecount(G_ff,s_ff,t_ff);
    %     out_degree_fb=edgecount(G_fb,s_fb,t_fb);

    num_edge_ff_out=zeros(length(G_ff.Nodes.Name),1);
    %     num_edge_ff_in=zeros(length(G_ff.Nodes.Name),1);
    for i=1:length(G_ff.Nodes.Name)
        %         for j=1:length(G_ff.Edges.EndNodes(:,1))
        %             if any(G_ff.Edges.EndNodes(:,1)==G_ff.Nodes.Name(i))
        %                 num_edge_ff_out(i)=sum(edgecount(G_ff,i,1:numnodes(G_ff)));
        %             end
        %
        %             if any(G_ff.Edges.EndNodes(:,2)==G_ff.Nodes.Name(i))
        %                 num_edge_ff_in(i)=sum(edgecount(G_ff,1:numnodes(G_ff),i));
        %             end
        %         end
        num_edge_ff_out(i)=length(outedges(G_ff,i));
        %         num_edge_ff_in=inedges(G_ff,i);
    end

    num_edge_fb_out=zeros(length(G_fb.Nodes.Name),1);
    %     num_edge_fb_in=zeros(length(G_fb.Nodes.Name),1);
    for i=1:length(G_fb.Nodes.Name)
        %         for j=1:length(G_fb.Edges.EndNodes(:,1))
        %             if any(G_fb.Edges.EndNodes(:,1)==G_fb.Nodes.Name(i))
        %                 num_edge_fb_out(i)=sum(edgecount(G_fb,i,1:numnodes(G_fb)));
        %             end
        %
        %             if any(G_ff.Edges.EndNodes(:,2)==G_ff.Nodes.Name(i))
        %                 num_edge_fb_in(i)=sum(edgecount(G_fb,1:numnodes(G_fb),i));
        %             end
        %         end
        num_edge_fb_out(i)=length(outedges(G_fb,i));
        %         num_edge_fb_in=inedges(G_fb,i);
    end

    for regi=1:4
        mea=load(data_fldr+'\mea_'+regions(regi)+'.mat');
        if datatable.orientation(fi)=="CW"
            well_out_ff{fi,regi}=num_edge_ff_out(ismember(G_ff.Nodes.Name,...
                mea.cw.channel_names));
            well_out_fb{fi,regi}=num_edge_fb_out(ismember(G_fb.Nodes.Name,...
                mea.cw.channel_names));

            well_elec_ff{fi,regi}=G_ff.Nodes.Name(ismember(G_ff.Nodes.Name,...
                mea.cw.channel_names));
            well_elec_fb{fi,regi}=G_fb.Nodes.Name(ismember(G_fb.Nodes.Name,...
                mea.cw.channel_names));

            tunnels=load(data_fldr+'\cw_vars_tunnel.mat');
            load(temporal_fldr+'\tunnelChannelList.mat');
            ff_tunnel_List=circshift(tunnelChannelList.near_target.cw,1);
            tun_elec_ff=tunnels.channel_names(ff_tunnel_List{regi});
            fb_tunnel_List=circshift(tunnelChannelList.near_source.cw,0);
            tun_elec_fb=tunnels.channel_names(fb_tunnel_List{regi});

            tunnel_out_ff{fi,regi}=num_edge_ff_out(ismember(G_ff.Nodes.Name,...
                tun_elec_ff));
            tunnel_ff_elec{fi,regi}=tun_elec_ff;
            tunnel_ff_orientation{fi,regi}="CW";

            tunnel_out_fb{fi,regi}=num_edge_fb_out(ismember(G_fb.Nodes.Name,...
                tun_elec_fb));
            tunnel_fb_elec{fi,regi}=tun_elec_fb;
            tunnel_fb_orientation{fi,regi}="CW";
        else
            well_out_ff{fi,regi}=num_edge_ff_out(ismember(G_ff.Nodes.Name,...
                mea.ccw.channel_names));
            well_out_fb{fi,regi}=num_edge_fb_out(ismember(G_fb.Nodes.Name,...
                mea.ccw.channel_names));

            well_elec_ff{fi,regi}=G_ff.Nodes.Name(ismember(G_ff.Nodes.Name,...
                mea.ccw.channel_names));
            well_elec_fb{fi,regi}=G_fb.Nodes.Name(ismember(G_fb.Nodes.Name,...
                mea.ccw.channel_names));

            load(temporal_fldr+'\tunnelChannelList.mat');
            tunnels=load(data_fldr+'\cw_vars_tunnel.mat');
            load(temporal_fldr+'\tunnelChannelList.mat');
            ff_tunnel_List=circshift(tunnelChannelList.near_target.ccw,1);
            tun_elec_ff=tunnels.channel_names(ff_tunnel_List{regi});
            fb_tunnel_List=circshift(tunnelChannelList.near_source.ccw,0);
            tun_elec_fb=tunnels.channel_names(fb_tunnel_List{regi});

            tunnel_out_ff{fi,regi}=num_edge_ff_out(ismember(G_ff.Nodes.Name,...
                tun_elec_ff));
            tunnel_ff_elec{fi,regi}=tun_elec_ff;
            tunnel_ff_orientation{fi,regi}="CCW";

            tunnel_out_fb{fi,regi}=num_edge_fb_out(ismember(G_fb.Nodes.Name,...
                tun_elec_fb));
            tunnel_fb_elec{fi,regi}=tun_elec_fb;
            tunnel_fb_orientation{fi,regi}="CCW";
        end
    end
end

well_out_ff=cell2table([well_out_ff,well_elec_ff,tunnel_ff_orientation(:,1)],...
    "VariableNames",["EC","DG","CA3","CA1","EC Elec","DG Elec","CA3 Elec","CA1 Elec","Orientation"]);
well_out_fb=cell2table([well_out_fb,well_elec_fb,tunnel_fb_orientation(:,1)],...
    "VariableNames",["EC","DG","CA3","CA1","EC Elec","DG Elec","CA3 Elec","CA1 Elec","Orientation"]);
tunnel_out_ff=cell2table([tunnel_out_ff,tunnel_ff_elec,tunnel_ff_orientation(:,1)],...
    "VariableNames",["EC-DG","DG-CA3","CA3-CA1","CA1-EC",...
    "EC-DG Elec","DG-CA3 Elec","CA3-CA1 Elec","CA1-EC Elec","Orientation"]);
tunnel_out_fb=cell2table([tunnel_out_fb,tunnel_fb_elec,tunnel_fb_orientation(:,1)],...
    "VariableNames",["EC-DG","DG-CA3","CA3-CA1","CA1-EC",...
    "EC-DG Elec","DG-CA3 Elec","CA3-CA1 Elec","CA1-EC Elec","Orientation"]);
end

function [well_out_ff,well_out_fb,tunnel_out_ff,tunnel_out_fb]...
    =get_edgecount_in(G,datatable,data_fldr,temporal_fldr)
well_out_ff=[];
tunnel_out_ff=[];
well_out_fb=[];
tunnel_out_fb=[];
regions=["EC","DG","CA3","CA1"];
% tunnel_reg=["EC","DG","CA3","CA1"];
for fi=1:length(G)
    Edges_ff=G{fi}.Edges(G{fi}.Edges.Code==0.1,:);
    Edges_fb=G{fi}.Edges(G{fi}.Edges.Code==0.6,:);

    G_ff=digraph(Edges_ff,G{fi}.Nodes);
    G_fb=digraph(Edges_fb,G{fi}.Nodes);

    %     s_ff=string(G_ff.Edges.EndNodes(:,1));
    %     t_ff=string(G_ff.Edges.EndNodes(:,2));
    %     s_fb=string(G_fb.Edges.EndNodes(:,1));
    %     t_fb=string(G_fb.Edges.EndNodes(:,2));

    %     [s_ff,t_ff]=findedge(G_ff);
    %     [s_fb,t_fb]=findedge(G_fb);

    %     out_degree_ff=edgecount(G_ff,s_ff,t_ff);
    %     out_degree_fb=edgecount(G_fb,s_fb,t_fb);

    num_edge_ff_out=zeros(length(G_ff.Nodes.Name),1);
    %     num_edge_ff_in=zeros(length(G_ff.Nodes.Name),1);
    for i=1:length(G_ff.Nodes.Name)
        %         for j=1:length(G_ff.Edges.EndNodes(:,1))
        %             if any(G_ff.Edges.EndNodes(:,1)==G_ff.Nodes.Name(i))
        %                 num_edge_ff_out(i)=sum(edgecount(G_ff,i,1:numnodes(G_ff)));
        %             end
        %
        %             if any(G_ff.Edges.EndNodes(:,2)==G_ff.Nodes.Name(i))
        %                 num_edge_ff_in(i)=sum(edgecount(G_ff,1:numnodes(G_ff),i));
        %             end
        %         end
        num_edge_ff_out(i)=length(inedges(G_ff,i));
        %         num_edge_ff_in=inedges(G_ff,i);
    end

    num_edge_fb_out=zeros(length(G_fb.Nodes.Name),1);
    %     num_edge_fb_in=zeros(length(G_fb.Nodes.Name),1);
    for i=1:length(G_fb.Nodes.Name)
        %         for j=1:length(G_fb.Edges.EndNodes(:,1))
        %             if any(G_fb.Edges.EndNodes(:,1)==G_fb.Nodes.Name(i))
        %                 num_edge_fb_out(i)=sum(edgecount(G_fb,i,1:numnodes(G_fb)));
        %             end
        %
        %             if any(G_ff.Edges.EndNodes(:,2)==G_ff.Nodes.Name(i))
        %                 num_edge_fb_in(i)=sum(edgecount(G_fb,1:numnodes(G_fb),i));
        %             end
        %         end
        num_edge_fb_out(i)=length(inedges(G_fb,i));
        %         num_edge_fb_in=inedges(G_fb,i);
    end

    for regi=1:4
        mea=load(data_fldr+'\mea_'+regions(regi)+'.mat');
        if datatable.orientation(fi)=="CW"
            well_out_ff{fi,regi}=num_edge_ff_out(ismember(G_ff.Nodes.Name,...
                mea.cw.channel_names));
            well_out_fb{fi,regi}=num_edge_fb_out(ismember(G_fb.Nodes.Name,...
                mea.cw.channel_names));

            well_elec_ff{fi,regi}=G_ff.Nodes.Name(ismember(G_ff.Nodes.Name,...
                mea.cw.channel_names));
            well_elec_fb{fi,regi}=G_fb.Nodes.Name(ismember(G_fb.Nodes.Name,...
                mea.cw.channel_names));

            tunnels=load(data_fldr+'\cw_vars_tunnel.mat');
            load(temporal_fldr+'\tunnelChannelList.mat');
            ff_tunnel_List=circshift(tunnelChannelList.near_source.cw,0);
            tun_elec_ff=tunnels.channel_names(ff_tunnel_List{regi});
            fb_tunnel_List=circshift(tunnelChannelList.near_target.cw,1);
            tun_elec_fb=tunnels.channel_names(fb_tunnel_List{regi});

            tunnel_out_ff{fi,regi}=num_edge_ff_out(ismember(G_ff.Nodes.Name,...
                tun_elec_ff));
            tunnel_ff_elec{fi,regi}=tun_elec_ff;
            tunnel_ff_orientation{fi,regi}="CW";

            tunnel_out_fb{fi,regi}=num_edge_fb_out(ismember(G_fb.Nodes.Name,...
                tun_elec_fb));
            tunnel_fb_elec{fi,regi}=tun_elec_fb;
            tunnel_fb_orientation{fi,regi}="CW";
        else
            well_out_ff{fi,regi}=num_edge_ff_out(ismember(G_ff.Nodes.Name,...
                mea.ccw.channel_names));
            well_out_fb{fi,regi}=num_edge_fb_out(ismember(G_fb.Nodes.Name,...
                mea.ccw.channel_names));

            well_elec_ff{fi,regi}=G_ff.Nodes.Name(ismember(G_ff.Nodes.Name,...
                mea.ccw.channel_names));
            well_elec_fb{fi,regi}=G_fb.Nodes.Name(ismember(G_fb.Nodes.Name,...
                mea.ccw.channel_names));

            load(temporal_fldr+'\tunnelChannelList.mat');
            tunnels=load(data_fldr+'\cw_vars_tunnel.mat');
            load(temporal_fldr+'\tunnelChannelList.mat');
            ff_tunnel_List=circshift(tunnelChannelList.near_source.ccw,0);
            tun_elec_ff=tunnels.channel_names(ff_tunnel_List{regi});
            fb_tunnel_List=circshift(tunnelChannelList.near_target.ccw,1);
            tun_elec_fb=tunnels.channel_names(fb_tunnel_List{regi});

            tunnel_out_ff{fi,regi}=num_edge_ff_out(ismember(G_ff.Nodes.Name,...
                tun_elec_ff));
            tunnel_ff_elec{fi,regi}=tun_elec_ff;
            tunnel_ff_orientation{fi,regi}="CCW";

            tunnel_out_fb{fi,regi}=num_edge_fb_out(ismember(G_fb.Nodes.Name,...
                tun_elec_fb));
            tunnel_fb_elec{fi,regi}=tun_elec_fb;
            tunnel_fb_orientation{fi,regi}="CCW";
        end
    end
end

well_out_ff=cell2table([well_out_ff,well_elec_ff,tunnel_ff_orientation(:,1)],...
    "VariableNames",["EC","DG","CA3","CA1","EC Elec","DG Elec","CA3 Elec","CA1 Elec","Orientation"]);
well_out_fb=cell2table([well_out_fb,well_elec_fb,tunnel_fb_orientation(:,1)],...
    "VariableNames",["EC","DG","CA3","CA1","EC Elec","DG Elec","CA3 Elec","CA1 Elec","Orientation"]);
tunnel_out_ff=cell2table([tunnel_out_ff,tunnel_ff_elec,tunnel_ff_orientation(:,1)],...
    "VariableNames",["EC-DG","DG-CA3","CA3-CA1","CA1-EC",...
    "EC-DG Elec","DG-CA3 Elec","CA3-CA1 Elec","CA1-EC Elec","Orientation"]);
tunnel_out_fb=cell2table([tunnel_out_fb,tunnel_fb_elec,tunnel_fb_orientation(:,1)],...
    "VariableNames",["EC-DG","DG-CA3","CA3-CA1","CA1-EC",...
    "EC-DG Elec","DG-CA3 Elec","CA3-CA1 Elec","CA1-EC Elec","Orientation"]);
end

function [well_out_ff,well_out_fb,tunnel_out_ff,tunnel_out_fb]...
    =get_edgecount(G,datatable,data_fldr,temporal_fldr)
well_out_ff=[];
tunnel_out_ff=[];
well_out_fb=[];
tunnel_out_fb=[];
regions=["EC","DG","CA3","CA1"];
% tunnel_reg=["EC","DG","CA3","CA1"];
for fi=1:length(G)
    Edges_ff=G{fi}.Edges(G{fi}.Edges.Code==0.1,:);
    Edges_fb=G{fi}.Edges(G{fi}.Edges.Code==0.6,:);

    G_ff=digraph(Edges_ff,G{fi}.Nodes);
    G_fb=digraph(Edges_fb,G{fi}.Nodes);

    out_degree_ff=edgecount(G_ff);
    out_degree_fb=edgecount(G_fb);

    for regi=1:4
        mea=load(data_fldr+'\mea_'+regions(regi)+'.mat');
        if datatable.orientation(fi)=="CW"
            well_out_ff{fi,regi}=out_degree_ff(ismember(G_ff.Nodes.Name,...
                mea.cw.channel_names));
            well_out_fb{fi,regi}=out_degree_fb(ismember(G_fb.Nodes.Name,...
                mea.cw.channel_names));

            well_elec_ff{fi,regi}=G_ff.Nodes.Name(ismember(G_ff.Nodes.Name,...
                mea.cw.channel_names));
            well_elec_fb{fi,regi}=G_fb.Nodes.Name(ismember(G_fb.Nodes.Name,...
                mea.cw.channel_names));

            tunnels=load(data_fldr+'\cw_vars_tunnel.mat');
            load(temporal_fldr+'\tunnelChannelList.mat');
            ff_tunnel_List=circshift(tunnelChannelList.near_source.cw,1);
            tun_elec_ff=tunnels.channel_names(ff_tunnel_List{regi});
            fb_tunnel_List=circshift(tunnelChannelList.near_target.cw,0);
            tun_elec_fb=tunnels.channel_names(fb_tunnel_List{regi});

            tunnel_out_ff{fi,regi}=out_degree_ff(ismember(G_ff.Nodes.Name,...
                tun_elec_ff));
            tunnel_ff_elec{fi,regi}=tun_elec_ff;
            tunnel_ff_orientation{fi,regi}="CW";

            tunnel_out_fb{fi,regi}=out_degree_fb(ismember(G_fb.Nodes.Name,...
                tun_elec_fb));
            tunnel_fb_elec{fi,regi}=tun_elec_fb;
            tunnel_fb_orientation{fi,regi}="CW";
        else
            well_out_ff{fi,regi}=out_degree_ff(ismember(G_ff.Nodes.Name,...
                mea.ccw.channel_names));
            well_out_fb{fi,regi}=out_degree_fb(ismember(G_fb.Nodes.Name,...
                mea.ccw.channel_names));

            well_elec_ff{fi,regi}=G_ff.Nodes.Name(ismember(G_ff.Nodes.Name,...
                mea.ccw.channel_names));
            well_elec_fb{fi,regi}=G_fb.Nodes.Name(ismember(G_fb.Nodes.Name,...
                mea.ccw.channel_names));

            load(temporal_fldr+'\tunnelChannelList.mat');
            tunnels=load(data_fldr+'\cw_vars_tunnel.mat');
            load(temporal_fldr+'\tunnelChannelList.mat');
            ff_tunnel_List=circshift(tunnelChannelList.near_source.ccw,1);
            tun_elec_ff=tunnels.channel_names(ff_tunnel_List{regi});
            fb_tunnel_List=circshift(tunnelChannelList.near_target.ccw,0);
            tun_elec_fb=tunnels.channel_names(fb_tunnel_List{regi});

            tunnel_out_ff{fi,regi}=out_degree_ff(ismember(G_ff.Nodes.Name,...
                tun_elec_ff));
            tunnel_ff_elec{fi,regi}=tun_elec_ff;
            tunnel_ff_orientation{fi,regi}="CCW";

            tunnel_out_fb{fi,regi}=out_degree_fb(ismember(G_fb.Nodes.Name,...
                tun_elec_fb));
            tunnel_fb_elec{fi,regi}=tun_elec_fb;
            tunnel_fb_orientation{fi,regi}="CCW";
        end
    end
end

well_out_ff=cell2table([well_out_ff,well_elec_ff,tunnel_ff_orientation(:,1)],...
    "VariableNames",["EC","DG","CA3","CA1","EC Elec","DG Elec","CA3 Elec","CA1 Elec","Orientation"]);
well_out_fb=cell2table([well_out_fb,well_elec_fb,tunnel_fb_orientation(:,1)],...
    "VariableNames",["EC","DG","CA3","CA1","EC Elec","DG Elec","CA3 Elec","CA1 Elec","Orientation"]);
tunnel_out_ff=cell2table([tunnel_out_ff,tunnel_ff_elec,tunnel_ff_orientation(:,1)],...
    "VariableNames",["EC-DG","DG-CA3","CA3-CA1","CA1-EC",...
    "EC-DG Elec","DG-CA3 Elec","CA3-CA1 Elec","CA1-EC Elec","Orientation"]);
tunnel_out_fb=cell2table([tunnel_out_fb,tunnel_fb_elec,tunnel_fb_orientation(:,1)],...
    "VariableNames",["EC-DG","DG-CA3","CA3-CA1","CA1-EC",...
    "EC-DG Elec","DG-CA3 Elec","CA3-CA1 Elec","CA1-EC Elec","Orientation"]);
end

function [well_cent_ff,well_cent_fb,tunnel_cent_ff,tunnel_cent_fb]...
    =get_centrality_in(G,datatable,data_fldr,temporal_fldr)
well_cent_ff=[];
tunnel_cent_ff=[];
well_cent_fb=[];
tunnel_cent_fb=[];
regions=["EC","DG","CA3","CA1"];
% tunnel_reg=["EC","DG","CA3","CA1"];
for fi=1:length(G)
    Edges_ff=G{fi}.Edges(G{fi}.Edges.Code==0.1,:);
    Edges_fb=G{fi}.Edges(G{fi}.Edges.Code==0.6,:);

    G_ff=digraph(Edges_ff,G{fi}.Nodes);
    G_fb=digraph(Edges_fb,G{fi}.Nodes);

    out_cent_ff=centrality(G_ff,'indegree','Importance',G_ff.Edges.Weight);
    out_cent_fb=centrality(G_fb,'indegree','Importance',G_fb.Edges.Weight);

    for regi=1:4
        mea=load(data_fldr+'\mea_'+regions(regi)+'.mat');
        if datatable.orientation(fi)=="CW"
            well_cent_ff{fi,regi}=out_cent_ff(ismember(G_ff.Nodes.Name,...
                mea.cw.channel_names));
            well_cent_fb{fi,regi}=out_cent_fb(ismember(G_fb.Nodes.Name,...
                mea.cw.channel_names));

            well_elec_ff{fi,regi}=G_ff.Nodes.Name(ismember(G_ff.Nodes.Name,...
                mea.cw.channel_names));
            well_elec_fb{fi,regi}=G_fb.Nodes.Name(ismember(G_fb.Nodes.Name,...
                mea.cw.channel_names));

            tunnels=load(data_fldr+'\cw_vars_tunnel.mat');
            load(temporal_fldr+'\tunnelChannelList.mat');
            ff_tunnel_List=circshift(tunnelChannelList.near_source.cw,1);
            tun_elec_ff=tunnels.channel_names(ff_tunnel_List{regi});
            fb_tunnel_List=circshift(tunnelChannelList.near_target.cw,0);
            tun_elec_fb=tunnels.channel_names(fb_tunnel_List{regi});

            tunnel_cent_ff{fi,regi}=out_cent_ff(ismember(G_ff.Nodes.Name,...
                tun_elec_ff));
            tunnel_ff_elec{fi,regi}=tun_elec_ff;
            tunnel_ff_orientation{fi,regi}="CW";

            tunnel_cent_fb{fi,regi}=out_cent_fb(ismember(G_fb.Nodes.Name,...
                tun_elec_fb));
            tunnel_fb_elec{fi,regi}=tun_elec_fb;
            tunnel_fb_orientation{fi,regi}="CW";
        else
            well_cent_ff{fi,regi}=out_cent_ff(ismember(G_ff.Nodes.Name,...
                mea.ccw.channel_names));
            well_cent_fb{fi,regi}=out_cent_fb(ismember(G_fb.Nodes.Name,...
                mea.ccw.channel_names));

            well_elec_ff{fi,regi}=G_ff.Nodes.Name(ismember(G_ff.Nodes.Name,...
                mea.ccw.channel_names));
            well_elec_fb{fi,regi}=G_fb.Nodes.Name(ismember(G_fb.Nodes.Name,...
                mea.ccw.channel_names));

            load(temporal_fldr+'\tunnelChannelList.mat');
            tunnels=load(data_fldr+'\cw_vars_tunnel.mat');
            load(temporal_fldr+'\tunnelChannelList.mat');
            ff_tunnel_List=circshift(tunnelChannelList.near_source.ccw,1);
            tun_elec_ff=tunnels.channel_names(ff_tunnel_List{regi});
            fb_tunnel_List=circshift(tunnelChannelList.near_target.ccw,0);
            tun_elec_fb=tunnels.channel_names(fb_tunnel_List{regi});

            tunnel_cent_ff{fi,regi}=out_cent_ff(ismember(G_ff.Nodes.Name,...
                tun_elec_ff));
            tunnel_ff_elec{fi,regi}=tun_elec_ff;
            tunnel_ff_orientation{fi,regi}="CCW";

            tunnel_cent_fb{fi,regi}=out_cent_fb(ismember(G_fb.Nodes.Name,...
                tun_elec_fb));
            tunnel_fb_elec{fi,regi}=tun_elec_fb;
            tunnel_fb_orientation{fi,regi}="CCW";
        end
    end
end

well_cent_ff=cell2table([well_cent_ff,well_elec_ff,tunnel_ff_orientation(:,1)],...
    "VariableNames",["EC","DG","CA3","CA1","EC Elec","DG Elec","CA3 Elec","CA1 Elec","Orientation"]);
well_cent_fb=cell2table([well_cent_fb,well_elec_fb,tunnel_fb_orientation(:,1)],...
    "VariableNames",["EC","DG","CA3","CA1","EC Elec","DG Elec","CA3 Elec","CA1 Elec","Orientation"]);
tunnel_cent_ff=cell2table([tunnel_cent_ff,tunnel_ff_elec,tunnel_ff_orientation(:,1)],...
    "VariableNames",["EC-DG","DG-CA3","CA3-CA1","CA1-EC",...
    "EC-DG Elec","DG-CA3 Elec","CA3-CA1 Elec","CA1-EC Elec","Orientation"]);
tunnel_cent_fb=cell2table([tunnel_cent_fb,tunnel_fb_elec,tunnel_fb_orientation(:,1)],...
    "VariableNames",["EC-DG","DG-CA3","CA3-CA1","CA1-EC",...
    "EC-DG Elec","DG-CA3 Elec","CA3-CA1 Elec","CA1-EC Elec","Orientation"]);
end

function [well_cent_ff,well_cent_fb,tunnel_cent_ff,tunnel_cent_fb]...
    =get_centrality_out(G,datatable,data_fldr,temporal_fldr)
well_cent_ff=[];
tunnel_cent_ff=[];
well_cent_fb=[];
tunnel_cent_fb=[];
regions=["EC","DG","CA3","CA1"];
% tunnel_reg=["EC","DG","CA3","CA1"];
for fi=1:length(G)
    Edges_ff=G{fi}.Edges(G{fi}.Edges.Code==0.1,:);
    Edges_fb=G{fi}.Edges(G{fi}.Edges.Code==0.6,:);

    G_ff=digraph(Edges_ff,G{fi}.Nodes);
    G_fb=digraph(Edges_fb,G{fi}.Nodes);

    out_cent_ff=centrality(G_ff,'outdegree','Importance',G_ff.Edges.Weight);
    out_cent_fb=centrality(G_fb,'outdegree','Importance',G_fb.Edges.Weight);

    for regi=1:4
        mea=load(data_fldr+'\mea_'+regions(regi)+'.mat');
        if datatable.orientation(fi)=="CW"
            well_cent_ff{fi,regi}=out_cent_ff(ismember(G_ff.Nodes.Name,...
                mea.cw.channel_names));
            well_cent_fb{fi,regi}=out_cent_fb(ismember(G_fb.Nodes.Name,...
                mea.cw.channel_names));

            well_elec_ff{fi,regi}=G_ff.Nodes.Name(ismember(G_ff.Nodes.Name,...
                mea.cw.channel_names));
            well_elec_fb{fi,regi}=G_fb.Nodes.Name(ismember(G_fb.Nodes.Name,...
                mea.cw.channel_names));

            tunnels=load(data_fldr+'\cw_vars_tunnel.mat');
            load(temporal_fldr+'\tunnelChannelList.mat');
            ff_tunnel_List=circshift(tunnelChannelList.near_target.cw,1);
            tun_elec_ff=tunnels.channel_names(ff_tunnel_List{regi});
            fb_tunnel_List=circshift(tunnelChannelList.near_source.cw,0);
            tun_elec_fb=tunnels.channel_names(fb_tunnel_List{regi});

            tunnel_cent_ff{fi,regi}=out_cent_ff(ismember(G_ff.Nodes.Name,...
                tun_elec_ff));
            tunnel_ff_elec{fi,regi}=tun_elec_ff;
            tunnel_ff_orientation{fi,regi}="CW";

            tunnel_cent_fb{fi,regi}=out_cent_fb(ismember(G_fb.Nodes.Name,...
                tun_elec_fb));
            tunnel_fb_elec{fi,regi}=tun_elec_fb;
            tunnel_fb_orientation{fi,regi}="CW";
        else
            well_cent_ff{fi,regi}=out_cent_ff(ismember(G_ff.Nodes.Name,...
                mea.ccw.channel_names));
            well_cent_fb{fi,regi}=out_cent_fb(ismember(G_fb.Nodes.Name,...
                mea.ccw.channel_names));

            well_elec_ff{fi,regi}=G_ff.Nodes.Name(ismember(G_ff.Nodes.Name,...
                mea.ccw.channel_names));
            well_elec_fb{fi,regi}=G_fb.Nodes.Name(ismember(G_fb.Nodes.Name,...
                mea.ccw.channel_names));

            load(temporal_fldr+'\tunnelChannelList.mat');
            tunnels=load(data_fldr+'\cw_vars_tunnel.mat');
            load(temporal_fldr+'\tunnelChannelList.mat');
            ff_tunnel_List=circshift(tunnelChannelList.near_target.ccw,1);
            tun_elec_ff=tunnels.channel_names(ff_tunnel_List{regi});
            fb_tunnel_List=circshift(tunnelChannelList.near_source.ccw,0);
            tun_elec_fb=tunnels.channel_names(fb_tunnel_List{regi});

            tunnel_cent_ff{fi,regi}=out_cent_ff(ismember(G_ff.Nodes.Name,...
                tun_elec_ff));
            tunnel_ff_elec{fi,regi}=tun_elec_ff;
            tunnel_ff_orientation{fi,regi}="CCW";

            tunnel_cent_fb{fi,regi}=out_cent_fb(ismember(G_fb.Nodes.Name,...
                tun_elec_fb));
            tunnel_fb_elec{fi,regi}=tun_elec_fb;
            tunnel_fb_orientation{fi,regi}="CCW";
        end
    end
end

well_cent_ff=cell2table([well_cent_ff,well_elec_ff,tunnel_ff_orientation(:,1)],...
    "VariableNames",["EC","DG","CA3","CA1","EC Elec","DG Elec","CA3 Elec","CA1 Elec","Orientation"]);
well_cent_fb=cell2table([well_cent_fb,well_elec_fb,tunnel_fb_orientation(:,1)],...
    "VariableNames",["EC","DG","CA3","CA1","EC Elec","DG Elec","CA3 Elec","CA1 Elec","Orientation"]);
tunnel_cent_ff=cell2table([tunnel_cent_ff,tunnel_ff_elec,tunnel_ff_orientation(:,1)],...
    "VariableNames",["EC-DG","DG-CA3","CA3-CA1","CA1-EC",...
    "EC-DG Elec","DG-CA3 Elec","CA3-CA1 Elec","CA1-EC Elec","Orientation"]);
tunnel_cent_fb=cell2table([tunnel_cent_fb,tunnel_fb_elec,tunnel_fb_orientation(:,1)],...
    "VariableNames",["EC-DG","DG-CA3","CA3-CA1","CA1-EC",...
    "EC-DG Elec","DG-CA3 Elec","CA3-CA1 Elec","CA1-EC Elec","Orientation"]);
end

function [well_out_ff,well_out_fb,tunnel_out_ff,tunnel_out_fb]...
    =get_weights(G,datatable,data_fldr,temporal_fldr)
well_out_ff=[];
tunnel_out_ff=[];
well_out_fb=[];
tunnel_out_fb=[];
regions=["EC","DG","CA3","CA1"];
% tunnel_reg=["EC","DG","CA3","CA1"];
for fi=1:length(G)
    Edges_ff=G{fi}.Edges(G{fi}.Edges.Code==0.1,:);
    Edges_fb=G{fi}.Edges(G{fi}.Edges.Code==0.6,:);

    G_ff=digraph(Edges_ff,G{fi}.Nodes);
    G_fb=digraph(Edges_fb,G{fi}.Nodes);

    for regi=1:4
        mea=load(data_fldr+'\mea_'+regions(regi)+'.mat');
        if datatable.orientation(fi)=="CW"

            well_elec_ff{fi,regi}=G_ff.Nodes.Name(ismember(G_ff.Nodes.Name,...
                mea.cw.channel_names));
            well_elec_fb{fi,regi}=G_fb.Nodes.Name(ismember(G_fb.Nodes.Name,...
                mea.cw.channel_names));

            well_out_ff{fi,regi}=G_ff.Edges.Weight(ismember(G_ff.Edges.EndNodes(:,1),...
                mea.cw.channel_names));
            well_out_fb{fi,regi}=G_fb.Edges.Weight(ismember(G_fb.Edges.EndNodes(:,1),...
                mea.cw.channel_names));

            tunnels=load(data_fldr+'\cw_vars_tunnel.mat');
            load(temporal_fldr+'\tunnelChannelList.mat');
            ff_tunnel_List=circshift(tunnelChannelList.near_target.cw,0);
            tun_elec_ff=tunnels.channel_names(ff_tunnel_List{regi});
            fb_tunnel_List=circshift(tunnelChannelList.near_source.cw,1);
            tun_elec_fb=tunnels.channel_names(fb_tunnel_List{regi});

            tunnel_out_ff{fi,regi}=G_ff.Edges.Weight(ismember(G_ff.Edges.EndNodes(:,1),...
                tun_elec_ff));
            tunnel_ff_elec{fi,regi}=tun_elec_ff;
            tunnel_ff_orientation{fi,regi}="CW";

            tunnel_out_fb{fi,regi}=G_fb.Edges.Weight(ismember(G_fb.Edges.EndNodes(:,1),...
                tun_elec_fb));
            tunnel_fb_elec{fi,regi}=tun_elec_fb;
            tunnel_fb_orientation{fi,regi}="CW";
        else
            well_out_ff{fi,regi}=G_ff.Edges.Weight(ismember(G_ff.Edges.EndNodes(:,1),...
                mea.ccw.channel_names));
            well_out_fb{fi,regi}=G_fb.Edges.Weight(ismember(G_fb.Edges.EndNodes(:,1),...
                mea.ccw.channel_names));

            well_elec_ff{fi,regi}=G_ff.Nodes.Name(ismember(G_ff.Nodes.Name,...
                mea.ccw.channel_names));
            well_elec_fb{fi,regi}=G_fb.Nodes.Name(ismember(G_fb.Nodes.Name,...
                mea.ccw.channel_names));

            load(temporal_fldr+'\tunnelChannelList.mat');
            tunnels=load(data_fldr+'\cw_vars_tunnel.mat');
            load(temporal_fldr+'\tunnelChannelList.mat');
            ff_tunnel_List=circshift(tunnelChannelList.near_target.ccw,0);
            tun_elec_ff=tunnels.channel_names(ff_tunnel_List{regi});
            fb_tunnel_List=circshift(tunnelChannelList.near_source.ccw,1);
            tun_elec_fb=tunnels.channel_names(fb_tunnel_List{regi});

            tunnel_out_ff{fi,regi}=G_ff.Edges.Weight(ismember(G_ff.Edges.EndNodes(:,1),...
                tun_elec_ff));
            tunnel_ff_elec{fi,regi}=tun_elec_ff;
            tunnel_ff_orientation{fi,regi}="CCW";

            tunnel_out_fb{fi,regi}=G_fb.Edges.Weight(ismember(G_fb.Edges.EndNodes(:,1),...
                tun_elec_fb));
            tunnel_fb_elec{fi,regi}=tun_elec_fb;
            tunnel_fb_orientation{fi,regi}="CCW";
        end
    end
end

well_out_ff=cell2table([well_out_ff,well_elec_ff,tunnel_ff_orientation(:,1)],...
    "VariableNames",["EC","DG","CA3","CA1","EC Elec","DG Elec","CA3 Elec","CA1 Elec","Orientation"]);
well_out_fb=cell2table([well_out_fb,well_elec_fb,tunnel_fb_orientation(:,1)],...
    "VariableNames",["EC","DG","CA3","CA1","EC Elec","DG Elec","CA3 Elec","CA1 Elec","Orientation"]);
tunnel_out_ff=cell2table([tunnel_out_ff,tunnel_ff_elec,tunnel_ff_orientation(:,1)],...
    "VariableNames",["EC-DG","DG-CA3","CA3-CA1","CA1-EC",...
    "EC-DG Elec","DG-CA3 Elec","CA3-CA1 Elec","CA1-EC Elec","Orientation"]);
tunnel_out_fb=cell2table([tunnel_out_fb,tunnel_fb_elec,tunnel_fb_orientation(:,1)],...
    "VariableNames",["EC-DG","DG-CA3","CA3-CA1","CA1-EC",...
    "EC-DG Elec","DG-CA3 Elec","CA3-CA1 Elec","CA1-EC Elec","Orientation"]);
end

% function meanedcol=mean_col(table)
%     meanedcol=mean(vertcat(table{:}));
% end

function edges_per_subregion(table_cells,labels)
%ensure table cells have same number of variables in table

num_vars=width(table_cells{1});
figure
to_plot=[];
err=[];
for i=1:num_vars
    for j=1:length(table_cells)
        to_plot(i,j)=mean(table_cells{j}{:,i});
        err(i,j)=std(table_cells{j}{:,i})/sqrt(length(table_cells{j}{:,i}));
    end
end
hold on
b=bar(labels,to_plot);

% error bars
[ngroups,nbars]=size(to_plot);
x=nan(nbars,ngroups);
for j=1:nbars
    x(j,:)=b(j).XEndPoints;
end
errorbar(x',to_plot,err,'k','LineStyle','none')
hold off
end

function edges_per_subregion_array(table_cells,labels,array_num)
%ensure table cells have same number of variables in table

num_vars=width(table_cells{1});
figure
to_plot=[];
err=[];
for i=1:num_vars
    for j=1:length(table_cells)
        to_plot(i,j)=mean(table_cells{j}{array_num,i});
        err(i,j)=std(table_cells{j}{array_num,i})/sqrt(length(table_cells{j}{array_num,i}));
    end
end
hold on
b=bar(labels,to_plot);

% error bars
[ngroups,nbars]=size(to_plot);
x=nan(nbars,ngroups);
for j=1:nbars
    x(j,:)=b(j).XEndPoints;
end
errorbar(x',to_plot,err,'k','LineStyle','none')
hold off
end

function plot_horzbar(wt_ff_tab,wt_fb_tab,tw_ff_tab,tw_fb_tab,labels)

%get avg
num_vars=width(wt_ff_tab{1}); %ensure all tables have same num vars
figure

wt_ff_avg=[];
wtff_err=[];
for i=1:num_vars
    for j=1:length(wt_ff_tab)
        wt_ff_avg(i,j)=mean(wt_ff_tab{j}{:,i});
        wtff_err(i,j)=std(wt_ff_tab{j}{:,i})/sqrt(length(wt_ff_tab{j}{:,i}));
    end
end

wt_fb_avg=[];
wtfb_err=[];
for i=1:num_vars
    for j=1:length(wt_fb_tab)
        wt_fb_avg(i,j)=mean(wt_fb_tab{j}{:,i});
        wtfb_err(i,j)=std(wt_fb_tab{j}{:,i})/sqrt(length(wt_fb_tab{j}{:,i}));
    end
end

tw_ff_avg=[];
twff_err=[];
for i=1:num_vars
    for j=1:length(tw_ff_tab)
        tw_ff_avg(i,j)=mean(tw_ff_tab{j}{:,i});
        twff_err(i,j)=std(tw_ff_tab{j}{:,i})/sqrt(length(tw_ff_tab{j}{:,i}));
    end
end

tw_fb_avg=[];
twfb_err=[];
for i=1:num_vars
    for j=1:length(tw_fb_tab)
        tw_fb_avg(i,j)=mean(tw_fb_tab{j}{:,i});
        twfb_err(i,j)=std(tw_fb_tab{j}{:,i})/sqrt(length(tw_fb_tab{j}{:,i}));
    end
end

%make dat2plot
dat2plot_ff=[];
dat2plot_ff_err=[];
for i=1:4
    dat2plot_ff=[dat2plot_ff;wt_ff_avg(i,:);tw_ff_avg(i,:)];
    dat2plot_ff_err=[dat2plot_ff_err;wtff_err(i,:);twff_err(i,:)];
end

dat2plot_fb=[];
dat2plot_fb_err=[];
for i=1:4
    dat2plot_fb=[dat2plot_fb;wt_fb_avg(i,:);tw_fb_avg(i,:)];
    dat2plot_fb_err=[dat2plot_fb_err;wtfb_err(i,:);twfb_err(i,:)];
end

% labels_cat=categorical(labels);
% labels_cat=reordercats(labels_cat,flip(labels));
%
% bff=barh(labels_cat,dat2plot_ff);
% hold on
% bfb=barh(labels_cat,-dat2plot_fb);

% labels_cat=categorical(labels);
% labels_cat=reordercats(labels_cat,flip(labels));

dat2plot_ff=flip(dat2plot_ff,1);
dat2plot_ff_err=flip(dat2plot_ff_err,1);
dat2plot_fb=flip(dat2plot_fb,1);
dat2plot_fb_err=flip(dat2plot_fb_err,1);

% bff=barh(labels_cat,dat2plot_ff);
bff=barh(dat2plot_ff);
hold on
% bfb=barh(labels_cat,-dat2plot_fb);
bfb=barh(-dat2plot_fb);
yticklabels(flip(labels))

% error bars
[ngroups,nbars]=size(dat2plot_ff);
x=nan(nbars,ngroups);
for j=1:nbars
    x(j,:)=bff(j).XEndPoints;
end
errorbar(dat2plot_ff,x',dat2plot_ff_err,'horizontal','k','LineStyle','none')

[ngroups,nbars]=size(dat2plot_fb);
x=nan(nbars,ngroups);
for j=1:nbars
    x(j,:)=bfb(j).XEndPoints;
end
errorbar(-dat2plot_fb,x',dat2plot_fb_err,'horizontal','k','LineStyle','none')
hold off
end

function plot_horzbar_from_anova...
    (dat2plot_ff,dat2plot_ff_err,dat2plot_fb,dat2plot_fb_err,labels,ff_colors,fb_colors)

% labels_cat=categorical(labels);
% labels_cat=reordercats(labels_cat,flip(labels));

dat2plot_ff=flip(flip(dat2plot_ff,1),2);
dat2plot_ff_err=flip(flip(dat2plot_ff_err,1),2);
dat2plot_fb=flip(flip(dat2plot_fb,1),2);
dat2plot_fb_err=flip(flip(dat2plot_fb_err,1),2);

% bff=barh(labels_cat,dat2plot_ff);
bff=barh(dat2plot_ff);
change_horzbar_colors(bff,ff_colors)
temp_CDAT=bff(1).FaceColor;
bff(1).FaceColor=bff(3).FaceColor;
bff(3).FaceColor=temp_CDAT;
hold on
% bfb=barh(labels_cat,-dat2plot_fb);
bfb=barh(-dat2plot_fb);
change_horzbar_colors(bfb,fb_colors)
temp_CDAT=bfb(1).FaceColor;
bfb(1).FaceColor=bfb(3).FaceColor;
bfb(3).FaceColor=temp_CDAT;
yticklabels(flip(labels))

% error bars
[ngroups,nbars]=size(dat2plot_ff);
x=nan(nbars,ngroups);
for j=1:nbars
    x(j,:)=bff(j).XEndPoints;
end
errorbar(dat2plot_ff,x',dat2plot_ff_err,'horizontal','k','LineStyle','none')

[ngroups,nbars]=size(dat2plot_fb);
x=nan(nbars,ngroups);
for j=1:nbars
    x(j,:)=bfb(j).XEndPoints;
end
errorbar(-dat2plot_fb,x',dat2plot_fb_err,'horizontal','k','LineStyle','none')
hold off
end

function plot_horzbar_CohenD...
    (dat2plot,dat2plot_err_neg,dat2plot_err_pos,labels)

% labels_cat=categorical(labels);
% labels_cat=reordercats(labels_cat,flip(labels));

dat2plot=flip(flip(dat2plot,1),2);
dat2plot_err_neg=flip(flip(dat2plot_err_neg,1),2);
dat2plot_err_pos=flip(flip(dat2plot_err_pos,1),2);

% bff=barh(labels_cat,dat2plot_ff);
bff=barh(dat2plot);
temp_CDAT=bff(1).FaceColor;
bff(1).FaceColor=bff(3).FaceColor;
bff(3).FaceColor=temp_CDAT;
hold on
yticklabels(flip(labels))

% error bars
[ngroups,nbars]=size(dat2plot);
x=nan(nbars,ngroups);
for j=1:nbars
    x(j,:)=bff(j).XEndPoints;
end
errorbar(dat2plot,x',dat2plot_err_neg,dat2plot_err_pos,'horizontal','k','LineStyle','none')

hold off
end

function change_horzbar_colors(c_fig,colors)

for i=1:length(c_fig(:))

    c_fig(i).FaceColor=colors(i);

end

end

function newtab=exclude_FID(tab,datatable,fid_exclude)

row_exclude=datatable.s_no==fid_exclude;
tab(row_exclude,:)=[];
newtab=tab;


end

function effect=ff_fb_meanEffect(ff_1,fb_1,ff_2,fb_2,ff_3,fb_3)

effect{1}=meanEffectSize(ff_1,fb_1,'Effect','cohen');
effect{2}=meanEffectSize(ff_2,fb_2,'Effect','cohen');
effect{3}=meanEffectSize(ff_3,fb_3,'Effect','cohen');

effect=cell2table(effect,"VariableNames",["NoStim","5HFS","40HFS"]);
end