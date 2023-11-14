%% Create Graphs with strength of connections varying with time representing slopes of AFR between nodes
clear
clc
close all
%Loading prerequisites
% load("D:\Brewer lab data\HFS\Temporal Analysis\No Stim\dataInfo.mat")
% load("D:\Brewer lab data\HFS\Temporal Analysis\No Stim\full_idx_allregion_unit_matched_stim.mat")

% load("D:\Brewer lab data\HFS\Temporal Analysis\5 HFS\dataInfo.mat")
% load("D:\Brewer lab data\HFS\Temporal Analysis\5 HFS\full_idx_allregion_unit_matched_stim.mat")

load("D:\Brewer lab data\HFS\Temporal Analysis\40 HFS\dataInfo.mat")
load("D:\Brewer lab data\HFS\Temporal Analysis\40 HFS\full_idx_allregion_unit_matched_stim.mat")

allregion_unit_matched=convert_allregion_unit_matched_220222(allregion_unit_matched_stim);

% full spike index
% data_folder_addr="D:\Brewer lab data\HFS\No Stim\full_index_pseudo_times";
% data_folder_addr="D:\Brewer lab data\HFS\Theta Stim\full_index_pseudo_times";
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

%no stim
% source = load('D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output Full Index slide 10sLM\concatenated_source_table.mat');
% target = load('D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output Full Index slide 10sLM\concatenated_target_table.mat');

% 5 HFS
% source=load("D:\Brewer lab data\HFS\Temporal Analysis\5 HFS\AFR Output Full Index slide 10sLM\concatenated_source_table.mat");
% target=load("D:\Brewer lab data\HFS\Temporal Analysis\5 HFS\AFR Output Full Index slide 10sLM\concatenated_target_table.mat");

% 40 HFS
source=load("D:\Brewer lab data\HFS\Temporal Analysis\40 HFS\AFR Output Full Index slide 10sLM\concatenated_source_table.mat");
target=load("D:\Brewer lab data\HFS\Temporal Analysis\40 HFS\AFR Output Full Index slide 10sLM\concatenated_target_table.mat");

NodeTable = create_node_table_200820;
regList = ["EC","DG","CA3","CA1"];
source_regList_labels = ["EC->tunnels","DG->tunnels","CA3->tunnels","CA1->tunnels"];
target_regList_labels = ["tunnels->EC","tunnels->DG","tunnels->CA3","tunnels->CA1"];
% powerlawfitfun = @(b,x) 10.^(b(2)*log10(x) + b(1)); 
powerlawfitfunc= @(b,x) b(1) + b(2)*log(x);
%% Create Legend
figure('Position',[100,100,200,100])
x=[1 2 1 2 1 2 1 2 1 2];
y=[1 1 2 2 3 3 4 4 5 5];
weights_slope=[0.4:0.4:2];
weights_rsq=[0.2:0.2:1];
nodeTable_legend=table([1:length(x)]',x',y','VariableNames',{'No','X','Y'});
endNodes_legend=[1 2; 3 4; 5 6; 7 8; 9 10];
edgetable_slope=table(endNodes_legend,weights_slope','VariableNames',{'EndNodes','Weight'});

G_legend=graph(edgetable_slope,nodeTable_legend);
timedep_legend=plot(G_legend,'XData',x,'YData',y,'NodeLabel',...
    [],'NodeColor','none');

timedep_legend.LineWidth=(weights_slope*10)/2;
timedep_legend.EdgeColor=[1 0 0];
ylim([0.5,5.5])
set(gca,'Visible','off')
% set(gca,'XColor','none','YColor','none')
axis tight

% set(findall(gca,'type','text'),'visible','on')
saveas(gcf,"D:\Brewer lab data\HFS\Temporal Analysis\legend_widths.png")

%%
figure('Position',[100,100,200,100])
x=[1 2 1 2];
y=[1 1 1.2 1.2];
weights_slope=[2,2];
weights_rsq=[2,2];
nodeTable_legend=table([1:length(x)]',x',y','VariableNames',{'No','X','Y'});
endNodes_legend=[1 2; 3 4];
edgetable_slope=table(endNodes_legend,weights_slope','VariableNames',{'EndNodes','Weight'});

G_legend=graph(edgetable_slope,nodeTable_legend);
timedep_legend=plot(G_legend,'XData',x,'YData',y,'NodeLabel',...
    [],'NodeColor','none');

timedep_legend.LineWidth=(weights_slope*10)/2;
timedep_legend.EdgeColor=[0 0 1;1 0 0];
ylim([0.5,5.5])
set(gca,'Visible','off')
% set(gca,'XColor','none','YColor','none')
axis tight
saveas(gcf,"D:\Brewer lab data\HFS\Temporal Analysis\legend_direction.png")
%% Time dependent graphs plots
time_vec = [0.05:0.1:299.95];
no_ti = length(time_vec);
for graph_type = ["slope","rsq"]
    for fi=1:length(allregion_unit_matched)
        output_folder = "./time_depn_graph_plots/"+fi+"/"+graph_type+"/";
        if ~exist(output_folder,"dir"), mkdir(output_folder); end
        for ti = 1:no_ti
            source_table_full = source.AFR_table(source.AFR_table.fi == fi,:);
            target_table_full = target.AFR_table(target.AFR_table.fi == fi,:);
            
            % Compute time subset
            source_subset = time_subset(source_table_full, ti, 1e-4, 1e-4, 0);
            target_subset = time_subset(target_table_full, ti, 1e-4, 1e-4, 0);
            
            if isempty(source_subset) && isempty(target_subset)
                continue;
            end
            
            
            % Computing Edge Table with only positive slopes
            [EdgeTable] = create_edge_table_200820(NodeTable, target_subset, source_subset, graph_type);
            
            G = digraph(EdgeTable,NodeTable);
            % Plotting graphs
            figure('Renderer','painters','Position',[10 10 800 800])
            [gp,ax] = plot_graph(G,dataInfo.orientation{fi}, graph_type);
            plot_tunnels_in_graph_220309( allregion_unit_matched{fi},G,graph_type=="slope")
            axis square
            title("Time = "+time_vec(ti)+"s", "FontSize",20)
            saveas(gcf,char(output_folder+graph_type+"_graph_"+strrep(string(time_vec(ti)),".","p")),'png')
            
        end
        close all
    end
end

%% Time dependent graphs plots cumulative
time_vec = [0.05:0.1:299.95];
no_ti = length(time_vec);
for graph_type = ["slope","rsq"]
    for fi=1:length(allregion_unit_matched)
        output_folder = "./time_depn_graph_plots/"+fi+"/"+graph_type+"/";
        if ~exist(output_folder,"dir"), mkdir(output_folder); end
        figure('Renderer','painters','Position',[10 10 800 800])
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
            
            [gp,ax] = plot_graph(G_well_cat,dataInfo.orientation{fi}, graph_type);
            if strcmp(graph_type,"slope")
                plot_tunnels_in_graph_220309( allregion_unit_matched{fi},G_well_cat,1)
            else
                plot_tunnels_in_graph_220309( allregion_unit_matched{fi},G_well_cat,0)
            end
            %title(graph_type+" FID "+string(fi))
            set(ax,'position',[0 0 1 1])
            saveas(gcf,char(output_folder+graph_type+"_graph_full_time"),'png')
        end
        disp(string(fi))
    end
    close all
end
%% Time dependent graph video
time_vec = source.AFR_table.Properties.CustomProperties.mdl_par_time_vec;
no_ti = length(time_vec);
for fi=2
    output_folder = "./time_depn_graph_plots/"+fi+"/";
    if ~exist(output_folder,"dir"), mkdir(output_folder); end
    vid_file_name = char(output_folder+"vid_fi_"+fi);
    if exist(vid_file_name+".mp4","file"), delete(vid_file_name+".mp4"), end
    vidfile = VideoWriter(vid_file_name,'MPEG-4');
    open(vidfile);
    
    for ti = 1:no_ti
        source_table_full = source.AFR_table(source.AFR_table.fi == fi,:);
        target_table_full = target.AFR_table(target.AFR_table.fi == fi,:);
        
        % Compute time subset
        source_subset = time_subset(source_table_full, ti, 0.1, 0.2, 0);
        target_subset = time_subset(target_table_full, ti, 0.1, 0.2, 0);
        
        if isempty(source_subset) && isempty(target_subset)
            continue;
        end
        
        figure('Renderer','painters','Position',[10 10 1600 800])
        %%% plotting slopes
        % Computing Edge Table with only positive slopes
        subplot(1,2,1)
        graph_type = "slope";
        [EdgeTable] = create_edge_table_200820(NodeTable, target_subset, source_subset, graph_type);
        G = digraph(EdgeTable,NodeTable);
        % Plotting graphs
        [gp,~] = plot_graph(G,dataInfo.orientation{fi}, graph_type);
        plot_tunnels_in_graph_220309(allregion_unit_matched{fi},G,1)
        axis square
        title(["Slope","Time = "+time_vec(ti)+"s"], "FontSize",20)
        
        %%% plotting rsq with only positive slopes
        subplot(1,2,2)
        graph_type = "rsq";
        [EdgeTable] = create_edge_table_200820(NodeTable, target_subset, source_subset, graph_type);
        G = digraph(EdgeTable,NodeTable);
        % Plotting graphs
        [gp,ax] = plot_graph(G,dataInfo.orientation{fi}, graph_type);
        plot_tunnels_in_graph_220309(allregion_unit_matched{fi},G,0)
        axis square
        title(["R-sq","Time = "+time_vec(ti)+"s"], "FontSize",20)
        
        for i=1:20
            F = getframe(gcf);
            writeVideo(vidfile,F);
        end
        
    end
    close(vidfile)
    close all
end

%% Plot allregion no. of edges per channel bars
time_vec = source.AFR_table.Properties.CustomProperties.mdl_par_time_vec;
no_ti = length(time_vec);
for graph_type = ["source","target"]
    for fi=1
        output_folder = "./time_depn_graph_plots/"+fi+"/no_edge_distro/"+graph_type+"/";
        if ~exist(output_folder,"dir"), mkdir(output_folder); end
        for ti = 1:no_ti
            figure('Position',[100 100 1600 900]);
            for regi=1:4
                
                if graph_type == "source"
                    source_table_full = source.AFR_table(source.AFR_table.fi == fi & source.AFR_table.well_reg == regList(regi),: );
                    source_subset = time_subset(source_table_full, ti, 0.1, 0.2, 1);
                    if isempty(source_subset)
                        continue;
                    end
                    plot_edge_number_distro(source_subset, regi);
                    title(source_regList_labels(regi))
                else
                    target_table_full = target.AFR_table(target.AFR_table.fi == fi & target.AFR_table.well_reg == regList(regi),: );
                    target_subset = time_subset(target_table_full, ti, 0.1, 0.2, 1);
                    if isempty(target_subset)
                        continue;
                    end
                    plot_edge_number_distro(target_subset, regi);
                    title(target_regList_labels(regi))
                end
            end
            saveas(gcf,char(output_folder+time_vec(ti)),'png')
        end
    end
end

%% Plot time varying total no. of edges per subregion
time_vec = source.AFR_table.Properties.CustomProperties.mdl_par_time_vec;
no_ti = length(time_vec);
shift=[4,1,2,3];

for fi=1:length(allregion_unit_matched)
    figure('Position',[100 100 1600 350]);
    t=tiledlayout(4,2,'TileSpacing','tight','Padding','none');
    for graph_type = ["source","target"]
        output_folder = "./time_depn_graph_plots/"+fi+"/total_deges/";
        
        
        if ~exist(output_folder,"dir"), mkdir(output_folder); end
        for regi=1:4
            
            if graph_type=="source"
                loc=[1 3 5 7];
            else
                loc=[2 4 6 8];
            end
            
            vec2plot = zeros(2,no_ti);
            if graph_type=="source"
                %subplot(4,2,loc(regi))
                nexttile(loc(regi))
            else
                %subplot(4,2,loc(shift(regi)))
                nexttile(loc(shift(regi)))
            end
            ax=gca;
            ax.YScale="log";
            ylim([0.1,inf])
            hold on
            for ti = 1:no_ti
                if graph_type == "source"
                    source_table_full = source.AFR_table(source.AFR_table.fi == fi & source.AFR_table.well_reg == regList(regi),: );
                    subset = time_subset(source_table_full, ti, 0.05, 0, 1); % default 0.1 0.2 1
                else
                    target_table_full = target.AFR_table(target.AFR_table.fi == fi & target.AFR_table.well_reg == regList(regi),: );
                    subset = time_subset(target_table_full, ti, 0.05, 0.1, 1);
                end
                if ~isempty(subset)
                    vec2plot(1,ti)=sum(subset.if_ff);
                    vec2plot(2,ti)=sum(~subset.if_ff);
                end
            end
            vec2plot(vec2plot==0)=0.1;
            p=plot(time_vec, vec2plot(1,:),time_vec, vec2plot(2,:));
            p(1).Color = 'r';
            p(2).Color = 'b';
            %             xlim([200 220])
            if graph_type=="source"
                title(source_regList_labels(regi))
            else
                title(target_regList_labels(regi))
            end
            yticks(ax,[0, 1, 10, 20, 60, 100])
            yticklabels(ax,{0,1,10,20,50,100})
        end
        
        hold off
%         saveas(gcf,char(output_folder+datatype+"w200_220"),'png')
    end
    saveas(gcf,char(output_folder+"source_target_total_edges"),'png')
end

%% Plot r-squared histogram with thresh per fi

thresh=0.2;

for fi=1:max(source.AFR_table.fi)
    output_folder = "./time_depn_graph_plots/"+fi+"/rsq_dist/";
    if ~exist(output_folder,"dir"), mkdir(output_folder); end
    for regi=1:4
        figure
        t=tiledlayout(2,2);
        
        for graph_type=["source", "target"]
            for direction=["ff","fb"]
                
                
                bins=[thresh:0.05:1];
                %             bins=convert_edges_2_centers(bins);
                
                if graph_type=="source"
                    if direction=="ff"
                        rsq_vals=cell2mat(source.AFR_table.rsq);
                        rsq_vals(rsq_vals<thresh)=NaN;
                        rsq_vals=rsq_vals(logical(source.AFR_table.if_ff) & ...
                            logical(source.AFR_table.well_reg == regList(regi)) &...
                            logical(source.AFR_table.fi==fi),:);
                    else
                        rsq_vals=cell2mat(source.AFR_table.rsq);
                        rsq_vals(rsq_vals<thresh)=NaN;
                        rsq_vals=rsq_vals(~logical(source.AFR_table.if_ff)& ...
                            logical(source.AFR_table.well_reg == regList(regi)) &...
                            logical(source.AFR_table.fi==fi),:);
                    end
                else
                    if direction=="ff"
                        rsq_vals=cell2mat(target.AFR_table.rsq);
                        rsq_vals(rsq_vals<thresh)=NaN;
                        rsq_vals=rsq_vals(logical(target.AFR_table.if_ff)& ...
                            logical(target.AFR_table.well_reg == regList(regi)) &...
                            logical(target.AFR_table.fi==fi),:);
                    else
                        rsq_vals=cell2mat(target.AFR_table.rsq);
                        rsq_vals(rsq_vals<thresh)=NaN;
                        rsq_vals=rsq_vals(~logical(target.AFR_table.if_ff)& ...
                            logical(target.AFR_table.well_reg == regList(regi)) &...
                            logical(target.AFR_table.fi==fi),:);
                    end
                end
                
%                 figure
                nexttile
                hold on
                
                if direction=="ff"
                    h=histogram(rsq_vals(~isnan(rsq_vals)),bins,'FaceColor','r','Normalization','probability');
                else
                    h=histogram(rsq_vals(~isnan(rsq_vals)),bins,'FaceColor','b','Normalization','probability');
                end
                
                title(graph_type+" "+direction + " " + regList(regi))
                ylim([0,0.26])
                xlabel("R-squared")
                ylabel("Probability")
                
                cum_hist_sum=cumsum(h.Values);
                
                [~,fifty_perc_line]=min(abs(cum_hist_sum-0.5));
                
                xline(bins(fifty_perc_line+1))
                hold off
%                 saveas(gcf,output_folder+graph_type+"_"+direction+"_"+regList(regi)+"_"+"rsq_dist_thresh","png")
                
            end
        end
        title(t,"Rsq Dist. for "+regList(regi)+ " FID "+fi)
        saveas(gcf,output_folder+regList(regi)+"_"+"rsq_dist_thresh","png")

    end
end

%% Plot r-squared histogram with thresh concatonated pdf

thresh=0.2;

close all

for regi=1:4
%     figure('units','normalized','outerposition',[0 0 1 1])
    figure('Position',[500,50,800,800])
    t=tiledlayout(2,2);
    t.TileSpacing='compact';
    output_folder = "./time_depn_graph_plots/"+"all"+"/rsq_dist/";
    if ~exist(output_folder,"dir"), mkdir(output_folder); end
    
    for graph_type=["source", "target"]
        for direction=["ff","fb"]
            
%             output_folder = "./time_depn_graph_plots/"+"all"+"/rsq_dist/";
%             if ~exist(output_folder,"dir"), mkdir(output_folder); end
            bins=[thresh:0.05:1];
            bin_centers=convert_edges_2_centers(bins);
            %             bins=convert_edges_2_centers(bins);
            rsq_vals=[];
            hist_vals=[];
            hist_se=[];
            for fi=1:max(source.AFR_table.fi)
                if graph_type=="source"
                    if direction=="ff"
                        rsq_vals=cell2mat(source.AFR_table.rsq);
                        rsq_vals(rsq_vals<thresh)=NaN;
                        rsq_vals=rsq_vals(logical(source.AFR_table.if_ff) & ...
                            logical(source.AFR_table.well_reg == regList(regi)) &...
                            logical(source.AFR_table.fi==fi),:);
                    else
                        rsq_vals=cell2mat(source.AFR_table.rsq);
                        rsq_vals(rsq_vals<thresh)=NaN;
                        rsq_vals=rsq_vals(~logical(source.AFR_table.if_ff)& ...
                            logical(source.AFR_table.well_reg == regList(regi)) &...
                            logical(source.AFR_table.fi==fi),:);
                    end
                else
                    if direction=="ff"
                        rsq_vals=cell2mat(target.AFR_table.rsq);
                        rsq_vals(rsq_vals<thresh)=NaN;
                        rsq_vals=rsq_vals(logical(target.AFR_table.if_ff)& ...
                            logical(target.AFR_table.well_reg == regList(regi)) &...
                            logical(target.AFR_table.fi==fi),:);
                    else
                        rsq_vals=cell2mat(target.AFR_table.rsq);
                        rsq_vals(rsq_vals<thresh)=NaN;
                        rsq_vals=rsq_vals(~logical(target.AFR_table.if_ff)& ...
                            logical(target.AFR_table.well_reg == regList(regi)) &...
                            logical(target.AFR_table.fi==fi),:);
                    end
                end
                figure('units','normalized','outerposition',[0 0 1 1])
                h=histogram(rsq_vals(~isnan(rsq_vals)),bins,'Visible','off','Normalization','probability');
                hist_vals(fi,:)=h.Values;
                close(gcf)
            end
            
            hist_vals(~any(~isnan(hist_vals), 2),:)=[];
            hist_vec=mean(hist_vals,1); hist_se=stdErr(hist_vals);
            
            %                 figure
            nexttile
            hold on
            
            if direction=="ff"
                h=histogram('BinEdges',bins,'BinCounts',hist_vec,'FaceColor','r');
            else
                h=histogram('BinEdges',bins,'BinCounts',hist_vec,'FaceColor','b');
            end
            
            e=errorbar(bin_centers,hist_vec,hist_se,'k');
            e.LineStyle='none';
            
            title(graph_type+" "+direction + " " + regList(regi), 'FontSize', 16)
            subtitle("n= "+string(size(hist_vals,1)))
            ylim([0,0.2])
            xlabel("R-squared")
            ylabel("Probability")
            
            cum_hist_sum=cumsum(h.Values);
            
            [~,fifty_perc_line]=min(abs(cum_hist_sum-0.5));
            
            xline(bins(fifty_perc_line+1))
            
            ax=gca;
            ax.FontSize=16;
            
            hold off
%             saveas(gcf,output_folder+graph_type+"_"+direction+"_"+regList(regi)+"_"+"rsq_dist_thresh","png")
            
            
        end
    end
    title(t,"Rsq Dist. "+regList(regi), 'FontSize', 32)
    saveas(gcf,output_folder+graph_type+"_"+regList(regi)+"_"+"rsq_dist_thresh","png")
end
%% Plot r-squared histogram with thresh concatonated semilog pdf

thresh=0.2;

close all

for regi=1:4
    figure('units','normalized','outerposition',[0 0 1 1])
    t=tiledlayout(2,2);
    output_folder = "./time_depn_graph_plots/"+"all"+"/rsq_dist/";
    if ~exist(output_folder,"dir"), mkdir(output_folder); end
    
    for graph_type=["source", "target"]
        for direction=["ff","fb"]
            
%             output_folder = "./time_depn_graph_plots/"+"all"+"/rsq_dist/";
%             if ~exist(output_folder,"dir"), mkdir(output_folder); end
            bins=[thresh:0.05:1];
            bin_centers=convert_edges_2_centers(bins);
            %             bins=convert_edges_2_centers(bins);
            rsq_vals=[];
            hist_vals=[];
            hist_se=[];
            for fi=1:max(source.AFR_table.fi)
                if graph_type=="source"
                    if direction=="ff"
                        rsq_vals=cell2mat(source.AFR_table.rsq);
                        rsq_vals(rsq_vals<thresh)=NaN;
                        rsq_vals=rsq_vals(logical(source.AFR_table.if_ff) & ...
                            logical(source.AFR_table.well_reg == regList(regi)) &...
                            logical(source.AFR_table.fi==fi),:);
                    else
                        rsq_vals=cell2mat(source.AFR_table.rsq);
                        rsq_vals(rsq_vals<thresh)=NaN;
                        rsq_vals=rsq_vals(~logical(source.AFR_table.if_ff)& ...
                            logical(source.AFR_table.well_reg == regList(regi)) &...
                            logical(source.AFR_table.fi==fi),:);
                    end
                else
                    if direction=="ff"
                        rsq_vals=cell2mat(target.AFR_table.rsq);
                        rsq_vals(rsq_vals<thresh)=NaN;
                        rsq_vals=rsq_vals(logical(target.AFR_table.if_ff)& ...
                            logical(target.AFR_table.well_reg == regList(regi)) &...
                            logical(target.AFR_table.fi==fi),:);
                    else
                        rsq_vals=cell2mat(target.AFR_table.rsq);
                        rsq_vals(rsq_vals<thresh)=NaN;
                        rsq_vals=rsq_vals(~logical(target.AFR_table.if_ff)& ...
                            logical(target.AFR_table.well_reg == regList(regi)) &...
                            logical(target.AFR_table.fi==fi),:);
                    end
                end
                figure('units','normalized','outerposition',[0 0 1 1])
                h=histogram(rsq_vals(~isnan(rsq_vals)),bins,'Visible','off','Normalization','probability');
                hist_vals(fi,:)=h.Values;
                close(gcf)
            end
            
            hist_vals(~any(~isnan(hist_vals), 2),:)=[];
            hist_vec=mean(hist_vals,1); hist_se=stdErr(hist_vals);
            
            %                 figure
            nexttile
            hold on
            
            if direction=="ff"
                h=histogram('BinEdges',bins,'BinCounts',hist_vec,'FaceColor','r');
            else
                h=histogram('BinEdges',bins,'BinCounts',hist_vec,'FaceColor','b');
            end
            
            e=errorbar(bin_centers,hist_vec,hist_se,'k');
            e.LineStyle='none';
            
            title(graph_type+" "+direction + " " + regList(regi))
            subtitle("n= "+string(size(hist_vals,1)))
            ylim([0,0.5])
            xlabel("R-squared")
            ylabel("Probability")
            
            ax=gca;
            ax.XScale="log";
            
            cum_hist_sum=cumsum(h.Values);
            
            [~,fifty_perc_line]=min(abs(cum_hist_sum-0.5));
            
            xline(bins(fifty_perc_line+1))
            hold off
%             saveas(gcf,output_folder+graph_type+"_"+direction+"_"+regList(regi)+"_"+"rsq_dist_thresh","png")
            
            
        end
    end
    title(t,"Rsq Dist. "+regList(regi))
    saveas(gcf,output_folder+graph_type+"_"+regList(regi)+"_"+"rsq_dist_thresh_semilog","png")
end
%% Plot slope histogram with thresh per fi

thresh=0.1;

for fi=1:max(source.AFR_table.fi)
    output_folder = "./time_depn_graph_plots/"+fi+"/slope_dist/";
    if ~exist(output_folder,"dir"), mkdir(output_folder); end
    for regi=1:4
        figure
        t=tiledlayout(2,2);
        
        for graph_type=["source", "target"]
            for direction=["ff","fb"]
                
                
                bins=[thresh:0.1:2];
                %             bins=convert_edges_2_centers(bins);
                
                if graph_type=="source"
                    if direction=="ff"
                        slope_vals=cell2mat(source.AFR_table.slope);
                        slope_vals(slope_vals<thresh)=NaN;
                        slope_vals=slope_vals(logical(source.AFR_table.if_ff) & ...
                            logical(source.AFR_table.well_reg == regList(regi)) &...
                            logical(source.AFR_table.fi==fi),:);
                    else
                        slope_vals=cell2mat(source.AFR_table.slope);
                        slope_vals(slope_vals<thresh)=NaN;
                        slope_vals=slope_vals(~logical(source.AFR_table.if_ff)& ...
                            logical(source.AFR_table.well_reg == regList(regi)) &...
                            logical(source.AFR_table.fi==fi),:);
                    end
                else
                    if direction=="ff"
                        slope_vals=cell2mat(target.AFR_table.slope);
                        slope_vals(slope_vals<thresh)=NaN;
                        slope_vals=slope_vals(logical(target.AFR_table.if_ff)& ...
                            logical(target.AFR_table.well_reg == regList(regi)) &...
                            logical(target.AFR_table.fi==fi),:);
                    else
                        slope_vals=cell2mat(target.AFR_table.slope);
                        slope_vals(slope_vals<thresh)=NaN;
                        slope_vals=slope_vals(~logical(target.AFR_table.if_ff)& ...
                            logical(target.AFR_table.well_reg == regList(regi)) &...
                            logical(target.AFR_table.fi==fi),:);
                    end
                end
                
%                 figure
                nexttile
                hold on
                
                if direction=="ff"
                    h=histogram(slope_vals(~isnan(slope_vals)),bins,'FaceColor','r','Normalization','probability');
                else
                    h=histogram(slope_vals(~isnan(slope_vals)),bins,'FaceColor','b','Normalization','probability');
                end
                
                title(graph_type+" "+direction + " " + regList(regi))
                ylim([0,0.26])
                xlabel("Slope")
                ylabel("Probability")
                
                cum_hist_sum=cumsum(h.Values);
                
                [~,fifty_perc_line]=min(abs(cum_hist_sum-0.5));
                
                xline(bins(fifty_perc_line+1))
                hold off
%                 saveas(gcf,output_folder+graph_type+"_"+direction+"_"+regList(regi)+"_"+"rsq_dist_thresh","png")
                
            end
        end
        title(t,"Slope Dist. for "+regList(regi)+ " FID "+fi)
        saveas(gcf,output_folder+regList(regi)+"_"+"slope_dist_thresh","png")

    end
end

%% Plot slope histogram with thresh concatonated pdf

thresh=0.1;

close all

for regi=1:4
%     figure('units','normalized','outerposition',[0 0 1 1])
    figure('Position',[500,50,800,800])
    t=tiledlayout(2,2);
    t.TileSpacing='compact';
    output_folder = "./time_depn_graph_plots/"+"all"+"/slope_dist/";
    if ~exist(output_folder,"dir"), mkdir(output_folder); end
    
    for graph_type=["source", "target"]
        for direction=["ff","fb"]
            
%             output_folder = "./time_depn_graph_plots/"+"all"+"/rsq_dist/";
%             if ~exist(output_folder,"dir"), mkdir(output_folder); end
            bins=[thresh:0.1:3];
            bin_centers=convert_edges_2_centers(bins);
            %             bins=convert_edges_2_centers(bins);
            slope_vals=[];
            hist_vals=[];
            hist_se=[];
            for fi=1:max(source.AFR_table.fi)
                if graph_type=="source"
                    if direction=="ff"
                        slope_vals=cell2mat(source.AFR_table.slope);
                        slope_vals(slope_vals<thresh)=NaN;
                        slope_vals=slope_vals(logical(source.AFR_table.if_ff) & ...
                            logical(source.AFR_table.well_reg == regList(regi)) &...
                            logical(source.AFR_table.fi==fi),:);
                    else
                        slope_vals=cell2mat(source.AFR_table.slope);
                        slope_vals(slope_vals<thresh)=NaN;
                        slope_vals=slope_vals(~logical(source.AFR_table.if_ff)& ...
                            logical(source.AFR_table.well_reg == regList(regi)) &...
                            logical(source.AFR_table.fi==fi),:);
                    end
                else
                    if direction=="ff"
                        slope_vals=cell2mat(target.AFR_table.slope);
                        slope_vals(slope_vals<thresh)=NaN;
                        slope_vals=slope_vals(logical(target.AFR_table.if_ff)& ...
                            logical(target.AFR_table.well_reg == regList(regi)) &...
                            logical(target.AFR_table.fi==fi),:);
                    else
                        slope_vals=cell2mat(target.AFR_table.slope);
                        slope_vals(slope_vals<thresh)=NaN;
                        slope_vals=slope_vals(~logical(target.AFR_table.if_ff)& ...
                            logical(target.AFR_table.well_reg == regList(regi)) &...
                            logical(target.AFR_table.fi==fi),:);
                    end
                end
                figure('units','normalized','outerposition',[0 0 1 1])
                h=histogram(slope_vals(~isnan(slope_vals)),bins,'Visible','off','Normalization','probability');
                hist_vals(fi,:)=h.Values;
                close(gcf)
            end
            
            hist_vals(~any(~isnan(hist_vals), 2),:)=[];
            hist_vec=mean(hist_vals,1); hist_se=stdErr(hist_vals);
            
            %                 figure
            nexttile
            hold on
            
            if direction=="ff"
                h=histogram('BinEdges',bins,'BinCounts',hist_vec,'FaceColor','r');
            else
                h=histogram('BinEdges',bins,'BinCounts',hist_vec,'FaceColor','b');
            end
            
            e=errorbar(bin_centers,hist_vec,hist_se,'k');
            e.LineStyle='none';
            
            title(graph_type+" "+direction + " " + regList(regi), 'FontSize',16)
            subtitle("n= "+string(size(hist_vals,1)))
            ylim([0,0.2])
            xlim([thresh,3])
            xlabel("Slope")
            ylabel("Probability")
%             xticks([0:0.1:3])
            xticks([0.2:0.2:3])
            
            cum_hist_sum=cumsum(h.Values);
            
            [~,fifty_perc_line]=min(abs(cum_hist_sum-0.5));
            
            xline(bins(fifty_perc_line+1))
            
            ax=gca;
            ax.FontSize=16;
            
            hold off
%             saveas(gcf,output_folder+graph_type+"_"+direction+"_"+regList(regi)+"_"+"rsq_dist_thresh","png")
            
            
        end
    end
    title(t,"Slope Dist. "+regList(regi),'FontSize',32)
    saveas(gcf,output_folder+graph_type+"_"+regList(regi)+"_"+"slope_dist_thresh","png")
end
%% Plot slope histogram with thresh concatonated semilog pdf

thresh=0.1;

close all

for regi=1:4
%     figure('units','normalized','outerposition',[0 0 1 1])
    figure('Position',[500,50,800,800])
    t=tiledlayout(2,2);
    t.TileSpacing='compact';
    output_folder = "./time_depn_graph_plots/"+"all"+"/slope_dist/";
    if ~exist(output_folder,"dir"), mkdir(output_folder); end
    
    for graph_type=["source", "target"]
        for direction=["ff","fb"]
            
%             output_folder = "./time_depn_graph_plots/"+"all"+"/rsq_dist/";
%             if ~exist(output_folder,"dir"), mkdir(output_folder); end
%             bins=[thresh:0.1:3];
            bins=logspace(-1,1,20);
            bin_centers=convert_edges_2_centers(bins);
            %             bins=convert_edges_2_centers(bins);
            slope_vals=[];
            hist_vals=[];
            hist_se=[];
            coli=1;
            histcount=[];
            for fi=1:max(source.AFR_table.fi)
                if graph_type=="source"
                    if direction=="ff"
                        slope_vals=cell2mat(source.AFR_table.slope);
                        slope_vals(slope_vals<thresh)=NaN;
                        slope_vals=slope_vals(logical(source.AFR_table.if_ff) & ...
                            logical(source.AFR_table.well_reg == regList(regi)) &...
                            logical(source.AFR_table.fi==fi),:);
                    else
                        slope_vals=cell2mat(source.AFR_table.slope);
                        slope_vals(slope_vals<thresh)=NaN;
                        slope_vals=slope_vals(~logical(source.AFR_table.if_ff)& ...
                            logical(source.AFR_table.well_reg == regList(regi)) &...
                            logical(source.AFR_table.fi==fi),:);
                    end
                else
                    if direction=="ff"
                        slope_vals=cell2mat(target.AFR_table.slope);
                        slope_vals(slope_vals<thresh)=NaN;
                        slope_vals=slope_vals(logical(target.AFR_table.if_ff)& ...
                            logical(target.AFR_table.well_reg == regList(regi)) &...
                            logical(target.AFR_table.fi==fi),:);
                    else
                        slope_vals=cell2mat(target.AFR_table.slope);
                        slope_vals(slope_vals<thresh)=NaN;
                        slope_vals=slope_vals(~logical(target.AFR_table.if_ff)& ...
                            logical(target.AFR_table.well_reg == regList(regi)) &...
                            logical(target.AFR_table.fi==fi),:);
                    end
                end
                figure('units','normalized','outerposition',[0 0 1 1])
%                 h=histogram(slope_vals(~isnan(slope_vals)),bins,'Visible','off','Normalization','probability');
%                 histcount(fi,:)=h.Values;
                [hist_object,hist_plot,histcount(coli,:)] = log_binned_histogram( slope_vals(~isnan(slope_vals)), bins, "pdf", 0);
                coli=coli+1;
                close(gcf)
            end
            
            histcount(~any(~isnan(histcount), 2),:)=[];
            hist_vec=mean(histcount,1); hist_se=stdErr(histcount);
            
            %                 figure
            nexttile
            hold on
            
            if direction=="ff"
                h=histogram('BinEdges',bins,'BinCounts',hist_vec,'Visible','off');
                plot(bin_centers,hist_vec,'r')
            else
                h=histogram('BinEdges',bins,'BinCounts',hist_vec,'Visible','off');
                plot(bin_centers,hist_vec,'b')
            end
            
            e=errorbar(bin_centers,hist_vec,hist_se,'k');
            e.LineStyle='none';
            
            title(graph_type+" "+direction + " " + regList(regi), 'FontSize',16)
            subtitle("n= "+string(size(histcount,1)))
            ylim([0,0.2])
            xlim([-Inf,3])
            xlabel("Slope")
            ylabel("Probability")
            xticks([0.125,0.25, 0.5, 0.75,1.0,1.5, 2.0, 3.0])
            
            ax=gca;
            ax.XScale="log";
            ax.FontSize=16;
            
            cum_hist_sum=cumsum(h.Values);
            
            [~,fifty_perc_line]=min(abs(cum_hist_sum-0.5));
            
            xline(bins(fifty_perc_line+1))
            
            % fit line
            
            
            
            hold off
%             saveas(gcf,output_folder+graph_type+"_"+direction+"_"+regList(regi)+"_"+"rsq_dist_thresh","png")
            
            
        end
    end
    title(t,"Slope Dist. "+regList(regi),'FontSize',32)
    saveas(gcf,output_folder+graph_type+"_"+regList(regi)+"_"+"slope_dist_thresh_pdf","png")
end

%% Plot slope histogram with thresh concatonated semilog cdf

thresh=0.1;

close all

for regi=1:4
%     figure('units','normalized','outerposition',[0 0 1 1])
    figure('Position',[500,50,800,800])
    t=tiledlayout(2,2);
    t.TileSpacing='compact';
%     t.Padding='compact';
    output_folder = "./time_depn_graph_plots/"+"all"+"/slope_dist/";
    if ~exist(output_folder,"dir"), mkdir(output_folder); end
    
    for graph_type=["source", "target"]
        for direction=["ff","fb"]
            
%             output_folder = "./time_depn_graph_plots/"+"all"+"/rsq_dist/";
%             if ~exist(output_folder,"dir"), mkdir(output_folder); end
%             bins=[thresh:0.1:3];
            bins=logspace(-1,1,20);
            bin_centers=convert_edges_2_centers(bins);
            %             bins=convert_edges_2_centers(bins);
            slope_vals=[];
            hist_vals=[];
            hist_se=[];
            coli=1;
            histcount=[];
            %             for fi=1:max(source.AFR_table.fi)
            if graph_type=="source"
                if direction=="ff"
                    slope_vals=cell2mat(source.AFR_table.slope);
                    slope_vals(slope_vals<thresh)=NaN;
                    slope_vals=slope_vals(logical(source.AFR_table.if_ff) & ...
                        logical(source.AFR_table.well_reg == regList(regi)),:);
                else
                    slope_vals=cell2mat(source.AFR_table.slope);
                    slope_vals(slope_vals<thresh)=NaN;
                    slope_vals=slope_vals(~logical(source.AFR_table.if_ff)& ...
                        logical(source.AFR_table.well_reg == regList(regi)),:);
                end
            else
                if direction=="ff"
                    slope_vals=cell2mat(target.AFR_table.slope);
                    slope_vals(slope_vals<thresh)=NaN;
                    slope_vals=slope_vals(logical(target.AFR_table.if_ff)& ...
                        logical(target.AFR_table.well_reg == regList(regi)),:);
                else
                    slope_vals=cell2mat(target.AFR_table.slope);
                    slope_vals(slope_vals<thresh)=NaN;
                    slope_vals=slope_vals(~logical(target.AFR_table.if_ff)& ...
                        logical(target.AFR_table.well_reg == regList(regi)),:);
                end
            end
%             figure('units','normalized','outerposition',[0 0 1 1])
            %h=histogram(slope_vals(~isnan(slope_vals)),bins,'Visible','off','Normalization','probability');
            %histcount(fi,:)=h.Values;
            nexttile
            hold on
            [hist_object,hist_plot,prob, xbin] = log_binned_histogram( slope_vals(~isnan(slope_vals)), bins, "cdf");
            %coli=coli+1;
%             close(gcf)
            %             end
            
%             histcount(~any(~isnan(histcount), 2),:)=[];
%             hist_vec=mean(histcount,1); hist_se=stdErr(histcount);
            
            %                 figure
%             nexttile
%             hold on
            
            if direction=="ff"
%                 h=histogram('BinEdges',bins,'BinCounts',hist_vec,'Visible','off');
%                 plot(bin_centers,hist_vec,'r')
                hist_plot.Color='r';
            else
%                 h=histogram('BinEdges',bins,'BinCounts',hist_vec,'Visible','off');
%                 plot(bin_centers,hist_vec,'b')
                hist_plot.Color='b';
            end
            
%             e=errorbar(bin_centers,hist_vec,hist_se,'k');
%             e.LineStyle='none';
            
            title(graph_type+" "+direction + " " + regList(regi), 'FontSize', 16)
%             subtitle("n= "+string(size(histcount,1)))
%             ylim([0,0.3])
            xlim([-Inf,3])
            xlabel("Slope")
            ylabel("Probability")
%             xticks([0:0.1:3])
            xticks([0.125,0.25, 0.5, 0.75,1.0,1.5, 2.0, 3.0])
            
            ax=gca;
%             ax.YScale="log";
            ax.FontSize=16;
            
%             cum_hist_sum=cumsum(h.Values);
            
%             [~,fifty_perc_line]=min(abs(cum_hist_sum-0.5));
            
%             xline(bins(fifty_perc_line+1))
            
            % fit line
            
            
            
            hold off
%             saveas(gcf,output_folder+graph_type+"_"+direction+"_"+regList(regi)+"_"+"rsq_dist_thresh","png")
            
            
        end
    end
    title(t,"Slope Dist. "+regList(regi), 'FontSize', 32)
    saveas(gcf,output_folder+graph_type+"_"+regList(regi)+"_"+"slope_dist_thresh_cdf","png")
end
%% Plot scatter plots
time_vec = [0.05:0.1:299.95];
no_ti = length(time_vec);
for graph_type = ["source","target"]
    for fi=1:length(allregion_unit_matched)
        for if_ff=[0,1]
            output_folder = "./time_depn_graph_plots/"+fi+"/"+"scatter"+"/"+graph_type+"/";
            if ~exist(output_folder,"dir"), mkdir(output_folder); end
            for ti = 1:no_ti
                
            end
        end
        
    end
    
end
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
