%% 8 pt graph analysis average

%% Load parameters

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

ff_tunnels=["EC-DG","DG-CA3","CA3-CA1","CA1-EC"];
fb_tunnels=["DG-EC","CA3-DG","CA1-CA3","EC-CA1"];
tunnels=ff_tunnels;

% powerlawfitfun = @(b,x) 10.^(b(2)*log10(x) + b(1)); 
powerlawfitfunc= @(b,x) b(1) + b(2)*log(x);

colors=[0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;...
    0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.3010 0.7450 0.9330;0.6350 0.0780 0.1840];
colormap(colors)
close all
%% Load Graph Analysis

load("D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output Full Index slide 10sLM\nostim_slope_G.mat")
load("D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output Full Index slide 10sLM\nostim_rsq_G.mat")

load("D:\Brewer lab data\HFS\Temporal Analysis\5 HFS\AFR Output Full Index slide 10sLM\hfs5_slope_G.mat")
load("D:\Brewer lab data\HFS\Temporal Analysis\5 HFS\AFR Output Full Index slide 10sLM\hfs5_rsq_G.mat")

load("D:\Brewer lab data\HFS\Temporal Analysis\40 HFS\AFR Output Full Index slide 10sLM\hfs40_slope_G.mat")
load("D:\Brewer lab data\HFS\Temporal Analysis\40 HFS\AFR Output Full Index slide 10sLM\hfs40_rsq_G.mat")

nostim_data_fldr="D:\Brewer lab data\HFS\No Stim\full_index_pseudo_times";
nostim_temporal_fldr="D:\Brewer lab data\HFS\Temporal Analysis\No Stim";
hfs5_data_fldr="D:\Brewer lab data\HFS\Theta Stim\full_index_pseudo_times";
hfs5_temporal_fldr="D:\Brewer lab data\HFS\Temporal Analysis\5 HFS";
hfs40_data_fldr="D:\Brewer lab data\HFS\HFS Stim\full_index_pseudo_times";
hfs40_temporal_fldr="D:\Brewer lab data\HFS\Temporal Analysis\40 HFS";

%% Slope Averages

%get weights
[nostim_well_ff,nostim_well_fb,nostim_tunnel_ff,nostim_tunnel_fb]=...
    get_weights(nostim_slope_G,nostim_datatable,nostim_data_fldr,nostim_temporal_fldr);
[hfs5_well_ff,hfs5_well_fb,hfs5_tunnel_ff,hfs5_tunnel_fb]=...
    get_weights(hfs5_slope_G,hfs5_datatable,hfs5_data_fldr,hfs5_temporal_fldr);
[hfs40_well_ff,hfs40_well_fb,hfs40_tunnel_ff,hfs40_tunnel_fb]=...
    get_weights(hfs40_slope_G,hfs40_datatable,hfs40_data_fldr,hfs40_temporal_fldr);

% get per array average
nostim_well_ff_per_array=per_array_avg(nostim_well_ff,regList);
nostim_well_fb_per_array=per_array_avg(nostim_well_fb,regList);
nostim_tunnel_ff_per_array=per_array_avg(nostim_tunnel_ff,regList);
nostim_tunnel_fb_per_array=per_array_avg(nostim_tunnel_fb,regList);

hfs5_well_ff_per_array=per_array_avg(hfs5_well_ff,regList);
hfs5_well_fb_per_array=per_array_avg(hfs5_well_fb,regList);
hfs5_tunnel_ff_per_array=per_array_avg(hfs5_tunnel_ff,regList);
hfs5_tunnel_fb_per_array=per_array_avg(hfs5_tunnel_fb,regList);

hfs40_well_ff_per_array=per_array_avg(hfs40_well_ff,regList);
hfs40_well_fb_per_array=per_array_avg(hfs40_well_fb,regList);
hfs40_tunnel_ff_per_array=per_array_avg(hfs40_tunnel_ff,regList);
hfs40_tunnel_fb_per_array=per_array_avg(hfs40_tunnel_fb,regList);

% find total average
nostim_well_ff_total_avg=mean_reg(nostim_well_ff_per_array,regList,tunnels,1,1);
nostim_well_fb_total_avg=mean_reg(nostim_well_fb_per_array,regList,tunnels,1,0);
nostim_tunnel_ff_total_avg=mean_reg(nostim_tunnel_ff_per_array,regList,tunnels,0,1);
nostim_tunnel_fb_total_avg=mean_reg(nostim_tunnel_fb_per_array,regList,tunnels,0,0);

hfs5_well_ff_total_avg=mean_reg(hfs5_well_ff_per_array,regList,tunnels,1,1);
hfs5_well_fb_total_avg=mean_reg(hfs5_well_fb_per_array,regList,tunnels,1,0);
hfs5_tunnel_ff_total_avg=mean_reg(hfs5_tunnel_ff_per_array,regList,tunnels,0,1);
hfs5_tunnel_fb_total_avg=mean_reg(hfs5_tunnel_fb_per_array,regList,tunnels,0,0);

hfs40_well_ff_total_avg=mean_reg(hfs40_well_ff_per_array,regList,tunnels,1,1);
hfs40_well_fb_total_avg=mean_reg(hfs40_well_fb_per_array,regList,tunnels,1,0);
hfs40_tunnel_ff_total_avg=mean_reg(hfs40_tunnel_ff_per_array,regList,tunnels,0,1);
hfs40_tunnel_fb_total_avg=mean_reg(hfs40_tunnel_fb_per_array,regList,tunnels,0,0);

% ff/fb diff
perc_diff=@(f,i)(f-i)./(i).*100;

nostim_wt_diff_slope=perc_diff(nostim_well_ff_total_avg.avg_wt,nostim_tunnel_fb_total_avg.avg_wt);
nostim_tw_diff_slope=perc_diff(nostim_tunnel_ff_total_avg.avg_wt,nostim_well_fb_total_avg.avg_wt);

hfs5_wt_diff_slope=perc_diff(hfs5_well_ff_total_avg.avg_wt,hfs5_tunnel_fb_total_avg.avg_wt);
hfs5_tw_diff_slope=perc_diff(hfs5_tunnel_ff_total_avg.avg_wt,hfs5_well_fb_total_avg.avg_wt);

hfs40_wt_diff_slope=perc_diff(hfs40_well_ff_total_avg.avg_wt,hfs40_tunnel_fb_total_avg.avg_wt);
hfs40_tw_diff_slope=perc_diff(hfs40_tunnel_ff_total_avg.avg_wt,hfs40_well_fb_total_avg.avg_wt);

%% Slope Construct Graph Table

avg_G_nostim=...
    construct_G(nostim_well_ff_total_avg,nostim_well_fb_total_avg,nostim_tunnel_ff_total_avg,nostim_tunnel_fb_total_avg);
plot_8pt_graph(avg_G_nostim,"slope")
% saveas(gcf,"D:\Brewer lab data\HFS\three stim share figs\nostim_slope_8pt_avg_graph.png")
ax=gca;
exportgraphics(ax,"D:\Brewer lab data\HFS\three stim share figs\nostim_slope_8pt_avg_graph.png",'Resolution',1500)

avg_G_5hfs=...
    construct_G(hfs5_well_ff_total_avg,hfs5_well_fb_total_avg,hfs5_tunnel_ff_total_avg,hfs5_tunnel_fb_total_avg);
plot_8pt_graph(avg_G_5hfs,"slope")
% saveas(gcf,"D:\Brewer lab data\HFS\three stim share figs\hfs5_slope_8pt_avg_graph.png")
ax=gca;
exportgraphics(ax,"D:\Brewer lab data\HFS\three stim share figs\hfs5_slope_8pt_avg_graph.png",'Resolution',1500)

avg_G_40hfs=...
    construct_G(hfs40_well_ff_total_avg,hfs40_well_fb_total_avg,hfs40_tunnel_ff_total_avg,hfs40_tunnel_fb_total_avg);
plot_8pt_graph(avg_G_40hfs,"slope")
% saveas(gcf,"D:\Brewer lab data\HFS\three stim share figs\hfs40_slope_8pt_avg_graph.png")
ax=gca;
exportgraphics(ax,"D:\Brewer lab data\HFS\three stim share figs\hfs40_slope_8pt_avg_graph.png",'Resolution',1500)
%% Anova Slope

%well ff
[well_ff_c,well_ff_m,well_ff_se]=calc_anova(nostim_well_ff,hfs5_well_ff,hfs40_well_ff,regList);

%well fb
[well_fb_c,well_fb_m,well_fb_se]=calc_anova(nostim_well_fb,hfs5_well_fb,hfs40_well_fb,regList);

%tunnel ff
[tunnel_ff_c,tunnel_ff_m,tunnel_ff_se]=calc_anova(nostim_tunnel_ff,hfs5_tunnel_ff,hfs40_tunnel_ff,regList);

%tunnel fb
[tunnel_fb_c,tunnel_fb_m,tunnel_fb_se]=calc_anova(nostim_tunnel_fb,hfs5_tunnel_fb,hfs40_tunnel_fb,regList);

%plot
close all
figure('units','normalized','outerposition',[0 0 1 1])
t=tiledlayout(1,4,'TileSpacing','compact');

means={well_ff_m,well_fb_m,tunnel_ff_m,tunnel_fb_m};
err={well_ff_se,well_fb_se,tunnel_ff_se,tunnel_fb_se};

for i=1:4
    ax(i)=nexttile;
    hold on
    regs=["EC","DG","CA3","CA1"];
    regs_cat=categorical(regs);
    regs_cat=reordercats(regs_cat,regs);
    regs_tunnel_ff=["EC-DG","DG-CA3","CA3-CA1","CA1-EC"];
    regs_tunnel_fb=["DG-EC","CA3-DG","CA1-CA3","EC-CA1"];
    
    regs_tunnel_ff_cat=reordercats(categorical(regs_tunnel_ff),regs_tunnel_ff);
    regs_tunnel_fb_cat=reordercats(categorical(regs_tunnel_fb),regs_tunnel_fb);
    
    if i==3
        title("FF Tunnel to Well")
        b=bar(regs_tunnel_ff_cat,means{i});
    elseif i==4
        title("FB Tunnel to Well")
        b=bar(regs_tunnel_fb_cat,means{i});
    elseif i==1
        title("FF Well to Tunnel")
        b=bar(regs_cat,means{i});
    elseif i==2
        title("FB Well to Tunnel")
        b=bar(regs_cat,means{i});
    end
    ax=gca;
    ax.FontSize=24;
    % error bars
    [ngroups,nbars]=size(means{i});
    x=nan(nbars,ngroups);
    for j=1:nbars
        x(j,:)=b(j).XEndPoints;
    end
    errorbar(x',means{i},err{i},'k','LineStyle','none')
    ylim([0.4,1.25])
    hold off
end
linkaxes(ax(:),'y')

saveas(gcf,"D:\Brewer lab data\HFS\three stim share figs\slope_bar.png")

%% Per array bar slope

for fid=1:10
    % graphing degrees averaged
    node_degrees_ff_well=[];
    node_degrees_ff_well_se=[];
    node_degrees_fb_well=[];
    node_degrees_fb_well_se=[];
    node_degrees_ff_tunnel=[];
    node_degrees_ff_tunnel_se=[];
    node_degrees_fb_tunnel=[];
    node_degrees_fb_tunnel_se=[];
    
    %relative idx
    nostim_idx=find(nostim_datatable.s_no==fid);
    hfs5_idx=find(hfs5_datatable.s_no==fid);
    hfs40_idx=find(hfs40_datatable.s_no==fid);
    
    %well in/out
    for regi=1:4
        
        %out
        %ff wells out
        nostim_well_ff_degs=nonzeros(cell2mat(nostim_well_ff{nostim_idx,regi}));
        hfs5_well_ff_degs=nonzeros(cell2mat(hfs5_well_ff{hfs5_idx,regi}));
        hfs40_well_ff_degs=nonzeros(cell2mat(hfs40_well_ff{hfs40_idx,regi}));
        node_degrees_ff_well=[node_degrees_ff_well;...
            [mean(nostim_well_ff_degs),mean(hfs5_well_ff_degs),mean(hfs40_well_ff_degs)]];
        node_degrees_ff_well_se=[node_degrees_ff_well_se;...
            [std(nostim_well_ff_degs)/sqrt(length(nostim_well_ff_degs)),...
            std(hfs5_well_ff_degs)/sqrt(length(hfs5_well_ff_degs)),...
            std(hfs40_well_ff_degs)/sqrt(length(hfs40_well_ff_degs))]];
        
        %fb wells out
        nostim_well_fb_degs=nonzeros(cell2mat(nostim_well_fb{nostim_idx,regi}));
        hfs5_well_fb_degs=nonzeros(cell2mat(hfs5_well_fb{hfs5_idx,regi}));
        hfs40_well_fb_degs=nonzeros(cell2mat(hfs40_well_fb{hfs40_idx,regi}));
        node_degrees_fb_well=[node_degrees_fb_well;...
            [mean(nostim_well_fb_degs),mean(hfs5_well_fb_degs),mean(hfs40_well_fb_degs)]];
        node_degrees_fb_well_se=[node_degrees_fb_well_se;...
            [std(nostim_well_fb_degs)/sqrt(length(nostim_well_fb_degs)),...
            std(hfs5_well_fb_degs)/sqrt(length(hfs5_well_fb_degs)),...
            std(hfs40_well_fb_degs)/sqrt(length(hfs40_well_fb_degs))]];
        
        %in
        %ff in
        nostim_tunnel_ff_degs=nonzeros(cell2mat(nostim_tunnel_ff{nostim_idx,regi}));
        hfs5_tunnel_ff_degs=nonzeros(cell2mat(hfs5_tunnel_ff{hfs5_idx,regi}));
        hfs40_tunnel_ff_degs=nonzeros(cell2mat(hfs40_tunnel_ff{hfs40_idx,regi}));
        node_degrees_ff_tunnel=[node_degrees_ff_tunnel;...
            [mean(nostim_tunnel_ff_degs),mean(hfs5_tunnel_ff_degs),mean(hfs40_tunnel_ff_degs)]];
        node_degrees_ff_tunnel_se=[node_degrees_ff_tunnel_se;...
            [std(nostim_tunnel_ff_degs)/sqrt(length(nostim_tunnel_ff_degs)),...
            std(hfs5_tunnel_ff_degs)/sqrt(length(hfs5_tunnel_ff_degs)),...
            std(hfs40_tunnel_ff_degs)/sqrt(length(hfs40_tunnel_ff_degs))]];
        
        %fb in
        nostim_tunnel_fb_degs=nonzeros(cell2mat(nostim_tunnel_fb{nostim_idx,regi}));
        hfs5_tunnel_fb_degs=nonzeros(cell2mat(hfs5_tunnel_fb{hfs5_idx,regi}));
        hfs40_tunnel_fb_degs=nonzeros(cell2mat(hfs40_tunnel_fb{hfs40_idx,regi}));
        node_degrees_fb_tunnel=[node_degrees_fb_tunnel;...
            [mean(nostim_tunnel_fb_degs),mean(hfs5_tunnel_fb_degs),mean(hfs40_tunnel_fb_degs)]];
        node_degrees_fb_tunnel_se=[node_degrees_fb_tunnel_se;...
            [std(nostim_tunnel_fb_degs)/sqrt(length(nostim_tunnel_fb_degs)),...
            std(hfs5_tunnel_fb_degs)/sqrt(length(hfs5_tunnel_fb_degs)),...
            std(hfs40_tunnel_fb_degs)/sqrt(length(hfs40_tunnel_fb_degs))]];
        
        % ff in anova
        name_mat=[repmat("No Stim",length(nostim_well_ff_degs),1);repmat("5 HFS",length(hfs5_well_ff_degs),1);...
            repmat("40 HFS",length(hfs40_well_ff_degs),1)];
%         [~,~,stats]=anova1([nostim_well_ff_degs;hfs5_well_ff_degs;hfs40_well_ff_degs],name_mat);
%         [c,~]=multcompare(stats);
%         ff_in_well_c{regi}=c;
        
        % fb in anova
        name_mat=[repmat("No Stim",length(nostim_well_fb_degs),1);repmat("5 HFS",length(hfs5_well_fb_degs),1);...
            repmat("40 HFS",length(hfs40_well_fb_degs),1)];
%         [~,~,stats]=anova1([nostim_well_fb_degs;hfs5_well_fb_degs;hfs40_well_fb_degs],name_mat);
%         [c,~]=multcompare(stats);
%         fb_in_well_c{regi}=c;
        
        % ff out anova
        name_mat=[repmat("No Stim",length(nostim_tunnel_ff_degs),1);repmat("5 HFS",length(hfs5_tunnel_ff_degs),1);...
            repmat("40 HFS",length(hfs40_tunnel_ff_degs),1)];
%         [~,~,stats]=anova1([nostim_tunnel_ff_degs;hfs5_tunnel_ff_degs;hfs40_tunnel_ff_degs],name_mat);
%         [c,~]=multcompare(stats);
%         ff_out_well_c{regi}=c;
        
        % fb in anova
        name_mat=[repmat("No Stim",length(nostim_tunnel_fb_degs),1);repmat("5 HFS",length(hfs5_tunnel_fb_degs),1);...
            repmat("40 HFS",length(hfs40_tunnel_fb_degs),1)];
%         [~,~,stats]=anova1([nostim_tunnel_fb_degs;hfs5_tunnel_fb_degs;hfs40_tunnel_fb_degs],name_mat);
%         [c,~]=multcompare(stats);
%         fb_out_well_c{regi}=c;
    end
    
    degs={node_degrees_ff_well,node_degrees_fb_well,node_degrees_ff_tunnel,node_degrees_fb_tunnel};
    degs_se={node_degrees_ff_well_se,node_degrees_fb_well_se,node_degrees_ff_tunnel_se,node_degrees_fb_tunnel_se};
    
    figure('units','normalized','outerposition',[0 0 1 1])
    t=tiledlayout(1,4,'TileSpacing','compact');
    ylabel(t,"Degrees Per Node",'FontSize',24)
    xlabel(t,"Subregion",'FontSize',24)
    title(t,"Slope Weights FID "+string(fid),'FontSize',24)
    for i=1:4
        ax(i)=nexttile;
        hold on
        regs=["EC","DG","CA3","CA1"];
        regs_cat=categorical(regs);
        regs_cat=reordercats(regs_cat,regs);
        regs_tunnel_ff=["EC-DG","DG-CA3","CA3-CA1","CA1-EC"];
        regs_tunnel_fb=["DG-EC","CA3-DG","CA1-CA3","EC-CA1"];
        
        regs_tunnel_ff_cat=reordercats(categorical(regs_tunnel_ff),regs_tunnel_ff);
        regs_tunnel_fb_cat=reordercats(categorical(regs_tunnel_fb),regs_tunnel_fb);
        
        if i==3
            title("Feedforward Tunnel",'FontSize',20)
            b=bar(regs_cat,degs{i});
        elseif i==4
            title("Feedback Tunnel",'FontSize',20)
            b=bar(regs_cat,degs{i});
        elseif i==1
            title("Feedforward Well",'FontSize',20)
            b=bar(regs_cat,degs{i});
        elseif i==2
            title("Feedback Well",'FontSize',20)
            b=bar(regs_cat,degs{i});
        end
        % error bars
        [ngroups,nbars]=size(degs{i});
        x=nan(nbars,ngroups);
        for j=1:nbars
            x(j,:)=b(j).XEndPoints;
        end
        errorbar(x',degs{i},degs_se{i},'k','LineStyle','none')
        
        hold off
    end
    linkaxes(ax(:),'y')
    colororder(colors)
    saveas(gcf,"D:\Brewer lab data\HFS\Temporal Analysis\Graph Analysis\slope weights nonzero FID "+string(fid)+".png")
end

%% rsq Averages

%get weights
[nostim_well_ff,nostim_well_fb,nostim_tunnel_ff,nostim_tunnel_fb]=...
    get_weights(nostim_rsq_G,nostim_datatable,nostim_data_fldr,nostim_temporal_fldr);
[hfs5_well_ff,hfs5_well_fb,hfs5_tunnel_ff,hfs5_tunnel_fb]=...
    get_weights(hfs5_rsq_G,hfs5_datatable,hfs5_data_fldr,hfs5_temporal_fldr);
[hfs40_well_ff,hfs40_well_fb,hfs40_tunnel_ff,hfs40_tunnel_fb]=...
    get_weights(hfs40_rsq_G,hfs40_datatable,hfs40_data_fldr,hfs40_temporal_fldr);

% get per array average
nostim_well_ff_per_array=per_array_avg(nostim_well_ff,regList);
nostim_well_fb_per_array=per_array_avg(nostim_well_fb,regList);
nostim_tunnel_ff_per_array=per_array_avg(nostim_tunnel_ff,regList);
nostim_tunnel_fb_per_array=per_array_avg(nostim_tunnel_fb,regList);

hfs5_well_ff_per_array=per_array_avg(hfs5_well_ff,regList);
hfs5_well_fb_per_array=per_array_avg(hfs5_well_fb,regList);
hfs5_tunnel_ff_per_array=per_array_avg(hfs5_tunnel_ff,regList);
hfs5_tunnel_fb_per_array=per_array_avg(hfs5_tunnel_fb,regList);

hfs40_well_ff_per_array=per_array_avg(hfs40_well_ff,regList);
hfs40_well_fb_per_array=per_array_avg(hfs40_well_fb,regList);
hfs40_tunnel_ff_per_array=per_array_avg(hfs40_tunnel_ff,regList);
hfs40_tunnel_fb_per_array=per_array_avg(hfs40_tunnel_fb,regList);

% find total average
nostim_well_ff_total_avg=mean_reg(nostim_well_ff_per_array,regList,tunnels,1,1);
nostim_well_fb_total_avg=mean_reg(nostim_well_fb_per_array,regList,tunnels,1,0);
nostim_tunnel_ff_total_avg=mean_reg(nostim_tunnel_ff_per_array,regList,tunnels,0,1);
nostim_tunnel_fb_total_avg=mean_reg(nostim_tunnel_fb_per_array,regList,tunnels,0,0);

hfs5_well_ff_total_avg=mean_reg(hfs5_well_ff_per_array,regList,tunnels,1,1);
hfs5_well_fb_total_avg=mean_reg(hfs5_well_fb_per_array,regList,tunnels,1,0);
hfs5_tunnel_ff_total_avg=mean_reg(hfs5_tunnel_ff_per_array,regList,tunnels,0,1);
hfs5_tunnel_fb_total_avg=mean_reg(hfs5_tunnel_fb_per_array,regList,tunnels,0,0);

hfs40_well_ff_total_avg=mean_reg(hfs40_well_ff_per_array,regList,tunnels,1,1);
hfs40_well_fb_total_avg=mean_reg(hfs40_well_fb_per_array,regList,tunnels,1,0);
hfs40_tunnel_ff_total_avg=mean_reg(hfs40_tunnel_ff_per_array,regList,tunnels,0,1);
hfs40_tunnel_fb_total_avg=mean_reg(hfs40_tunnel_fb_per_array,regList,tunnels,0,0);

% ff/fb diff
perc_diff=@(f,i)(f-i)./(i).*100;

nostim_wt_diff_rsq=perc_diff(nostim_well_ff_total_avg.avg_wt,nostim_tunnel_fb_total_avg.avg_wt);
nostim_tw_diff_rsq=perc_diff(nostim_tunnel_ff_total_avg.avg_wt,nostim_well_fb_total_avg.avg_wt);

hfs5_wt_diff_rsq=perc_diff(hfs5_well_ff_total_avg.avg_wt,hfs5_tunnel_fb_total_avg.avg_wt);
hfs5_tw_diff_rsq=perc_diff(hfs5_tunnel_ff_total_avg.avg_wt,hfs5_well_fb_total_avg.avg_wt);

hfs40_wt_diff_rsq=perc_diff(hfs40_well_ff_total_avg.avg_wt,hfs40_tunnel_fb_total_avg.avg_wt);
hfs40_tw_diff_rsq=perc_diff(hfs40_tunnel_ff_total_avg.avg_wt,hfs40_well_fb_total_avg.avg_wt);
%% rsq Construct Graph Table

avg_G_nostim=...
    construct_G(nostim_well_ff_total_avg,nostim_well_fb_total_avg,nostim_tunnel_ff_total_avg,nostim_tunnel_fb_total_avg);
plot_8pt_graph(avg_G_nostim,"rsq")
% saveas(gcf,"D:\Brewer lab data\HFS\three stim share figs\nostim_rsq_8pt_avg_graph.png")
ax=gca;
exportgraphics(ax,"D:\Brewer lab data\HFS\three stim share figs\nostim_rsq_8pt_avg_graph.png",'Resolution',1500)

avg_G_5hfs=...
    construct_G(hfs5_well_ff_total_avg,hfs5_well_fb_total_avg,hfs5_tunnel_ff_total_avg,hfs5_tunnel_fb_total_avg);
plot_8pt_graph(avg_G_5hfs,"rsq")
% saveas(gcf,"D:\Brewer lab data\HFS\three stim share figs\hfs5_rsq_8pt_avg_graph.png")
ax=gca;
exportgraphics(ax,"D:\Brewer lab data\HFS\three stim share figs\hfs5_rsq_8pt_avg_graph.png",'Resolution',1500)

avg_G_40hfs=...
    construct_G(hfs40_well_ff_total_avg,hfs40_well_fb_total_avg,hfs40_tunnel_ff_total_avg,hfs40_tunnel_fb_total_avg);
plot_8pt_graph(avg_G_40hfs,"rsq")
% saveas(gcf,"D:\Brewer lab data\HFS\three stim share figs\hfs40_rsq_8pt_avg_graph.png")
ax=gca;
exportgraphics(ax,"D:\Brewer lab data\HFS\three stim share figs\hfs40_rsq_8pt_avg_graph.png",'Resolution',1500)

%% Anova rsq

%well ff
[well_ff_c,well_ff_m,well_ff_se]=calc_anova(nostim_well_ff,hfs5_well_ff,hfs40_well_ff,regList);

%well fb
[well_fb_c,well_fb_m,well_fb_se]=calc_anova(nostim_well_fb,hfs5_well_fb,hfs40_well_fb,regList);

%tunnel ff
[tunnel_ff_c,tunnel_ff_m,tunnel_ff_se]=calc_anova(nostim_tunnel_ff,hfs5_tunnel_ff,hfs40_tunnel_ff,regList);

%tunnel fb
[tunnel_fb_c,tunnel_fb_m,tunnel_fb_se]=calc_anova(nostim_tunnel_fb,hfs5_tunnel_fb,hfs40_tunnel_fb,regList);

%plot
close all
figure('units','normalized','outerposition',[0 0 1 1])
t=tiledlayout(1,4,'TileSpacing','compact');

means={well_ff_m,well_fb_m,tunnel_ff_m,tunnel_fb_m};
err={well_ff_se,well_fb_se,tunnel_ff_se,tunnel_fb_se};

for i=1:4
    ax(i)=nexttile;
    hold on
    regs=["EC","DG","CA3","CA1"];
    regs_cat=categorical(regs);
    regs_cat=reordercats(regs_cat,regs);
    regs_tunnel_ff=["EC-DG","DG-CA3","CA3-CA1","CA1-EC"];
    regs_tunnel_fb=["DG-EC","CA3-DG","CA1-CA3","EC-CA1"];
    
    regs_tunnel_ff_cat=reordercats(categorical(regs_tunnel_ff),regs_tunnel_ff);
    regs_tunnel_fb_cat=reordercats(categorical(regs_tunnel_fb),regs_tunnel_fb);
    
    if i==3
        title("FF Tunnel to Well")
        b=bar(regs_tunnel_ff_cat,means{i});
    elseif i==4
        title("FB Tunnel to Well")
        b=bar(regs_tunnel_fb_cat,means{i});
    elseif i==1
        title("FF Well to Tunnel")
        b=bar(regs_cat,means{i});
    elseif i==2
        title("FB Well to Tunnel")
        b=bar(regs_cat,means{i});
    end
    ax=gca;
    ax.FontSize=24;
    ylim([0.2,0.5])
    % error bars
    [ngroups,nbars]=size(means{i});
    x=nan(nbars,ngroups);
    for j=1:nbars
        x(j,:)=b(j).XEndPoints;
    end
    errorbar(x',means{i},err{i},'k','LineStyle','none')
    
    hold off
end

linkaxes(ax(:),'y')


saveas(gcf,"D:\Brewer lab data\HFS\three stim share figs\rsq_bar.png")
%% Per array bar Rsq

for fid=1:10
    % graphing degrees averaged
    node_degrees_ff_well=[];
    node_degrees_ff_well_se=[];
    node_degrees_fb_well=[];
    node_degrees_fb_well_se=[];
    node_degrees_ff_tunnel=[];
    node_degrees_ff_tunnel_se=[];
    node_degrees_fb_tunnel=[];
    node_degrees_fb_tunnel_se=[];
    
    %relative idx
    nostim_idx=find(nostim_datatable.s_no==fid);
    hfs5_idx=find(hfs5_datatable.s_no==fid);
    hfs40_idx=find(hfs40_datatable.s_no==fid);
    
    %well in/out
    for regi=1:4
        
        %out
        %ff wells out
        nostim_well_ff_degs=nonzeros(cell2mat(nostim_well_ff{nostim_idx,regi}));
        hfs5_well_ff_degs=nonzeros(cell2mat(hfs5_well_ff{hfs5_idx,regi}));
        hfs40_well_ff_degs=nonzeros(cell2mat(hfs40_well_ff{hfs40_idx,regi}));
        node_degrees_ff_well=[node_degrees_ff_well;...
            [mean(nostim_well_ff_degs),mean(hfs5_well_ff_degs),mean(hfs40_well_ff_degs)]];
        node_degrees_ff_well_se=[node_degrees_ff_well_se;...
            [std(nostim_well_ff_degs)/sqrt(length(nostim_well_ff_degs)),...
            std(hfs5_well_ff_degs)/sqrt(length(hfs5_well_ff_degs)),...
            std(hfs40_well_ff_degs)/sqrt(length(hfs40_well_ff_degs))]];
        
        %fb wells out
        nostim_well_fb_degs=nonzeros(cell2mat(nostim_well_fb{nostim_idx,regi}));
        hfs5_well_fb_degs=nonzeros(cell2mat(hfs5_well_fb{hfs5_idx,regi}));
        hfs40_well_fb_degs=nonzeros(cell2mat(hfs40_well_fb{hfs40_idx,regi}));
        node_degrees_fb_well=[node_degrees_fb_well;...
            [mean(nostim_well_fb_degs),mean(hfs5_well_fb_degs),mean(hfs40_well_fb_degs)]];
        node_degrees_fb_well_se=[node_degrees_fb_well_se;...
            [std(nostim_well_fb_degs)/sqrt(length(nostim_well_fb_degs)),...
            std(hfs5_well_fb_degs)/sqrt(length(hfs5_well_fb_degs)),...
            std(hfs40_well_fb_degs)/sqrt(length(hfs40_well_fb_degs))]];
        
        %in
        %ff in
        nostim_tunnel_ff_degs=nonzeros(cell2mat(nostim_tunnel_ff{nostim_idx,regi}));
        hfs5_tunnel_ff_degs=nonzeros(cell2mat(hfs5_tunnel_ff{hfs5_idx,regi}));
        hfs40_tunnel_ff_degs=nonzeros(cell2mat(hfs40_tunnel_ff{hfs40_idx,regi}));
        node_degrees_ff_tunnel=[node_degrees_ff_tunnel;...
            [mean(nostim_tunnel_ff_degs),mean(hfs5_tunnel_ff_degs),mean(hfs40_tunnel_ff_degs)]];
        node_degrees_ff_tunnel_se=[node_degrees_ff_tunnel_se;...
            [std(nostim_tunnel_ff_degs)/sqrt(length(nostim_tunnel_ff_degs)),...
            std(hfs5_tunnel_ff_degs)/sqrt(length(hfs5_tunnel_ff_degs)),...
            std(hfs40_tunnel_ff_degs)/sqrt(length(hfs40_tunnel_ff_degs))]];
        
        %fb in
        nostim_tunnel_fb_degs=nonzeros(cell2mat(nostim_tunnel_fb{nostim_idx,regi}));
        hfs5_tunnel_fb_degs=nonzeros(cell2mat(hfs5_tunnel_fb{hfs5_idx,regi}));
        hfs40_tunnel_fb_degs=nonzeros(cell2mat(hfs40_tunnel_fb{hfs40_idx,regi}));
        node_degrees_fb_tunnel=[node_degrees_fb_tunnel;...
            [mean(nostim_tunnel_fb_degs),mean(hfs5_tunnel_fb_degs),mean(hfs40_tunnel_fb_degs)]];
        node_degrees_fb_tunnel_se=[node_degrees_fb_tunnel_se;...
            [std(nostim_tunnel_fb_degs)/sqrt(length(nostim_tunnel_fb_degs)),...
            std(hfs5_tunnel_fb_degs)/sqrt(length(hfs5_tunnel_fb_degs)),...
            std(hfs40_tunnel_fb_degs)/sqrt(length(hfs40_tunnel_fb_degs))]];
        
        % ff in anova
        name_mat=[repmat("No Stim",length(nostim_well_ff_degs),1);repmat("5 HFS",length(hfs5_well_ff_degs),1);...
            repmat("40 HFS",length(hfs40_well_ff_degs),1)];
%         [~,~,stats]=anova1([nostim_well_ff_degs;hfs5_well_ff_degs;hfs40_well_ff_degs],name_mat);
%         [c,~]=multcompare(stats);
%         ff_in_well_c{regi}=c;
        
        % fb in anova
        name_mat=[repmat("No Stim",length(nostim_well_fb_degs),1);repmat("5 HFS",length(hfs5_well_fb_degs),1);...
            repmat("40 HFS",length(hfs40_well_fb_degs),1)];
%         [~,~,stats]=anova1([nostim_well_fb_degs;hfs5_well_fb_degs;hfs40_well_fb_degs],name_mat);
%         [c,~]=multcompare(stats);
%         fb_in_well_c{regi}=c;
        
        % ff out anova
        name_mat=[repmat("No Stim",length(nostim_tunnel_ff_degs),1);repmat("5 HFS",length(hfs5_tunnel_ff_degs),1);...
            repmat("40 HFS",length(hfs40_tunnel_ff_degs),1)];
%         [~,~,stats]=anova1([nostim_tunnel_ff_degs;hfs5_tunnel_ff_degs;hfs40_tunnel_ff_degs],name_mat);
%         [c,~]=multcompare(stats);
%         ff_out_well_c{regi}=c;
        
        % fb in anova
        name_mat=[repmat("No Stim",length(nostim_tunnel_fb_degs),1);repmat("5 HFS",length(hfs5_tunnel_fb_degs),1);...
            repmat("40 HFS",length(hfs40_tunnel_fb_degs),1)];
%         [~,~,stats]=anova1([nostim_tunnel_fb_degs;hfs5_tunnel_fb_degs;hfs40_tunnel_fb_degs],name_mat);
%         [c,~]=multcompare(stats);
%         fb_out_well_c{regi}=c;
    end
    
    degs={node_degrees_ff_well,node_degrees_fb_well,node_degrees_ff_tunnel,node_degrees_fb_tunnel};
    degs_se={node_degrees_ff_well_se,node_degrees_fb_well_se,node_degrees_ff_tunnel_se,node_degrees_fb_tunnel_se};
    
    figure('units','normalized','outerposition',[0 0 1 1])
    t=tiledlayout(1,4,'TileSpacing','compact');
    ylabel(t,"Degrees Per Node",'FontSize',24)
    xlabel(t,"Subregion",'FontSize',24)
    title(t,"R^2 Weights FID "+string(fid),'FontSize',24)
    for i=1:4
        ax(i)=nexttile;
        hold on
        regs=["EC","DG","CA3","CA1"];
        regs_cat=categorical(regs);
        regs_cat=reordercats(regs_cat,regs);
        regs_tunnel_ff=["EC-DG","DG-CA3","CA3-CA1","CA1-EC"];
        regs_tunnel_fb=["DG-EC","CA3-DG","CA1-CA3","EC-CA1"];
        
        regs_tunnel_ff_cat=reordercats(categorical(regs_tunnel_ff),regs_tunnel_ff);
        regs_tunnel_fb_cat=reordercats(categorical(regs_tunnel_fb),regs_tunnel_fb);
        
        if i==3
            title("Feedforward Tunnel",'FontSize',20)
            b=bar(regs_cat,degs{i});
        elseif i==4
            title("Feedback Tunnel",'FontSize',20)
            b=bar(regs_cat,degs{i});
        elseif i==1
            title("Feedforward Well",'FontSize',20)
            b=bar(regs_cat,degs{i});
        elseif i==2
            title("Feedback Well",'FontSize',20)
            b=bar(regs_cat,degs{i});
        end
        % error bars
        [ngroups,nbars]=size(degs{i});
        x=nan(nbars,ngroups);
        for j=1:nbars
            x(j,:)=b(j).XEndPoints;
        end
        errorbar(x',degs{i},degs_se{i},'k','LineStyle','none')
        
        hold off
    end
    linkaxes(ax(:),'y')
    colororder(colors)
    saveas(gcf,"D:\Brewer lab data\HFS\Temporal Analysis\Graph Analysis\rsq weights nonzero FID "+string(fid)+".png")
end
%% Create Legend Slope
create_legend_slope([0.5:(1.2-0.5)/4:1.2],2,50,0.49)
exportgraphics(gcf,"D:\Brewer lab data\HFS\Temporal Analysis\Graph Analysis\slope 8pt legend.png",'Resolution',1500)

%% Create Legend Rsq
create_legend_slope([0.25:(0.45-0.25)/4:0.45],1,100,0.24)
exportgraphics(gcf,"D:\Brewer lab data\HFS\Temporal Analysis\Graph Analysis\rsq 8pt legend.png",'Resolution',1500)

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

function avg=per_array_avg(weight_tbl,reg)

avg=weight_tbl;
nanmean=@(x)mean(x,'omitnan');
avg{:,[1:length(reg)]}=num2cell(cellfun(@nanmean,table2cell(avg(:,[1:length(reg)]))));

% for i=1:length(reg)
%     avg{:,i}=num2cell(mean(table2cell(avg(:,i))));
% end

end

function [avg_table]=mean_reg(weight_tbl,reg,tun_reg,is_well_source,is_ff)
if ~is_well_source
    is_ff=~is_ff;
end

avg_table=[];
avg=[];
node_order=[];
if ~is_ff
    reg=circshift(reg,-1);
end

for i=1:length(reg)
    avg(1,i)=mean(cell2mat(weight_tbl{:,i}),1,'omitnan');
    
    if is_well_source
        node_order=[node_order;[reg(i),tun_reg(i)]];
    else
        node_order=[node_order;[tun_reg(i),reg(i)]];
    end
end

avg_table=table(avg',node_order,'VariableNames',{'avg_wt','node_order'});

end

function G=construct_G(well_ff,well_fb,tunnel_ff,tunnel_fb)

G=[];
No=[1:8];
Name={'EC','EC-DG','DG','DG-CA3','CA3','CA3-CA1','CA1','CA1-EC'};
% X=[1 2 3 3 3 2 1 1];
X=[1 2 3 2.75 3 2 1 1.25];
% Y=[3 3 3 2 1 1 1 2];
Y=[3 2.75 3 2 1 1.25 1 2];
NodeTable=table(No',Name',X',Y','VariableNames',{'No','Name','X','Y'});

s=[1 2 3 4 5 6 7 8 2 3 4 5 6 7 8 1];
t=[2 3 4 5 6 7 8 1 1 2 3 4 5 6 7 8];

G=digraph(s,t,1,NodeTable);
% digraph_order=string(G.Edges.EndNodes);
digraph_order=[string(Name(s))',string(Name(t))'];

wts=[well_ff;well_fb;tunnel_ff;tunnel_fb];
wts_nodes=string(wts.node_order);
[~,order_wts]=ismember(digraph_order,wts_nodes,'rows');

% well_wts=[well_ff.avg_wt;well_fb.avg_wt];
% well_wts=well_wts(:);
% tunnel_wts=[tunnel_ff.avg_wt;tunnel_fb.avg_wt];
% tunnel_wts=tunnel_wts(:);

weights=wts.avg_wt(order_wts);

G=digraph(s,t,weights,NodeTable);

code=[repmat(0.1,length(well_ff.avg_wt),1);repmat(0.6,length(well_fb.avg_wt),1)];
code=[code;code];
digraph_order=string(G.Edges.EndNodes);
[~,order_wts]=ismember(digraph_order,wts_nodes,'rows');
code=code(order_wts);
code=table(code,'VariableNames',{'Code'});
Edges=[G.Edges,code];

G=digraph(Edges,NodeTable);
end

function plot_8pt_graph(G,graph_type)
figure('Renderer','painters','Position',[10 10 800 800])

if graph_type=='slope'
    max_lwidth=2;
    %SCALING FACTOR FOR LINE WEIGHTS IN FRONT
    LWidths=50*(G.Edges.Weight-0.5)/max_lwidth;
else
    max_lwidth=1;
    %SCALING FACTOR FOR LINE WEIGHTS IN FRONT
    LWidths=100*(G.Edges.Weight-0.25)/max_lwidth;
end




% Plotting edges with color and linewidths
edge_obj = plot(G,'XData',G.Nodes.X,'YData',G.Nodes.Y,...
    'NodeColor','black');
%graph_obj.EdgeLabel = round(G.Edges.Weight,2);
% edge_obj.NodeLabel = G.Nodes.Name;
edge_obj.NodeLabel = {};
edge_obj.NodeFontSize=16;
edge_obj.MarkerSize=10;
%     arrow_obj.ShowArrows = 'off';
edge_obj.LineWidth = LWidths;
edge_obj.EdgeCData = G.Edges.Code;

hold on
colormap(flag(4))
caxis([0 1])
axis square
% xlim([0.5 12.5]), ylim([0.5 12.5])

% Plotting Arrows with blacks
arrow_obj = plot(G,'XData',G.Nodes.X,'YData',G.Nodes.Y,...
    'NodeColor','black');
arrow_obj.LineStyle = ':';
arrow_obj.EdgeColor = 'black';
arrow_obj.ArrowSize = 10;
arrow_obj.NodeFontSize=16;
arrow_obj.NodeLabel={};

ax = gca;
set(gca,'XTick',[],'YTick',[])
end

function [c_all,means_all,se_all]=calc_anova(w1,w2,w3,reg)

means_all=[];
se_all=[];

for i=1:length(reg)
    w=[cell2mat(w1{:,i});cell2mat(w2{:,i});cell2mat(w3{:,i})];
    g=[repmat("g1",length(cell2mat(w1{:,i})),1);...
        repmat("g2",length(cell2mat(w2{:,i})),1);repmat("g3",length(cell2mat(w3{:,i})),1)];
    [~,~,stats]=anova1(w,g);
    [c,means]=multcompare(stats);
    c_all{i}=c;
    means_all=[means_all;means(:,1)'];
    se=[means(1,2),means(2,2),means(3,2)];
    se_all=[se_all;se];
end

end

function create_legend_slope(weights_vec,max_wt,wt_mult,scale)
figure('Position',[100,100,200,200])
x=[1 2 1 2 1 2 1 2 1 2];
y=[1 1 2 2 3 3 4 4 5 5];
% weights_slope=[0.4:0.4:2];
nodeTable_legend=table([1:length(x)]',x',y','VariableNames',{'No','X','Y'});
endNodes_legend=[1 2; 3 4; 5 6; 7 8; 9 10];
edgetable=table(endNodes_legend,weights_vec','VariableNames',{'EndNodes','Weight'});

G_legend=graph(edgetable,nodeTable_legend);
timedep_legend=plot(G_legend,'XData',x,'YData',y,'NodeLabel',...
    [],'NodeColor','none');

timedep_legend.LineWidth=((weights_vec-scale)*wt_mult)/max_wt;
timedep_legend.EdgeColor=[1 0 0];
ylim([0.5,5.5])
set(gca,'Visible','off')
% set(gca,'XColor','none','YColor','none')
axis tight
end