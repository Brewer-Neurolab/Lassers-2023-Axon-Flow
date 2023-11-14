%% Create Graphs with strength of connections varying with time representing slopes of AFR between nodes
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
%% Graph Analysis Calc
% No need to run this section if these data structures already exist

% [nostim_slope_G,nostim_rsq_G]=cat_edges(nostim_allregion_unit_matched,nostim_source,nostim_target,NodeTable,nostim_datatable);
% disp("No stim computed.")
% [hfs5_slope_G,hfs5_rsq_G]=cat_edges(hfs5_allregion_unit_matched,hfs5_source,hfs5_target,NodeTable,hfs5_datatable);
% disp("5 HFS computed.")
% [hfs40_slope_G,hfs40_rsq_G]=cat_edges(hfs5_allregion_unit_matched,hfs40_source,hfs40_target,NodeTable,hfs40_datatable);
% disp("40 HFS computed.")
% 
% save("D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output Full Index slide 10sLM\nostim_slope_G.mat","nostim_slope_G")
% save("D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output Full Index slide 10sLM\nostim_rsq_G.mat","nostim_rsq_G")
% 
% save("D:\Brewer lab data\HFS\Temporal Analysis\5 HFS\AFR Output Full Index slide 10sLM\hfs5_slope_G.mat","hfs5_slope_G")
% save("D:\Brewer lab data\HFS\Temporal Analysis\5 HFS\AFR Output Full Index slide 10sLM\hfs5_rsq_G.mat","hfs5_rsq_G")
% 
% save("D:\Brewer lab data\HFS\Temporal Analysis\40 HFS\AFR Output Full Index slide 10sLM\hfs40_slope_G.mat","hfs40_slope_G")
% save("D:\Brewer lab data\HFS\Temporal Analysis\40 HFS\AFR Output Full Index slide 10sLM\hfs40_rsq_G.mat","hfs40_rsq_G")

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

% [ff_source_nostim,ff_target_nostim,fb_source_nostim,fb_target_nostim] = ...
%     split_dir_nodes(nostim_slope_G,nostim_source,nostim_target,...
%     nostim_data_fldr,nostim_datatable);

%% Graph analysis

%% get degrees
[nostim_well_out_ff,nostim_well_out_fb,nostim_tunnel_out_ff,nostim_tunnel_out_fb]=...
    get_degrees(nostim_slope_G,nostim_datatable,nostim_data_fldr,nostim_temporal_fldr);

[hfs5_well_out_ff,hfs5_well_out_fb,hfs5_tunnel_out_ff,hfs5_tunnel_out_fb]=...
    get_degrees(hfs5_slope_G,hfs5_datatable,hfs5_data_fldr,hfs5_temporal_fldr);

[hfs40_well_out_ff,hfs40_well_out_fb,hfs40_tunnel_out_ff,hfs40_tunnel_out_fb]=...
    get_degrees(hfs40_slope_G,hfs40_datatable,hfs40_data_fldr,hfs40_temporal_fldr);


% graphing degrees averaged
node_degrees_ff_well=[];
node_degrees_ff_well_se=[];
node_degrees_fb_well=[];
node_degrees_fb_well_se=[];
node_degrees_ff_tunnel=[];
node_degrees_ff_tunnel_se=[];
node_degrees_fb_tunnel=[];
node_degrees_fb_tunnel_se=[];

for regi=1:4
    
    %ff wells
    nostim_well_ff_degs=cell2mat(nostim_well_out_ff{:,regi});
    hfs5_well_ff_degs=cell2mat(hfs5_well_out_ff{:,regi});
    hfs40_well_ff_degs=cell2mat(hfs40_well_out_ff{:,regi});
    node_degrees_ff_well=[node_degrees_ff_well;...
        [mean(nostim_well_ff_degs),mean(hfs5_well_ff_degs),mean(hfs40_well_ff_degs)]];
    node_degrees_ff_well_se=[node_degrees_ff_well_se;...
        [std(nostim_well_ff_degs)/sqrt(length(nostim_well_ff_degs)),...
        std(hfs5_well_ff_degs)/sqrt(length(hfs5_well_ff_degs)),...
        std(hfs40_well_ff_degs)/sqrt(length(hfs40_well_ff_degs))]];
    
    %fb wells
    nostim_well_fb_degs=cell2mat(nostim_well_out_fb{:,regi});
    hfs5_well_fb_degs=cell2mat(hfs5_well_out_fb{:,regi});
    hfs40_well_fb_degs=cell2mat(hfs40_well_out_fb{:,regi});
    node_degrees_fb_well=[node_degrees_fb_well;...
        [mean(nostim_well_fb_degs),mean(hfs5_well_fb_degs),mean(hfs40_well_fb_degs)]];
    node_degrees_fb_well_se=[node_degrees_fb_well_se;...
        [std(nostim_well_fb_degs)/sqrt(length(nostim_well_fb_degs)),...
        std(hfs5_well_fb_degs)/sqrt(length(hfs5_well_fb_degs)),...
        std(hfs40_well_fb_degs)/sqrt(length(hfs40_well_fb_degs))]];
    
    %ff tunnel
    nostim_tunnel_ff_degs=cell2mat(nostim_tunnel_out_ff{:,regi});
    hfs5_tunnel_ff_degs=cell2mat(hfs5_tunnel_out_ff{:,regi});
    hfs40_tunnel_ff_degs=cell2mat(hfs40_tunnel_out_ff{:,regi});
    node_degrees_ff_tunnel=[node_degrees_ff_tunnel;...
        [mean(nostim_tunnel_ff_degs),mean(hfs5_tunnel_ff_degs),mean(hfs40_tunnel_ff_degs)]];
    node_degrees_ff_tunnel_se=[node_degrees_ff_tunnel_se;...
        [std(nostim_tunnel_ff_degs)/sqrt(length(nostim_tunnel_ff_degs)),...
        std(hfs5_tunnel_ff_degs)/sqrt(length(hfs5_tunnel_ff_degs)),...
        std(hfs40_tunnel_ff_degs)/sqrt(length(hfs40_tunnel_ff_degs))]];
    
    %fb tunnel
    nostim_tunnel_fb_degs=cell2mat(nostim_tunnel_out_fb{:,regi});
    hfs5_tunnel_fb_degs=cell2mat(hfs5_tunnel_out_fb{:,regi});
    hfs40_tunnel_fb_degs=cell2mat(hfs40_tunnel_out_fb{:,regi});
    node_degrees_fb_tunnel=[node_degrees_fb_tunnel;...
        [mean(nostim_tunnel_fb_degs),mean(hfs5_tunnel_fb_degs),mean(hfs40_tunnel_fb_degs)]];
    node_degrees_fb_tunnel_se=[node_degrees_fb_tunnel_se;...
        [std(nostim_tunnel_fb_degs)/sqrt(length(nostim_tunnel_fb_degs)),...
        std(hfs5_tunnel_fb_degs)/sqrt(length(hfs5_tunnel_fb_degs)),...
        std(hfs40_tunnel_fb_degs)/sqrt(length(hfs40_tunnel_fb_degs))]];
end

degs={node_degrees_ff_well,node_degrees_fb_well,node_degrees_ff_tunnel,node_degrees_fb_tunnel};
degs_se={node_degrees_ff_well_se,node_degrees_fb_well_se,node_degrees_ff_tunnel_se,node_degrees_fb_tunnel_se};

figure('Position',[0,0,800,800])
t=tiledlayout(2,2,'TileSpacing','compact');
ylabel(t,"Out Degrees")
xlabel(t,"Subregion")
title(t,"Out Degrees")
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
        title("Feedforward Tunnel Degs")
        b=bar(regs_tunnel_ff_cat,degs{i});
    elseif i==4
        title("Feedback Tunnel Degs")
        b=bar(regs_tunnel_fb_cat,degs{i});
    elseif i==1
        title("Feedforward Well Degs")
        b=bar(regs_cat,degs{i});
    elseif i==2
        title("Feedback Well Degs")
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
% saveas(gcf,"D:\Brewer lab data\HFS\Temporal Analysis\Graph Analysis\degrees.png")
%% get degrees non-zeros
%out edges
[nostim_well_out_ff,nostim_well_out_fb,nostim_tunnel_out_ff,nostim_tunnel_out_fb]=...
    get_out_degrees(nostim_slope_G,nostim_datatable,nostim_data_fldr,nostim_temporal_fldr);

[hfs5_well_out_ff,hfs5_well_out_fb,hfs5_tunnel_out_ff,hfs5_tunnel_out_fb]=...
    get_out_degrees(hfs5_slope_G,hfs5_datatable,hfs5_data_fldr,hfs5_temporal_fldr);

[hfs40_well_out_ff,hfs40_well_out_fb,hfs40_tunnel_out_ff,hfs40_tunnel_out_fb]=...
    get_out_degrees(hfs40_slope_G,hfs40_datatable,hfs40_data_fldr,hfs40_temporal_fldr);
% in edges
[nostim_well_in_ff,nostim_well_in_fb,nostim_tunnel_in_ff,nostim_tunnel_in_fb]=...
    get_in_degrees(nostim_slope_G,nostim_datatable,nostim_data_fldr,nostim_temporal_fldr);

[hfs5_well_in_ff,hfs5_well_in_fb,hfs5_tunnel_in_ff,hfs5_tunnel_in_fb]=...
    get_in_degrees(hfs5_slope_G,hfs5_datatable,hfs5_data_fldr,hfs5_temporal_fldr);

[hfs40_well_in_ff,hfs40_well_in_fb,hfs40_tunnel_in_ff,hfs40_tunnel_in_fb]=...
    get_in_degrees(hfs40_slope_G,hfs40_datatable,hfs40_data_fldr,hfs40_temporal_fldr);

% graphing degrees averaged
node_degrees_ff_well=[];
node_degrees_ff_well_se=[];
node_degrees_fb_well=[];
node_degrees_fb_well_se=[];
node_degrees_ff_tunnel=[];
node_degrees_ff_tunnel_se=[];
node_degrees_fb_tunnel=[];
node_degrees_fb_tunnel_se=[];


%% well in/out
ff_in_well_m=[];
ff_in_well_se=[];
fb_in_well_m=[];
fb_in_well_se=[];
ff_out_well_m=[];
ff_out_well_se=[];
fb_out_well_m=[];
fb_out_well_se=[];

for regi=1:4
    
    %out
    %ff wells out
    nostim_well_ff_degs=nonzeros(cell2mat(nostim_well_in_ff{:,regi}));
    hfs5_well_ff_degs=nonzeros(cell2mat(hfs5_well_in_ff{:,regi}));
    hfs40_well_ff_degs=nonzeros(cell2mat(hfs40_well_in_ff{:,regi}));
    node_degrees_ff_well=[node_degrees_ff_well;...
        [mean(nostim_well_ff_degs),mean(hfs5_well_ff_degs),mean(hfs40_well_ff_degs)]];
    node_degrees_ff_well_se=[node_degrees_ff_well_se;...
        [std(nostim_well_ff_degs)/sqrt(length(nostim_well_ff_degs)),...
        std(hfs5_well_ff_degs)/sqrt(length(hfs5_well_ff_degs)),...
        std(hfs40_well_ff_degs)/sqrt(length(hfs40_well_ff_degs))]];
    
    %fb wells out
    nostim_well_fb_degs=nonzeros(cell2mat(nostim_well_in_fb{:,regi}));
    hfs5_well_fb_degs=nonzeros(cell2mat(hfs5_well_in_fb{:,regi}));
    hfs40_well_fb_degs=nonzeros(cell2mat(hfs40_well_in_fb{:,regi}));
    node_degrees_fb_well=[node_degrees_fb_well;...
        [mean(nostim_well_fb_degs),mean(hfs5_well_fb_degs),mean(hfs40_well_fb_degs)]];
    node_degrees_fb_well_se=[node_degrees_fb_well_se;...
        [std(nostim_well_fb_degs)/sqrt(length(nostim_well_fb_degs)),...
        std(hfs5_well_fb_degs)/sqrt(length(hfs5_well_fb_degs)),...
        std(hfs40_well_fb_degs)/sqrt(length(hfs40_well_fb_degs))]];
    
    %in
    %ff in
    nostim_tunnel_ff_degs=nonzeros(cell2mat(nostim_well_out_ff{:,regi}));
    hfs5_tunnel_ff_degs=nonzeros(cell2mat(hfs5_well_out_ff{:,regi}));
    hfs40_tunnel_ff_degs=nonzeros(cell2mat(hfs40_well_out_ff{:,regi}));
    node_degrees_ff_tunnel=[node_degrees_ff_tunnel;...
        [mean(nostim_tunnel_ff_degs),mean(hfs5_tunnel_ff_degs),mean(hfs40_tunnel_ff_degs)]];
    node_degrees_ff_tunnel_se=[node_degrees_ff_tunnel_se;...
        [std(nostim_tunnel_ff_degs)/sqrt(length(nostim_tunnel_ff_degs)),...
        std(hfs5_tunnel_ff_degs)/sqrt(length(hfs5_tunnel_ff_degs)),...
        std(hfs40_tunnel_ff_degs)/sqrt(length(hfs40_tunnel_ff_degs))]];
    
    %fb in
    nostim_tunnel_fb_degs=nonzeros(cell2mat(nostim_well_out_fb{:,regi}));
    hfs5_tunnel_fb_degs=nonzeros(cell2mat(hfs5_well_out_fb{:,regi}));
    hfs40_tunnel_fb_degs=nonzeros(cell2mat(hfs40_well_out_fb{:,regi}));
    node_degrees_fb_tunnel=[node_degrees_fb_tunnel;...
        [mean(nostim_tunnel_fb_degs),mean(hfs5_tunnel_fb_degs),mean(hfs40_tunnel_fb_degs)]];
    node_degrees_fb_tunnel_se=[node_degrees_fb_tunnel_se;...
        [std(nostim_tunnel_fb_degs)/sqrt(length(nostim_tunnel_fb_degs)),...
        std(hfs5_tunnel_fb_degs)/sqrt(length(hfs5_tunnel_fb_degs)),...
        std(hfs40_tunnel_fb_degs)/sqrt(length(hfs40_tunnel_fb_degs))]];
    
    % ff in anova
    name_mat=[repmat("No Stim",length(nostim_well_ff_degs),1);repmat("5 HFS",length(hfs5_well_ff_degs),1);...
        repmat("40 HFS",length(hfs40_well_ff_degs),1)];
    [~,~,stats]=anova1([nostim_well_ff_degs;hfs5_well_ff_degs;hfs40_well_ff_degs],name_mat);
    [c,m]=multcompare(stats);
    ff_in_well_c{regi}=c;
    ff_in_well_m=[ff_in_well_m;m(:,1)'];
    ff_in_well_se=[ff_in_well_se;m(:,2)'];
    
    % fb in anova
    name_mat=[repmat("No Stim",length(nostim_well_fb_degs),1);repmat("5 HFS",length(hfs5_well_fb_degs),1);...
        repmat("40 HFS",length(hfs40_well_fb_degs),1)];
    [~,~,stats]=anova1([nostim_well_fb_degs;hfs5_well_fb_degs;hfs40_well_fb_degs],name_mat);
    [c,m]=multcompare(stats);
    fb_in_well_c{regi}=c;
    fb_in_well_m=[fb_in_well_m;m(:,1)'];
    fb_in_well_se=[fb_in_well_se;m(:,2)'];
    
    % ff out anova
    name_mat=[repmat("No Stim",length(nostim_tunnel_ff_degs),1);repmat("5 HFS",length(hfs5_tunnel_ff_degs),1);...
        repmat("40 HFS",length(hfs40_tunnel_ff_degs),1)];
    [~,~,stats]=anova1([nostim_tunnel_ff_degs;hfs5_tunnel_ff_degs;hfs40_tunnel_ff_degs],name_mat);
    [c,m]=multcompare(stats);
    ff_out_well_c{regi}=c;
    ff_out_well_m=[ff_out_well_m;m(:,1)'];
    ff_out_well_se=[ff_out_well_se;m(:,2)'];
    
    % fb in anova
    name_mat=[repmat("No Stim",length(nostim_tunnel_fb_degs),1);repmat("5 HFS",length(hfs5_tunnel_fb_degs),1);...
        repmat("40 HFS",length(hfs40_tunnel_fb_degs),1)];
    [~,~,stats]=anova1([nostim_tunnel_fb_degs;hfs5_tunnel_fb_degs;hfs40_tunnel_fb_degs],name_mat);
    [c,m]=multcompare(stats);
    fb_out_well_c{regi}=c;
    fb_out_well_m=[fb_out_well_m;m(:,1)'];
    fb_out_well_se=[fb_out_well_se;m(:,2)'];
    
    %in ff fb anova
    v=[nostim_well_ff_degs;nostim_well_fb_degs];
    g=[repmat("ff",length(nostim_well_ff_degs),1);repmat("fb",length(nostim_well_fb_degs),1)];
    [~,~,stats]=anova1(v,g);
    [c,~]=multcompare(stats);
    ff_fb_in_c{regi,1}=c;
    
    v=[hfs5_well_ff_degs;hfs5_well_fb_degs];
    g=[repmat("ff",length(hfs5_well_ff_degs),1);repmat("fb",length(hfs5_well_fb_degs),1)];
    [~,~,stats]=anova1(v,g);
    [c,~]=multcompare(stats);
    ff_fb_in_c{regi,2}=c;
    
    v=[hfs40_well_ff_degs;hfs40_well_fb_degs];
    g=[repmat("ff",length(hfs40_well_ff_degs),1);repmat("fb",length(hfs40_well_fb_degs),1)];
    [~,~,stats]=anova1(v,g);
    [c,~]=multcompare(stats);
    ff_fb_in_c{regi,3}=c;
    
    %out ff fb anova
    v=[nostim_tunnel_ff_degs;nostim_tunnel_fb_degs];
    g=[repmat("ff",length(nostim_tunnel_ff_degs),1);repmat("fb",length(nostim_tunnel_fb_degs),1)];
    [~,~,stats]=anova1(v,g);
    [c,~]=multcompare(stats);
    ff_fb_out_c{regi,1}=c;
    
    v=[hfs5_tunnel_ff_degs;hfs5_tunnel_fb_degs];
    g=[repmat("ff",length(hfs5_tunnel_ff_degs),1);repmat("fb",length(hfs5_tunnel_fb_degs),1)];
    [~,~,stats]=anova1(v,g);
    [c,~]=multcompare(stats);
    ff_fb_out_c{regi,2}=c;
    
    v=[hfs40_tunnel_ff_degs;hfs40_tunnel_fb_degs];
    g=[repmat("ff",length(hfs40_tunnel_ff_degs),1);repmat("fb",length(hfs40_tunnel_fb_degs),1)];
    [~,~,stats]=anova1(v,g);
    [c,~]=multcompare(stats);
    ff_fb_out_c{regi,3}=c;
    
    effect_in{regi}=ff_fb_meanEffect(nostim_well_ff_degs,nostim_well_fb_degs,hfs5_well_ff_degs,hfs5_well_fb_degs,hfs40_well_ff_degs,hfs40_well_fb_degs);
    effect_out{regi}=ff_fb_meanEffect(nostim_tunnel_ff_degs,nostim_tunnel_fb_degs,hfs5_tunnel_ff_degs,hfs5_tunnel_fb_degs,hfs40_tunnel_ff_degs,hfs40_tunnel_fb_degs);
end

ff_fb_c = reshape([ff_fb_in_c(:) ff_fb_out_c(:)]',2*size(ff_fb_in_c,1), []);
get_last=@(x) x(1,end);
ff_fb_p = cellfun(get_last,ff_fb_c);

% degs={node_degrees_ff_well,node_degrees_fb_well,node_degrees_ff_tunnel,node_degrees_fb_tunnel};
% degs_se={node_degrees_ff_well_se,node_degrees_fb_well_se,node_degrees_ff_tunnel_se,node_degrees_fb_tunnel_se};

degs={ff_in_well_m,fb_in_well_m,ff_out_well_m,fb_out_well_m};
degs_se={ff_in_well_se,fb_in_well_se,ff_out_well_se,fb_out_well_se};

figure('Position',[0,0,800,800])
t=tiledlayout(2,2,'TileSpacing','compact');
ylabel(t,"Degrees Per Node",'FontSize',24)
xlabel(t,"Subregion",'FontSize',24)
title(t,"Well Degrees",'FontSize',24)
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
        title("Feedforward Out Degrees",'FontSize',20)
        b=bar(regs_cat,degs{i});
    elseif i==4
        title("Feedback Out Degrees",'FontSize',20)
        b=bar(regs_cat,degs{i});
    elseif i==1
        title("Feedforward In Degrees",'FontSize',20)
        b=bar(regs_cat,degs{i});
    elseif i==2
        title("Feedback In Degrees",'FontSize',20)
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
% saveas(gcf,"D:\Brewer lab data\HFS\Temporal Analysis\Graph Analysis\degrees well nonzero.png")
%% horzbar from anova
dat2plot_ff=[];
dat2plot_ff_se=[];
dat2plot_fb=[];
dat2plot_fb_se=[];
labels={'EC In','EC Out','DG In','DG Out','CA3 In','CA3 Out','CA1 In','CA1 Out'};
figure( 'Position', [100 100 1400 600])
for i=1:4
    dat2plot_ff=[dat2plot_ff;degs{1}(i,:);degs{3}(i,:)];
    dat2plot_ff_se=[dat2plot_ff_se;degs_se{1}(i,:);degs_se{3}(i,:)];
    dat2plot_fb=[dat2plot_fb;degs{2}(i,:);degs{4}(i,:)];
    dat2plot_fb_se=[dat2plot_fb_se;degs_se{2}(i,:);degs_se{4}(i,:)];
end

plot_horzbar_from_anova...
    (dat2plot_ff,dat2plot_ff_se,dat2plot_fb,dat2plot_fb_se,labels,ff_colors,fb_colors);
xlabel("FB Edges                  FF Edges")
xlim([-5,5])
xticks([-5:1:5])
ylim([0,9])
set(gca,'LineWidth',2,'TickLength',[0.02 0.02]);
set(gca,'FontSize',18)
% exportgraphics(gcf,".\three stim share figs\well_degree_comp.png",'Resolution',1500)
exportgraphics(gcf,".\Temporal Analysis\Graph Analysis\well_degree_comp.png",'Resolution',1500)

mean_diff=dat2plot_ff-dat2plot_fb;

%% FF/FB% degrees

ff_fb_perc=(dat2plot_ff./dat2plot_fb)*100;

% ff_fb_SD=dat2plot

ff_fb_SD=[];
ff_fb_lengths=[];
for regi=1:4
    
    %out
    %ff wells out
    nostim_well_ff_degs=nonzeros(cell2mat(nostim_well_in_ff{:,regi}));
    hfs5_well_ff_degs=nonzeros(cell2mat(hfs5_well_in_ff{:,regi}));
    hfs40_well_ff_degs=nonzeros(cell2mat(hfs40_well_in_ff{:,regi}));
    node_degrees_ff_well_m=[mean(nostim_well_ff_degs),mean(hfs5_well_ff_degs),mean(hfs40_well_ff_degs)];
    node_degrees_ff_well_sd=[std(nostim_well_ff_degs),std(hfs5_well_ff_degs),std(hfs40_well_ff_degs)];
    ff_out_lengths=[length(nostim_well_ff_degs),length(hfs5_well_ff_degs),length(hfs40_well_ff_degs)];

    %fb wells out
    nostim_well_fb_degs=nonzeros(cell2mat(nostim_well_in_fb{:,regi}));
    hfs5_well_fb_degs=nonzeros(cell2mat(hfs5_well_in_fb{:,regi}));
    hfs40_well_fb_degs=nonzeros(cell2mat(hfs40_well_in_fb{:,regi}));
    node_degrees_fb_well_m=[mean(nostim_well_fb_degs),mean(hfs5_well_fb_degs),mean(hfs40_well_fb_degs)];
    node_degrees_fb_well_sd=[std(nostim_well_fb_degs),std(hfs5_well_fb_degs),std(hfs40_well_fb_degs)];
    fb_out_lengths=[length(nostim_well_fb_degs),length(hfs5_well_fb_degs),length(hfs40_well_fb_degs)];
    
    %in
    %ff in
    nostim_tunnel_ff_degs=nonzeros(cell2mat(nostim_well_out_ff{:,regi}));
    hfs5_tunnel_ff_degs=nonzeros(cell2mat(hfs5_well_out_ff{:,regi}));
    hfs40_tunnel_ff_degs=nonzeros(cell2mat(hfs40_well_out_ff{:,regi}));
    node_degrees_ff_tunnel_m=[mean(nostim_tunnel_ff_degs),mean(hfs5_tunnel_ff_degs),mean(hfs40_tunnel_ff_degs)];
    node_degrees_ff_tunnel_sd=[std(nostim_tunnel_ff_degs),std(hfs5_tunnel_ff_degs), std(hfs40_tunnel_ff_degs)];
    ff_in_lengths=[length(nostim_tunnel_ff_degs),length(hfs5_tunnel_ff_degs),length(hfs40_tunnel_ff_degs)];

    %fb in
    nostim_tunnel_fb_degs=nonzeros(cell2mat(nostim_well_out_fb{:,regi}));
    hfs5_tunnel_fb_degs=nonzeros(cell2mat(hfs5_well_out_fb{:,regi}));
    hfs40_tunnel_fb_degs=nonzeros(cell2mat(hfs40_well_out_fb{:,regi}));
    node_degrees_fb_tunnel_m=[mean(nostim_tunnel_fb_degs),mean(hfs5_tunnel_fb_degs),mean(hfs40_tunnel_fb_degs)];
    node_degrees_fb_tunnel_sd=[std(nostim_tunnel_fb_degs),std(hfs5_tunnel_fb_degs),std(hfs40_tunnel_fb_degs)];
    fb_in_lengths=[length(nostim_tunnel_fb_degs),length(hfs5_tunnel_fb_degs),length(hfs40_tunnel_fb_degs)];

    ff_fb_SD=[ff_fb_SD;((node_degrees_ff_well_m./node_degrees_fb_well_m)*100).*...
        sqrt((node_degrees_ff_well_sd./node_degrees_fb_well_sd).^2+...
        (node_degrees_fb_well_sd./node_degrees_fb_well_sd).^2);...
        ((node_degrees_ff_tunnel_m./node_degrees_fb_tunnel_m)*100).*...
        sqrt((node_degrees_ff_tunnel_m./node_degrees_fb_tunnel_m).^2+...
        (node_degrees_fb_tunnel_sd./node_degrees_fb_tunnel_sd).^2)];
    ff_fb_lengths=[ff_fb_lengths;(ff_in_lengths+fb_in_lengths);(ff_out_lengths+fb_out_lengths)];
end

labels={'EC-DG In','EC-DG Out','DG-CA3 In','DG-CA3 Out','CA3-CA1 In','CA3-CA1 Out','CA1-EC In','CA1-EC Out'};
regs=categorical(labels);
regs=reordercats(regs,labels);

figure( 'Position', [100 100 1400 600])
b=bar(regs,ff_fb_perc);
b(1).BaseValue=100;
b(1).BaseLine.LineStyle = "--";
b(1).BaseLine.LineWidth=2;

ff_fb_SE=ff_fb_SD./sqrt(ff_fb_lengths);
hold on
% grouped error bars
[ngroups,nbars]=size(ff_fb_perc);
x=[];
for j=1:nbars
    x(j,:)=b(j).XEndPoints;
end
% errorbar(x',ff_fb_perc,ff_fb_SE,'k','LineStyle','none')
hold off
ylabel("% FF/FB")
set(gca,'LineWidth',2,'TickLength',[0.02 0.02]);
set(gca,'FontSize',18)
exportgraphics(gcf,".\Temporal Analysis\Graph Analysis\degrees_ff_fb_perc.png",'Resolution',1500)

%% Cohen's D Bar Graphs

CD_in=[];
CD_out=[];
CD_in_CE_neg=[];
CD_in_CE_pos=[];
CD_out_CE_neg=[];
CD_out_CE_pos=[];
labels={'EC In','EC Out','DG In','DG Out','CA3 In','CA3 Out','CA1 In','CA1 Out'};

for i=1:4
    CD_in=[CD_in;[effect_in{i}{1,1}{1}.Effect,effect_in{i}{1,2}{1}.Effect,effect_in{i}{1,3}{1}.Effect]];
    CD_in_CE_neg=[CD_in_CE_neg;[effect_in{i}{1,1}{1}.ConfidenceIntervals(1),effect_in{i}{1,2}{1}.ConfidenceIntervals(1),effect_in{i}{1,3}{1}.ConfidenceIntervals(1)]];
    CD_in_CE_pos=[CD_in_CE_pos;[effect_in{i}{1,1}{1}.ConfidenceIntervals(2),effect_in{i}{1,2}{1}.ConfidenceIntervals(2),effect_in{i}{1,3}{1}.ConfidenceIntervals(2)]];
    CD_out=[CD_out;[effect_out{i}{1,1}{1}.Effect,effect_out{i}{1,2}{1}.Effect,effect_out{i}{1,3}{1}.Effect]];
    CD_out_CE_neg=[CD_out_CE_neg;[effect_out{i}{1,1}{1}.ConfidenceIntervals(1),effect_out{i}{1,2}{1}.ConfidenceIntervals(1),effect_out{i}{1,3}{1}.ConfidenceIntervals(1)]];
    CD_out_CE_pos=[CD_out_CE_pos;[effect_out{i}{1,1}{1}.ConfidenceIntervals(2),effect_out{i}{1,2}{1}.ConfidenceIntervals(2),effect_out{i}{1,3}{1}.ConfidenceIntervals(2)]];
end

CD_m=interleave_rows(CD_in,CD_out);
CD_CE_neg=interleave_rows(CD_in_CE_neg,CD_out_CE_neg);
CD_CE_pos=interleave_rows(CD_in_CE_pos,CD_out_CE_pos);

figure( 'Position', [100 100 700 600])
plot_horzbar_CohenD...
    (CD_m,CD_CE_neg,CD_CE_pos,labels);
xlabel("Degree Directional Mean Effect Size")
xlim([-3,3])
ylim([0,9])
set(gca,'LineWidth',2,'TickLength',[0.02 0.02]);
set(gca,'FontSize',18)
exportgraphics(gcf,".\Temporal Analysis\Graph Analysis\well_degree_cohenD.png",'Resolution',1500)
%% tunnel in/out
ff_in_well_m=[];
ff_in_well_se=[];
fb_in_well_m=[];
fb_in_well_se=[];
ff_out_well_m=[];
ff_out_well_se=[];
fb_out_well_m=[];
fb_out_well_se=[];

% graphing degrees averaged
node_degrees_ff_well=[];
node_degrees_ff_well_se=[];
node_degrees_fb_well=[];
node_degrees_fb_well_se=[];
node_degrees_ff_tunnel=[];
node_degrees_ff_tunnel_se=[];
node_degrees_fb_tunnel=[];
node_degrees_fb_tunnel_se=[];

for regi=1:4
    
    %ff out tunnel
    nostim_well_ff_degs=nonzeros(cell2mat(nostim_tunnel_in_ff{:,regi}));
    hfs5_well_ff_degs=nonzeros(cell2mat(hfs5_tunnel_in_ff{:,regi}));
    hfs40_well_ff_degs=nonzeros(cell2mat(hfs40_tunnel_in_ff{:,regi}));
    node_degrees_ff_well=[node_degrees_ff_well;...
        [mean(nostim_well_ff_degs),mean(hfs5_well_ff_degs),mean(hfs40_well_ff_degs)]];
    node_degrees_ff_well_se=[node_degrees_ff_well_se;...
        [std(nostim_well_ff_degs)/sqrt(length(nostim_well_ff_degs)),...
        std(hfs5_well_ff_degs)/sqrt(length(hfs5_well_ff_degs)),...
        std(hfs40_well_ff_degs)/sqrt(length(hfs40_well_ff_degs))]];
    
    %fb out tunnel
    nostim_well_fb_degs=nonzeros(cell2mat(nostim_tunnel_in_fb{:,regi}));
    hfs5_well_fb_degs=nonzeros(cell2mat(hfs5_tunnel_in_fb{:,regi}));
    hfs40_well_fb_degs=nonzeros(cell2mat(hfs40_tunnel_in_fb{:,regi}));
    node_degrees_fb_well=[node_degrees_fb_well;...
        [mean(nostim_well_fb_degs),mean(hfs5_well_fb_degs),mean(hfs40_well_fb_degs)]];
    node_degrees_fb_well_se=[node_degrees_fb_well_se;...
        [std(nostim_well_fb_degs)/sqrt(length(nostim_well_fb_degs)),...
        std(hfs5_well_fb_degs)/sqrt(length(hfs5_well_fb_degs)),...
        std(hfs40_well_fb_degs)/sqrt(length(hfs40_well_fb_degs))]];
    
    %ff in tunnel
    nostim_tunnel_ff_degs=nonzeros(cell2mat(nostim_tunnel_out_ff{:,regi}));
    hfs5_tunnel_ff_degs=nonzeros(cell2mat(hfs5_tunnel_out_ff{:,regi}));
    hfs40_tunnel_ff_degs=nonzeros(cell2mat(hfs40_tunnel_out_ff{:,regi}));
    node_degrees_ff_tunnel=[node_degrees_ff_tunnel;...
        [mean(nostim_tunnel_ff_degs),mean(hfs5_tunnel_ff_degs),mean(hfs40_tunnel_ff_degs)]];
    node_degrees_ff_tunnel_se=[node_degrees_ff_tunnel_se;...
        [std(nostim_tunnel_ff_degs)/sqrt(length(nostim_tunnel_ff_degs)),...
        std(hfs5_tunnel_ff_degs)/sqrt(length(hfs5_tunnel_ff_degs)),...
        std(hfs40_tunnel_ff_degs)/sqrt(length(hfs40_tunnel_ff_degs))]];
    
    %fb in tunnel
    nostim_tunnel_fb_degs=nonzeros(cell2mat(nostim_tunnel_out_fb{:,regi}));
    hfs5_tunnel_fb_degs=nonzeros(cell2mat(hfs5_tunnel_out_fb{:,regi}));
    hfs40_tunnel_fb_degs=nonzeros(cell2mat(hfs40_tunnel_out_fb{:,regi}));
    node_degrees_fb_tunnel=[node_degrees_fb_tunnel;...
        [mean(nostim_tunnel_fb_degs),mean(hfs5_tunnel_fb_degs),mean(hfs40_tunnel_fb_degs)]];
    node_degrees_fb_tunnel_se=[node_degrees_fb_tunnel_se;...
        [std(nostim_tunnel_fb_degs)/sqrt(length(nostim_tunnel_fb_degs)),...
        std(hfs5_tunnel_fb_degs)/sqrt(length(hfs5_tunnel_fb_degs)),...
        std(hfs40_tunnel_fb_degs)/sqrt(length(hfs40_tunnel_fb_degs))]];
    
   % ff in anova
    name_mat=[repmat("No Stim",length(nostim_well_ff_degs),1);repmat("5 HFS",length(hfs5_well_ff_degs),1);...
        repmat("40 HFS",length(hfs40_well_ff_degs),1)];
    [~,~,stats]=anova1([nostim_well_ff_degs;hfs5_well_ff_degs;hfs40_well_ff_degs],name_mat);
    [c,m]=multcompare(stats);
    ff_in_well_c{regi}=c;
    ff_in_well_m=[ff_in_well_m;m(:,1)'];
    ff_in_well_se=[ff_in_well_se;m(:,2)'];
    
    % fb in anova
    name_mat=[repmat("No Stim",length(nostim_well_fb_degs),1);repmat("5 HFS",length(hfs5_well_fb_degs),1);...
        repmat("40 HFS",length(hfs40_well_fb_degs),1)];
    [~,~,stats]=anova1([nostim_well_fb_degs;hfs5_well_fb_degs;hfs40_well_fb_degs],name_mat);
    [c,m]=multcompare(stats);
    fb_in_well_c{regi}=c;
    fb_in_well_m=[fb_in_well_m;m(:,1)'];
    fb_in_well_se=[fb_in_well_se;m(:,2)'];
    
    % ff out anova
    name_mat=[repmat("No Stim",length(nostim_tunnel_ff_degs),1);repmat("5 HFS",length(hfs5_tunnel_ff_degs),1);...
        repmat("40 HFS",length(hfs40_tunnel_ff_degs),1)];
    [~,~,stats]=anova1([nostim_tunnel_ff_degs;hfs5_tunnel_ff_degs;hfs40_tunnel_ff_degs],name_mat);
    [c,m]=multcompare(stats);
    ff_out_well_c{regi}=c;
    ff_out_well_m=[ff_out_well_m;m(:,1)'];
    ff_out_well_se=[ff_out_well_se;m(:,2)'];
    
    % fb in anova
    name_mat=[repmat("No Stim",length(nostim_tunnel_fb_degs),1);repmat("5 HFS",length(hfs5_tunnel_fb_degs),1);...
        repmat("40 HFS",length(hfs40_tunnel_fb_degs),1)];
    [~,~,stats]=anova1([nostim_tunnel_fb_degs;hfs5_tunnel_fb_degs;hfs40_tunnel_fb_degs],name_mat);
    [c,m]=multcompare(stats);
    fb_out_well_c{regi}=c;
    fb_out_well_m=[fb_out_well_m;m(:,1)'];
    fb_out_well_se=[fb_out_well_se;m(:,2)'];
    
    %in ff fb anova
    v=[nostim_well_ff_degs;nostim_well_fb_degs];
    g=[repmat("ff",length(nostim_well_ff_degs),1);repmat("fb",length(nostim_well_fb_degs),1)];
    [~,~,stats]=anova1(v,g);
    [c,~]=multcompare(stats);
    ff_fb_in_c{regi,1}=c;
    
    v=[hfs5_well_ff_degs;hfs5_well_fb_degs];
    g=[repmat("ff",length(hfs5_well_ff_degs),1);repmat("fb",length(hfs5_well_fb_degs),1)];
    [~,~,stats]=anova1(v,g);
    [c,~]=multcompare(stats);
    ff_fb_in_c{regi,2}=c;
    
    v=[hfs40_well_ff_degs;hfs40_well_fb_degs];
    g=[repmat("ff",length(hfs40_well_ff_degs),1);repmat("fb",length(hfs40_well_fb_degs),1)];
    [~,~,stats]=anova1(v,g);
    [c,~]=multcompare(stats);
    ff_fb_in_c{regi,3}=c;
    
    %out ff fb anova
    v=[nostim_tunnel_ff_degs;nostim_tunnel_fb_degs];
    g=[repmat("ff",length(nostim_tunnel_ff_degs),1);repmat("fb",length(nostim_tunnel_fb_degs),1)];
    [~,~,stats]=anova1(v,g);
    [c,~]=multcompare(stats);
    ff_fb_out_c{regi,1}=c;
    
    v=[hfs5_tunnel_ff_degs;hfs5_tunnel_fb_degs];
    g=[repmat("ff",length(hfs5_tunnel_ff_degs),1);repmat("fb",length(hfs5_tunnel_fb_degs),1)];
    [~,~,stats]=anova1(v,g);
    [c,~]=multcompare(stats);
    ff_fb_out_c{regi,2}=c;
    
    v=[hfs40_tunnel_ff_degs;hfs40_tunnel_fb_degs];
    g=[repmat("ff",length(hfs40_tunnel_ff_degs),1);repmat("fb",length(hfs40_tunnel_fb_degs),1)];
    [~,~,stats]=anova1(v,g);
    [c,~]=multcompare(stats);
    ff_fb_out_c{regi,3}=c;

    effect_in{regi}=ff_fb_meanEffect(nostim_well_ff_degs,nostim_well_fb_degs,hfs5_well_ff_degs,hfs5_well_fb_degs,hfs40_well_ff_degs,hfs40_well_fb_degs);
    effect_out{regi}=ff_fb_meanEffect(nostim_tunnel_ff_degs,nostim_tunnel_fb_degs,hfs5_tunnel_ff_degs,hfs5_tunnel_fb_degs,hfs40_tunnel_ff_degs,hfs40_tunnel_fb_degs);
end

ff_fb_c = reshape([ff_fb_in_c(:) ff_fb_out_c(:)]',2*size(ff_fb_in_c,1), []);
get_last=@(x) x(1,end);
ff_fb_p = cellfun(get_last,ff_fb_c);

% degs={node_degrees_ff_well,node_degrees_fb_well,node_degrees_ff_tunnel,node_degrees_fb_tunnel};
% degs_se={node_degrees_ff_well_se,node_degrees_fb_well_se,node_degrees_ff_tunnel_se,node_degrees_fb_tunnel_se};

degs={ff_in_well_m,fb_in_well_m,ff_out_well_m,fb_out_well_m};
degs_se={ff_in_well_se,fb_in_well_se,ff_out_well_se,fb_out_well_se};

figure('Position',[0,0,800,800])
t=tiledlayout(2,2,'TileSpacing','compact');
ylabel(t,"Degrees Per Node",'FontSize',24)
xlabel(t,"Subregion",'FontSize',24)
title(t,"Tunnel Degrees",'FontSize',24)
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
        title("Feedforward Out Degrees",'FontSize',20)
        b=bar(regs_tunnel_ff_cat,degs{i});
    elseif i==4
        title("Feedback Out Degrees",'FontSize',20)
        b=bar(regs_tunnel_fb_cat,degs{i});
    elseif i==1
        title("Feedforward In Degrees",'FontSize',20)
        b=bar(regs_tunnel_ff_cat,degs{i});
    elseif i==2
        title("Feedback In Degrees",'FontSize',20)
        b=bar(regs_tunnel_fb_cat,degs{i});
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
% saveas(gcf,"D:\Brewer lab data\HFS\Temporal Analysis\Graph Analysis\degrees tunnel nonzero.png")
%% horzbar from anova
dat2plot_ff=[];
dat2plot_ff_se=[];
dat2plot_fb=[];
dat2plot_fb_se=[];
labels={'EC-DG In','EC-DG Out','DG-CA3 In','DG-CA3 Out','CA3-CA1 In','CA3-CA1 Out','CA1-EC In','CA1-EC Out'};
figure( 'Position', [100 100 1400 600])
for i=1:4
    dat2plot_ff=[dat2plot_ff;degs{1}(i,:);degs{3}(i,:)];
    dat2plot_ff_se=[dat2plot_ff_se;degs_se{1}(i,:);degs_se{3}(i,:)];
    dat2plot_fb=[dat2plot_fb;degs{2}(i,:);degs{4}(i,:)];
    dat2plot_fb_se=[dat2plot_fb_se;degs_se{2}(i,:);degs_se{4}(i,:)];
end

plot_horzbar_from_anova...
    (dat2plot_ff,dat2plot_ff_se,dat2plot_fb,dat2plot_fb_se,labels,ff_colors,fb_colors)
xlabel("FB Edges                  FF Edges")
xlim([-16,16])
% xticks([-5:1:5])
ylim([0,9])
set(gca,'LineWidth',2,'TickLength',[0.02 0.02]);
set(gca,'FontSize',18)
% exportgraphics(gcf,".\three stim share figs\tunnel_degree_comp.png",'Resolution',1500)
exportgraphics(gcf,".\Temporal Analysis\Graph Analysis\tunnel_degree_comp.png",'Resolution',1500)

% use ff_fb_p and mean diff to decide where to put the +/-FF on horz bars
mean_diff=dat2plot_ff-dat2plot_fb;
%% get degrees non-zeros per array
%out edges
[nostim_well_out_ff,nostim_well_out_fb,nostim_tunnel_out_ff,nostim_tunnel_out_fb]=...
    get_out_degrees(nostim_slope_G,nostim_datatable,nostim_data_fldr,nostim_temporal_fldr);

[hfs5_well_out_ff,hfs5_well_out_fb,hfs5_tunnel_out_ff,hfs5_tunnel_out_fb]=...
    get_out_degrees(hfs5_slope_G,hfs5_datatable,hfs5_data_fldr,hfs5_temporal_fldr);

[hfs40_well_out_ff,hfs40_well_out_fb,hfs40_tunnel_out_ff,hfs40_tunnel_out_fb]=...
    get_out_degrees(hfs40_slope_G,hfs40_datatable,hfs40_data_fldr,hfs40_temporal_fldr);
% in edges
[nostim_well_in_ff,nostim_well_in_fb,nostim_tunnel_in_ff,nostim_tunnel_in_fb]=...
    get_in_degrees(nostim_slope_G,nostim_datatable,nostim_data_fldr,nostim_temporal_fldr);

[hfs5_well_in_ff,hfs5_well_in_fb,hfs5_tunnel_in_ff,hfs5_tunnel_in_fb]=...
    get_in_degrees(hfs5_slope_G,hfs5_datatable,hfs5_data_fldr,hfs5_temporal_fldr);

[hfs40_well_in_ff,hfs40_well_in_fb,hfs40_tunnel_in_ff,hfs40_tunnel_in_fb]=...
    get_in_degrees(hfs40_slope_G,hfs40_datatable,hfs40_data_fldr,hfs40_temporal_fldr);

for i=1:10
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
    nostim_idx=find(nostim_datatable.s_no==i);
    hfs5_idx=find(hfs5_datatable.s_no==i);
    hfs40_idx=find(hfs40_datatable.s_no==i);
    
    %well in/out
    for regi=1:4
        
        %out
        %ff wells out
        nostim_well_ff_degs=nonzeros(cell2mat(nostim_well_in_ff{nostim_idx,regi}));
        hfs5_well_ff_degs=nonzeros(cell2mat(hfs5_well_in_ff{hfs5_idx,regi}));
        hfs40_well_ff_degs=nonzeros(cell2mat(hfs40_well_in_ff{hfs40_idx,regi}));
        node_degrees_ff_well=[node_degrees_ff_well;...
            [mean(nostim_well_ff_degs),mean(hfs5_well_ff_degs),mean(hfs40_well_ff_degs)]];
        node_degrees_ff_well_se=[node_degrees_ff_well_se;...
            [std(nostim_well_ff_degs)/sqrt(length(nostim_well_ff_degs)),...
            std(hfs5_well_ff_degs)/sqrt(length(hfs5_well_ff_degs)),...
            std(hfs40_well_ff_degs)/sqrt(length(hfs40_well_ff_degs))]];
        
        %fb wells out
        nostim_well_fb_degs=nonzeros(cell2mat(nostim_well_in_fb{nostim_idx,regi}));
        hfs5_well_fb_degs=nonzeros(cell2mat(hfs5_well_in_fb{hfs5_idx,regi}));
        hfs40_well_fb_degs=nonzeros(cell2mat(hfs40_well_in_fb{hfs40_idx,regi}));
        node_degrees_fb_well=[node_degrees_fb_well;...
            [mean(nostim_well_fb_degs),mean(hfs5_well_fb_degs),mean(hfs40_well_fb_degs)]];
        node_degrees_fb_well_se=[node_degrees_fb_well_se;...
            [std(nostim_well_fb_degs)/sqrt(length(nostim_well_fb_degs)),...
            std(hfs5_well_fb_degs)/sqrt(length(hfs5_well_fb_degs)),...
            std(hfs40_well_fb_degs)/sqrt(length(hfs40_well_fb_degs))]];
        
        %in
        %ff in
        nostim_tunnel_ff_degs=nonzeros(cell2mat(nostim_well_out_ff{nostim_idx,regi}));
        hfs5_tunnel_ff_degs=nonzeros(cell2mat(hfs5_well_out_ff{hfs5_idx,regi}));
        hfs40_tunnel_ff_degs=nonzeros(cell2mat(hfs40_well_out_ff{hfs40_idx,regi}));
        node_degrees_ff_tunnel=[node_degrees_ff_tunnel;...
            [mean(nostim_tunnel_ff_degs),mean(hfs5_tunnel_ff_degs),mean(hfs40_tunnel_ff_degs)]];
        node_degrees_ff_tunnel_se=[node_degrees_ff_tunnel_se;...
            [std(nostim_tunnel_ff_degs)/sqrt(length(nostim_tunnel_ff_degs)),...
            std(hfs5_tunnel_ff_degs)/sqrt(length(hfs5_tunnel_ff_degs)),...
            std(hfs40_tunnel_ff_degs)/sqrt(length(hfs40_tunnel_ff_degs))]];
        
        %fb in
        nostim_tunnel_fb_degs=nonzeros(cell2mat(nostim_well_out_fb{nostim_idx,regi}));
        hfs5_tunnel_fb_degs=nonzeros(cell2mat(hfs5_well_out_fb{hfs5_idx,regi}));
        hfs40_tunnel_fb_degs=nonzeros(cell2mat(hfs40_well_out_fb{hfs40_idx,regi}));
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
    
    figure('Position',[0,0,800,800])
    t=tiledlayout(2,2,'TileSpacing','compact');
    ylabel(t,"Degrees Per Node",'FontSize',24)
    xlabel(t,"Subregion",'FontSize',24)
    title(t,"Well Degrees",'FontSize',24)
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
            title("Feedforward Out Degrees",'FontSize',20)
            b=bar(regs_cat,degs{i});
        elseif i==4
            title("Feedback Out Degrees",'FontSize',20)
            b=bar(regs_cat,degs{i});
        elseif i==1
            title("Feedforward In Degrees",'FontSize',20)
            b=bar(regs_cat,degs{i});
        elseif i==2
            title("Feedback In Degrees",'FontSize',20)
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
%     saveas(gcf,"D:\Brewer lab data\HFS\Temporal Analysis\Graph Analysis\degrees well nonzero FID "+string(i)+".png")
    
    %tunnel in/out
    
    % graphing degrees averaged
    node_degrees_ff_well=[];
    node_degrees_ff_well_se=[];
    node_degrees_fb_well=[];
    node_degrees_fb_well_se=[];
    node_degrees_ff_tunnel=[];
    node_degrees_ff_tunnel_se=[];
    node_degrees_fb_tunnel=[];
    node_degrees_fb_tunnel_se=[];
    
    for regi=1:4
        
        %ff out tunnel
        nostim_well_ff_degs=nonzeros(cell2mat(nostim_tunnel_in_ff{nostim_idx,regi}));
        hfs5_well_ff_degs=nonzeros(cell2mat(hfs5_tunnel_in_ff{hfs5_idx,regi}));
        hfs40_well_ff_degs=nonzeros(cell2mat(hfs40_tunnel_in_ff{hfs40_idx,regi}));
        node_degrees_ff_well=[node_degrees_ff_well;...
            [mean(nostim_well_ff_degs),mean(hfs5_well_ff_degs),mean(hfs40_well_ff_degs)]];
        node_degrees_ff_well_se=[node_degrees_ff_well_se;...
            [std(nostim_well_ff_degs)/sqrt(length(nostim_well_ff_degs)),...
            std(hfs5_well_ff_degs)/sqrt(length(hfs5_well_ff_degs)),...
            std(hfs40_well_ff_degs)/sqrt(length(hfs40_well_ff_degs))]];
        
        %fb out tunnel
        nostim_well_fb_degs=nonzeros(cell2mat(nostim_tunnel_in_fb{nostim_idx,regi}));
        hfs5_well_fb_degs=nonzeros(cell2mat(hfs5_tunnel_in_fb{hfs5_idx,regi}));
        hfs40_well_fb_degs=nonzeros(cell2mat(hfs40_tunnel_in_fb{hfs40_idx,regi}));
        node_degrees_fb_well=[node_degrees_fb_well;...
            [mean(nostim_well_fb_degs),mean(hfs5_well_fb_degs),mean(hfs40_well_fb_degs)]];
        node_degrees_fb_well_se=[node_degrees_fb_well_se;...
            [std(nostim_well_fb_degs)/sqrt(length(nostim_well_fb_degs)),...
            std(hfs5_well_fb_degs)/sqrt(length(hfs5_well_fb_degs)),...
            std(hfs40_well_fb_degs)/sqrt(length(hfs40_well_fb_degs))]];
        
        %ff in tunnel
        nostim_tunnel_ff_degs=nonzeros(cell2mat(nostim_tunnel_out_ff{nostim_idx,regi}));
        hfs5_tunnel_ff_degs=nonzeros(cell2mat(hfs5_tunnel_out_ff{hfs5_idx,regi}));
        hfs40_tunnel_ff_degs=nonzeros(cell2mat(hfs40_tunnel_out_ff{hfs40_idx,regi}));
        node_degrees_ff_tunnel=[node_degrees_ff_tunnel;...
            [mean(nostim_tunnel_ff_degs),mean(hfs5_tunnel_ff_degs),mean(hfs40_tunnel_ff_degs)]];
        node_degrees_ff_tunnel_se=[node_degrees_ff_tunnel_se;...
            [std(nostim_tunnel_ff_degs)/sqrt(length(nostim_tunnel_ff_degs)),...
            std(hfs5_tunnel_ff_degs)/sqrt(length(hfs5_tunnel_ff_degs)),...
            std(hfs40_tunnel_ff_degs)/sqrt(length(hfs40_tunnel_ff_degs))]];
        
        %fb in tunnel
        nostim_tunnel_fb_degs=nonzeros(cell2mat(nostim_tunnel_out_fb{nostim_idx,regi}));
        hfs5_tunnel_fb_degs=nonzeros(cell2mat(hfs5_tunnel_out_fb{hfs5_idx,regi}));
        hfs40_tunnel_fb_degs=nonzeros(cell2mat(hfs40_tunnel_out_fb{hfs40_idx,regi}));
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
%         ff_in_tunnel_c{regi}=c;
        
        % fb in anova
        name_mat=[repmat("No Stim",length(nostim_well_fb_degs),1);repmat("5 HFS",length(hfs5_well_fb_degs),1);...
            repmat("40 HFS",length(hfs40_well_fb_degs),1)];
%         [~,~,stats]=anova1([nostim_well_fb_degs;hfs5_well_fb_degs;hfs40_well_fb_degs],name_mat);
%         [c,~]=multcompare(stats);
%         fb_in_tunnel_c{regi}=c;
        
        % ff out anova
        name_mat=[repmat("No Stim",length(nostim_tunnel_ff_degs),1);repmat("5 HFS",length(hfs5_tunnel_ff_degs),1);...
            repmat("40 HFS",length(hfs40_tunnel_ff_degs),1)];
%         [~,~,stats]=anova1([nostim_tunnel_ff_degs;hfs5_tunnel_ff_degs;hfs40_tunnel_ff_degs],name_mat);
%         [c,~]=multcompare(stats);
%         ff_out_tunnel_c{regi}=c;
        
        % fb in anova
        name_mat=[repmat("No Stim",length(nostim_tunnel_fb_degs),1);repmat("5 HFS",length(hfs5_tunnel_fb_degs),1);...
            repmat("40 HFS",length(hfs40_tunnel_fb_degs),1)];
%         [~,~,stats]=anova1([nostim_tunnel_fb_degs;hfs5_tunnel_fb_degs;hfs40_tunnel_fb_degs],name_mat);
%         [c,~]=multcompare(stats);
%         fb_out_tunnel_c{regi}=c;
    end
    
    degs={node_degrees_ff_well,node_degrees_fb_well,node_degrees_ff_tunnel,node_degrees_fb_tunnel};
    degs_se={node_degrees_ff_well_se,node_degrees_fb_well_se,node_degrees_ff_tunnel_se,node_degrees_fb_tunnel_se};
    
    figure('Position',[0,0,800,800])
    t=tiledlayout(2,2,'TileSpacing','compact');
    ylabel(t,"Degrees Per Node",'FontSize',24)
    xlabel(t,"Subregion",'FontSize',24)
    title(t,"Tunnel Degrees",'FontSize',24)
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
            title("Feedforward Out Degrees",'FontSize',20)
            b=bar(regs_tunnel_ff_cat,degs{i});
        elseif i==4
            title("Feedback Out Degrees",'FontSize',20)
            b=bar(regs_tunnel_fb_cat,degs{i});
        elseif i==1
            title("Feedforward In Degrees",'FontSize',20)
            b=bar(regs_tunnel_ff_cat,degs{i});
        elseif i==2
            title("Feedback In Degrees",'FontSize',20)
            b=bar(regs_tunnel_fb_cat,degs{i});
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
%     saveas(gcf,"D:\Brewer lab data\HFS\Temporal Analysis\Graph Analysis\degrees tunnel nonzero FID"+string(i)+".png")
end
%% well get edge counts non-zeros

%out edges
[nostim_well_out_ff,nostim_well_out_fb,nostim_tunnel_out_ff,nostim_tunnel_out_fb]=...
    get_edgecount_out(nostim_slope_G,nostim_datatable,nostim_data_fldr,nostim_temporal_fldr);

[hfs5_well_out_ff,hfs5_well_out_fb,hfs5_tunnel_out_ff,hfs5_tunnel_out_fb]=...
    get_edgecount_out(hfs5_slope_G,hfs5_datatable,hfs5_data_fldr,hfs5_temporal_fldr);

[hfs40_well_out_ff,hfs40_well_out_fb,hfs40_tunnel_out_ff,hfs40_tunnel_out_fb]=...
    get_edgecount_out(hfs40_slope_G,hfs40_datatable,hfs40_data_fldr,hfs40_temporal_fldr);
% in edges
[nostim_well_in_ff,nostim_well_in_fb,nostim_tunnel_in_ff,nostim_tunnel_in_fb]=...
    get_edgecount_in(nostim_slope_G,nostim_datatable,nostim_data_fldr,nostim_temporal_fldr);

[hfs5_well_in_ff,hfs5_well_in_fb,hfs5_tunnel_in_ff,hfs5_tunnel_in_fb]=...
    get_edgecount_in(hfs5_slope_G,hfs5_datatable,hfs5_data_fldr,hfs5_temporal_fldr);

[hfs40_well_in_ff,hfs40_well_in_fb,hfs40_tunnel_in_ff,hfs40_tunnel_in_fb]=...
    get_edgecount_in(hfs40_slope_G,hfs40_datatable,hfs40_data_fldr,hfs40_temporal_fldr);

% graphing degrees averaged
node_degrees_ff_well=[];
node_degrees_ff_well_se=[];
node_degrees_fb_well=[];
node_degrees_fb_well_se=[];
node_degrees_ff_tunnel=[];
node_degrees_ff_tunnel_se=[];
node_degrees_fb_tunnel=[];
node_degrees_fb_tunnel_se=[];

%well in/out
for regi=1:4
    
    %out
    %ff wells out
    nostim_well_ff_degs=nonzeros(cell2mat(nostim_well_in_ff{:,regi}));
    hfs5_well_ff_degs=nonzeros(cell2mat(hfs5_well_in_ff{:,regi}));
    hfs40_well_ff_degs=nonzeros(cell2mat(hfs40_well_in_ff{:,regi}));
    node_degrees_ff_well=[node_degrees_ff_well;...
        [mean(nostim_well_ff_degs),mean(hfs5_well_ff_degs),mean(hfs40_well_ff_degs)]];
    node_degrees_ff_well_se=[node_degrees_ff_well_se;...
        [std(nostim_well_ff_degs)/sqrt(length(nostim_well_ff_degs)),...
        std(hfs5_well_ff_degs)/sqrt(length(hfs5_well_ff_degs)),...
        std(hfs40_well_ff_degs)/sqrt(length(hfs40_well_ff_degs))]];
    
    %fb wells out
    nostim_well_fb_degs=nonzeros(cell2mat(nostim_well_in_fb{:,regi}));
    hfs5_well_fb_degs=nonzeros(cell2mat(hfs5_well_in_fb{:,regi}));
    hfs40_well_fb_degs=nonzeros(cell2mat(hfs40_well_in_fb{:,regi}));
    node_degrees_fb_well=[node_degrees_fb_well;...
        [mean(nostim_well_fb_degs),mean(hfs5_well_fb_degs),mean(hfs40_well_fb_degs)]];
    node_degrees_fb_well_se=[node_degrees_fb_well_se;...
        [std(nostim_well_fb_degs)/sqrt(length(nostim_well_fb_degs)),...
        std(hfs5_well_fb_degs)/sqrt(length(hfs5_well_fb_degs)),...
        std(hfs40_well_fb_degs)/sqrt(length(hfs40_well_fb_degs))]];
    
    %in
    %ff in
    nostim_tunnel_ff_degs=nonzeros(cell2mat(nostim_well_out_ff{:,regi}));
    hfs5_tunnel_ff_degs=nonzeros(cell2mat(hfs5_well_out_ff{:,regi}));
    hfs40_tunnel_ff_degs=nonzeros(cell2mat(hfs40_well_out_ff{:,regi}));
    node_degrees_ff_tunnel=[node_degrees_ff_tunnel;...
        [mean(nostim_tunnel_ff_degs),mean(hfs5_tunnel_ff_degs),mean(hfs40_tunnel_ff_degs)]];
    node_degrees_ff_tunnel_se=[node_degrees_ff_tunnel_se;...
        [std(nostim_tunnel_ff_degs)/sqrt(length(nostim_tunnel_ff_degs)),...
        std(hfs5_tunnel_ff_degs)/sqrt(length(hfs5_tunnel_ff_degs)),...
        std(hfs40_tunnel_ff_degs)/sqrt(length(hfs40_tunnel_ff_degs))]];
    
    %fb in
    nostim_tunnel_fb_degs=nonzeros(cell2mat(nostim_well_out_fb{:,regi}));
    hfs5_tunnel_fb_degs=nonzeros(cell2mat(hfs5_well_out_fb{:,regi}));
    hfs40_tunnel_fb_degs=nonzeros(cell2mat(hfs40_well_out_fb{:,regi}));
    node_degrees_fb_tunnel=[node_degrees_fb_tunnel;...
        [mean(nostim_tunnel_fb_degs),mean(hfs5_tunnel_fb_degs),mean(hfs40_tunnel_fb_degs)]];
    node_degrees_fb_tunnel_se=[node_degrees_fb_tunnel_se;...
        [std(nostim_tunnel_fb_degs)/sqrt(length(nostim_tunnel_fb_degs)),...
        std(hfs5_tunnel_fb_degs)/sqrt(length(hfs5_tunnel_fb_degs)),...
        std(hfs40_tunnel_fb_degs)/sqrt(length(hfs40_tunnel_fb_degs))]];
    
    % ff in anova
    name_mat=[repmat("No Stim",length(nostim_well_ff_degs),1);repmat("5 HFS",length(hfs5_well_ff_degs),1);...
        repmat("40 HFS",length(hfs40_well_ff_degs),1)];
    [~,~,stats]=anova1([nostim_well_ff_degs;hfs5_well_ff_degs;hfs40_well_ff_degs],name_mat);
    [c,~]=multcompare(stats);
    ff_in_well_c{regi}=c;
    
    % fb in anova
    name_mat=[repmat("No Stim",length(nostim_well_fb_degs),1);repmat("5 HFS",length(hfs5_well_fb_degs),1);...
        repmat("40 HFS",length(hfs40_well_fb_degs),1)];
    [~,~,stats]=anova1([nostim_well_fb_degs;hfs5_well_fb_degs;hfs40_well_fb_degs],name_mat);
    [c,~]=multcompare(stats);
    fb_in_well_c{regi}=c;
    
    % ff out anova
    name_mat=[repmat("No Stim",length(nostim_tunnel_ff_degs),1);repmat("5 HFS",length(hfs5_tunnel_ff_degs),1);...
        repmat("40 HFS",length(hfs40_tunnel_ff_degs),1)];
    [~,~,stats]=anova1([nostim_tunnel_ff_degs;hfs5_tunnel_ff_degs;hfs40_tunnel_ff_degs],name_mat);
    [c,~]=multcompare(stats);
    ff_out_well_c{regi}=c;
    
    % fb in anova
    name_mat=[repmat("No Stim",length(nostim_tunnel_fb_degs),1);repmat("5 HFS",length(hfs5_tunnel_fb_degs),1);...
        repmat("40 HFS",length(hfs40_tunnel_fb_degs),1)];
    [~,~,stats]=anova1([nostim_tunnel_fb_degs;hfs5_tunnel_fb_degs;hfs40_tunnel_fb_degs],name_mat);
    [c,~]=multcompare(stats);
    fb_out_well_c{regi}=c;
end

degs={node_degrees_ff_well,node_degrees_fb_well,node_degrees_ff_tunnel,node_degrees_fb_tunnel};
degs_se={node_degrees_ff_well_se,node_degrees_fb_well_se,node_degrees_ff_tunnel_se,node_degrees_fb_tunnel_se};

figure('Position',[0,0,800,800])
t=tiledlayout(2,2,'TileSpacing','compact');
ylabel(t,"Edges Per Node",'FontSize',24)
xlabel(t,"Subregion",'FontSize',24)
title(t,"Well Edge Count",'FontSize',24)
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
        title("Feedforward Out Edges",'FontSize',20)
        b=bar(regs_cat,degs{i});
    elseif i==4
        title("Feedback Out Edges",'FontSize',20)
        b=bar(regs_cat,degs{i});
    elseif i==1
        title("Feedforward In Edges",'FontSize',20)
        b=bar(regs_cat,degs{i});
    elseif i==2
        title("Feedback In Edges",'FontSize',20)
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
% saveas(gcf,"D:\Brewer lab data\HFS\Temporal Analysis\Graph Analysis\edgecount well nonzero.png")

%% tunnel in/out

% graphing degrees averaged
node_degrees_ff_well=[];
node_degrees_ff_well_se=[];
node_degrees_fb_well=[];
node_degrees_fb_well_se=[];
node_degrees_ff_tunnel=[];
node_degrees_ff_tunnel_se=[];
node_degrees_fb_tunnel=[];
node_degrees_fb_tunnel_se=[];

for regi=1:4
    
    %ff out tunnel
    nostim_well_ff_degs=nonzeros(cell2mat(nostim_tunnel_in_ff{:,regi}));
    hfs5_well_ff_degs=nonzeros(cell2mat(hfs5_tunnel_in_ff{:,regi}));
    hfs40_well_ff_degs=nonzeros(cell2mat(hfs40_tunnel_in_ff{:,regi}));
    node_degrees_ff_well=[node_degrees_ff_well;...
        [mean(nostim_well_ff_degs),mean(hfs5_well_ff_degs),mean(hfs40_well_ff_degs)]];
    node_degrees_ff_well_se=[node_degrees_ff_well_se;...
        [std(nostim_well_ff_degs)/sqrt(length(nostim_well_ff_degs)),...
        std(hfs5_well_ff_degs)/sqrt(length(hfs5_well_ff_degs)),...
        std(hfs40_well_ff_degs)/sqrt(length(hfs40_well_ff_degs))]];
    
    %fb out tunnel
    nostim_well_fb_degs=nonzeros(cell2mat(nostim_tunnel_in_fb{:,regi}));
    hfs5_well_fb_degs=nonzeros(cell2mat(hfs5_tunnel_in_fb{:,regi}));
    hfs40_well_fb_degs=nonzeros(cell2mat(hfs40_tunnel_in_fb{:,regi}));
    node_degrees_fb_well=[node_degrees_fb_well;...
        [mean(nostim_well_fb_degs),mean(hfs5_well_fb_degs),mean(hfs40_well_fb_degs)]];
    node_degrees_fb_well_se=[node_degrees_fb_well_se;...
        [std(nostim_well_fb_degs)/sqrt(length(nostim_well_fb_degs)),...
        std(hfs5_well_fb_degs)/sqrt(length(hfs5_well_fb_degs)),...
        std(hfs40_well_fb_degs)/sqrt(length(hfs40_well_fb_degs))]];
    
    %ff in tunnel
    nostim_tunnel_ff_degs=nonzeros(cell2mat(nostim_tunnel_out_ff{:,regi}));
    hfs5_tunnel_ff_degs=nonzeros(cell2mat(hfs5_tunnel_out_ff{:,regi}));
    hfs40_tunnel_ff_degs=nonzeros(cell2mat(hfs40_tunnel_out_ff{:,regi}));
    node_degrees_ff_tunnel=[node_degrees_ff_tunnel;...
        [mean(nostim_tunnel_ff_degs),mean(hfs5_tunnel_ff_degs),mean(hfs40_tunnel_ff_degs)]];
    node_degrees_ff_tunnel_se=[node_degrees_ff_tunnel_se;...
        [std(nostim_tunnel_ff_degs)/sqrt(length(nostim_tunnel_ff_degs)),...
        std(hfs5_tunnel_ff_degs)/sqrt(length(hfs5_tunnel_ff_degs)),...
        std(hfs40_tunnel_ff_degs)/sqrt(length(hfs40_tunnel_ff_degs))]];
    
    %fb in tunnel
    nostim_tunnel_fb_degs=nonzeros(cell2mat(nostim_tunnel_out_fb{:,regi}));
    hfs5_tunnel_fb_degs=nonzeros(cell2mat(hfs5_tunnel_out_fb{:,regi}));
    hfs40_tunnel_fb_degs=nonzeros(cell2mat(hfs40_tunnel_out_fb{:,regi}));
    node_degrees_fb_tunnel=[node_degrees_fb_tunnel;...
        [mean(nostim_tunnel_fb_degs),mean(hfs5_tunnel_fb_degs),mean(hfs40_tunnel_fb_degs)]];
    node_degrees_fb_tunnel_se=[node_degrees_fb_tunnel_se;...
        [std(nostim_tunnel_fb_degs)/sqrt(length(nostim_tunnel_fb_degs)),...
        std(hfs5_tunnel_fb_degs)/sqrt(length(hfs5_tunnel_fb_degs)),...
        std(hfs40_tunnel_fb_degs)/sqrt(length(hfs40_tunnel_fb_degs))]];
    
    % ff in anova
    name_mat=[repmat("No Stim",length(nostim_well_ff_degs),1);repmat("5 HFS",length(hfs5_well_ff_degs),1);...
        repmat("40 HFS",length(hfs40_well_ff_degs),1)];
    [~,~,stats]=anova1([nostim_well_ff_degs;hfs5_well_ff_degs;hfs40_well_ff_degs],name_mat);
    [c,~]=multcompare(stats);
    ff_in_tunnel_c{regi}=c;
    
    % fb in anova
    name_mat=[repmat("No Stim",length(nostim_well_fb_degs),1);repmat("5 HFS",length(hfs5_well_fb_degs),1);...
        repmat("40 HFS",length(hfs40_well_fb_degs),1)];
    [~,~,stats]=anova1([nostim_well_fb_degs;hfs5_well_fb_degs;hfs40_well_fb_degs],name_mat);
    [c,~]=multcompare(stats);
    fb_in_tunnel_c{regi}=c;
    
    % ff out anova
    name_mat=[repmat("No Stim",length(nostim_tunnel_ff_degs),1);repmat("5 HFS",length(hfs5_tunnel_ff_degs),1);...
        repmat("40 HFS",length(hfs40_tunnel_ff_degs),1)];
    [~,~,stats]=anova1([nostim_tunnel_ff_degs;hfs5_tunnel_ff_degs;hfs40_tunnel_ff_degs],name_mat);
    [c,~]=multcompare(stats);
    ff_out_tunnel_c{regi}=c;
    
    % fb in anova
    name_mat=[repmat("No Stim",length(nostim_tunnel_fb_degs),1);repmat("5 HFS",length(hfs5_tunnel_fb_degs),1);...
        repmat("40 HFS",length(hfs40_tunnel_fb_degs),1)];
    [~,~,stats]=anova1([nostim_tunnel_fb_degs;hfs5_tunnel_fb_degs;hfs40_tunnel_fb_degs],name_mat);
    [c,~]=multcompare(stats);
    fb_out_tunnel_c{regi}=c;
end

degs={node_degrees_ff_well,node_degrees_fb_well,node_degrees_ff_tunnel,node_degrees_fb_tunnel};
degs_se={node_degrees_ff_well_se,node_degrees_fb_well_se,node_degrees_ff_tunnel_se,node_degrees_fb_tunnel_se};

figure('Position',[0,0,800,800])
t=tiledlayout(2,2,'TileSpacing','compact');
ylabel(t,"Edges Per Node",'FontSize',24)
xlabel(t,"Subregion",'FontSize',24)
title(t,"Tunnel Edge Count",'FontSize',24)
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
        title("Feedforward Out Edges",'FontSize',20)
        b=bar(regs_tunnel_ff_cat,degs{i});
    elseif i==4
        title("Feedback Out Edges",'FontSize',20)
        b=bar(regs_tunnel_fb_cat,degs{i});
    elseif i==1
        title("Feedforward In Edges",'FontSize',20)
        b=bar(regs_tunnel_ff_cat,degs{i});
    elseif i==2
        title("Feedback In Edges",'FontSize',20)
        b=bar(regs_tunnel_fb_cat,degs{i});
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
% saveas(gcf,"D:\Brewer lab data\HFS\Temporal Analysis\Graph Analysis\edgecount tunnel nonzero.png")
%% get centrality
[nostim_well_in_ff,nostim_well_in_fb,~,~]=...
    get_centrality_in(nostim_slope_G,nostim_datatable,nostim_data_fldr,nostim_temporal_fldr);

[hfs5_well_in_ff,hfs5_well_in_fb,~,~]=...
    get_centrality_in(hfs5_slope_G,hfs5_datatable,hfs5_data_fldr,hfs5_temporal_fldr);

[hfs40_well_in_ff,hfs40_well_in_fb,~,~]=...
    get_centrality_in(hfs40_slope_G,hfs40_datatable,hfs40_data_fldr,hfs40_temporal_fldr);

[nostim_well_out_ff,nostim_well_out_fb,~,~]=...
    get_centrality_out(nostim_slope_G,nostim_datatable,nostim_data_fldr,nostim_temporal_fldr);

[hfs5_well_out_ff,hfs5_well_out_fb,~,~]=...
    get_centrality_out(hfs5_slope_G,hfs5_datatable,hfs5_data_fldr,hfs5_temporal_fldr);

[hfs40_well_out_ff,hfs40_well_out_fb,~,~]=...
    get_centrality_out(hfs40_slope_G,hfs40_datatable,hfs40_data_fldr,hfs40_temporal_fldr);

% graphing degrees averaged
node_in_ff=[];
node_in_ff_se=[];
node_in_fb=[];
node_in_fb_se=[];
node_out_ff=[];
node_out_ff_se=[];
node_out_fb=[];
node_out_fb_se=[];

for regi=1:4
    
    %ff in
    nostim_well_ff_degs=cell2mat(nostim_well_in_ff{:,regi});
    hfs5_well_ff_degs=cell2mat(hfs5_well_in_ff{:,regi});
    hfs40_well_ff_degs=cell2mat(hfs40_well_in_ff{:,regi});
    node_in_ff=[node_in_ff;...
        [mean(nostim_well_ff_degs),mean(hfs5_well_ff_degs),mean(hfs40_well_ff_degs)]];
    node_in_ff_se=[node_in_ff_se;...
        [std(nostim_well_ff_degs)/sqrt(length(nostim_well_ff_degs)),...
        std(hfs5_well_ff_degs)/sqrt(length(hfs5_well_ff_degs)),...
        std(hfs40_well_ff_degs)/sqrt(length(hfs40_well_ff_degs))]];
    
    %fb in
    nostim_well_fb_degs=cell2mat(nostim_well_in_fb{:,regi});
    hfs5_well_fb_degs=cell2mat(hfs5_well_in_fb{:,regi});
    hfs40_well_fb_degs=cell2mat(hfs40_well_in_fb{:,regi});
    node_in_fb=[node_in_fb;...
        [mean(nostim_well_fb_degs),mean(hfs5_well_fb_degs),mean(hfs40_well_fb_degs)]];
    node_in_fb_se=[node_in_fb_se;...
        [std(nostim_well_fb_degs)/sqrt(length(nostim_well_fb_degs)),...
        std(hfs5_well_fb_degs)/sqrt(length(hfs5_well_fb_degs)),...
        std(hfs40_well_fb_degs)/sqrt(length(hfs40_well_fb_degs))]];
    
    %ff out
    nostim_tunnel_ff_degs=cell2mat(nostim_well_out_ff{:,regi});
    hfs5_tunnel_ff_degs=cell2mat(hfs5_well_out_ff{:,regi});
    hfs40_tunnel_ff_degs=cell2mat(hfs40_well_out_ff{:,regi});
    node_out_ff=[node_out_ff;...
        [mean(nostim_tunnel_ff_degs),mean(hfs5_tunnel_ff_degs),mean(hfs40_tunnel_ff_degs)]];
    node_out_ff_se=[node_out_ff_se;...
        [std(nostim_tunnel_ff_degs)/sqrt(length(nostim_tunnel_ff_degs)),...
        std(hfs5_tunnel_ff_degs)/sqrt(length(hfs5_tunnel_ff_degs)),...
        std(hfs40_tunnel_ff_degs)/sqrt(length(hfs40_tunnel_ff_degs))]];
    
    %fb out
    nostim_tunnel_fb_degs=cell2mat(nostim_well_out_fb{:,regi});
    hfs5_tunnel_fb_degs=cell2mat(hfs5_well_out_fb{:,regi});
    hfs40_tunnel_fb_degs=cell2mat(hfs40_well_out_fb{:,regi});
    node_out_fb=[node_out_fb;...
        [mean(nostim_tunnel_fb_degs),mean(hfs5_tunnel_fb_degs),mean(hfs40_tunnel_fb_degs)]];
    node_out_fb_se=[node_out_fb_se;...
        [std(nostim_tunnel_fb_degs)/sqrt(length(nostim_tunnel_fb_degs)),...
        std(hfs5_tunnel_fb_degs)/sqrt(length(hfs5_tunnel_fb_degs)),...
        std(hfs40_tunnel_fb_degs)/sqrt(length(hfs40_tunnel_fb_degs))]];
end

degs={node_in_ff,node_in_fb,node_out_ff,node_out_fb};
degs_se={node_in_ff_se,node_in_ff_se,node_out_ff_se,node_out_fb_se};

figure('Position',[0,0,800,800])
t=tiledlayout(2,2,'TileSpacing','compact');
ylabel(t,"Centrality")
xlabel(t,"Subregion")
title(t,"Centrality")
for i=1:4
    ax(i)=nexttile;
    hold on
    regs=["EC","DG","CA3","CA1"];
    regs_cat=categorical(regs);
    regs_cat=reordercats(regs_cat,regs);
%     regs_tunnel_ff=["EC-DG","DG-CA3","CA3-CA1","CA1-EC"];
%     regs_tunnel_fb=["DG-EC","CA3-DG","CA1-CA3","EC-CA1"];
    
    regs_tunnel_ff_cat=reordercats(categorical(regs_tunnel_ff),regs_tunnel_ff);
    regs_tunnel_fb_cat=reordercats(categorical(regs_tunnel_fb),regs_tunnel_fb);
    
    if i==3
        title("Feedforward Out Degree Centrality")
        b=bar(regs_cat,degs{i});
    elseif i==4
        title("Feedback Out Degree Centrality")
        b=bar(regs_cat,degs{i});
    elseif i==1
        title("Feedforward In Degree Centrality")
        b=bar(regs_cat,degs{i});
    elseif i==2
        title("Feedback In Degree Centrality")
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
% saveas(gcf,"D:\Brewer lab data\HFS\Temporal Analysis\Graph Analysis\centrality.png")

%% get centrality nonzero slope
%out edges
[nostim_well_out_ff,nostim_well_out_fb,nostim_tunnel_out_ff,nostim_tunnel_out_fb]=...
    get_centrality_out(nostim_slope_G,nostim_datatable,nostim_data_fldr,nostim_temporal_fldr);

[hfs5_well_out_ff,hfs5_well_out_fb,hfs5_tunnel_out_ff,hfs5_tunnel_out_fb]=...
    get_centrality_out(hfs5_slope_G,hfs5_datatable,hfs5_data_fldr,hfs5_temporal_fldr);

[hfs40_well_out_ff,hfs40_well_out_fb,hfs40_tunnel_out_ff,hfs40_tunnel_out_fb]=...
    get_centrality_out(hfs40_slope_G,hfs40_datatable,hfs40_data_fldr,hfs40_temporal_fldr);
% in edges
[nostim_well_in_ff,nostim_well_in_fb,nostim_tunnel_in_ff,nostim_tunnel_in_fb]=...
    get_centrality_in(nostim_slope_G,nostim_datatable,nostim_data_fldr,nostim_temporal_fldr);

[hfs5_well_in_ff,hfs5_well_in_fb,hfs5_tunnel_in_ff,hfs5_tunnel_in_fb]=...
    get_centrality_in(hfs5_slope_G,hfs5_datatable,hfs5_data_fldr,hfs5_temporal_fldr);

[hfs40_well_in_ff,hfs40_well_in_fb,hfs40_tunnel_in_ff,hfs40_tunnel_in_fb]=...
    get_centrality_in(hfs40_slope_G,hfs40_datatable,hfs40_data_fldr,hfs40_temporal_fldr);

% graphing degrees averaged
node_degrees_ff_well=[];
node_degrees_ff_well_se=[];
node_degrees_fb_well=[];
node_degrees_fb_well_se=[];
node_degrees_ff_tunnel=[];
node_degrees_ff_tunnel_se=[];
node_degrees_fb_tunnel=[];
node_degrees_fb_tunnel_se=[];

%% well in/out
ff_in_well_m=[];
ff_in_well_se=[];
fb_in_well_m=[];
fb_in_well_se=[];
ff_out_well_m=[];
ff_out_well_se=[];
fb_out_well_m=[];
fb_out_well_se=[];

for regi=1:4
    
    %out
    %ff wells out
    nostim_well_ff_degs=nonzeros(cell2mat(nostim_well_in_ff{:,regi}));
    hfs5_well_ff_degs=nonzeros(cell2mat(hfs5_well_in_ff{:,regi}));
    hfs40_well_ff_degs=nonzeros(cell2mat(hfs40_well_in_ff{:,regi}));
    node_degrees_ff_well=[node_degrees_ff_well;...
        [mean(nostim_well_ff_degs),mean(hfs5_well_ff_degs),mean(hfs40_well_ff_degs)]];
    node_degrees_ff_well_se=[node_degrees_ff_well_se;...
        [std(nostim_well_ff_degs)/sqrt(length(nostim_well_ff_degs)),...
        std(hfs5_well_ff_degs)/sqrt(length(hfs5_well_ff_degs)),...
        std(hfs40_well_ff_degs)/sqrt(length(hfs40_well_ff_degs))]];
    
    %fb wells out
    nostim_well_fb_degs=nonzeros(cell2mat(nostim_well_in_fb{:,regi}));
    hfs5_well_fb_degs=nonzeros(cell2mat(hfs5_well_in_fb{:,regi}));
    hfs40_well_fb_degs=nonzeros(cell2mat(hfs40_well_in_fb{:,regi}));
    node_degrees_fb_well=[node_degrees_fb_well;...
        [mean(nostim_well_fb_degs),mean(hfs5_well_fb_degs),mean(hfs40_well_fb_degs)]];
    node_degrees_fb_well_se=[node_degrees_fb_well_se;...
        [std(nostim_well_fb_degs)/sqrt(length(nostim_well_fb_degs)),...
        std(hfs5_well_fb_degs)/sqrt(length(hfs5_well_fb_degs)),...
        std(hfs40_well_fb_degs)/sqrt(length(hfs40_well_fb_degs))]];
    
    %in
    %ff in
    nostim_tunnel_ff_degs=nonzeros(cell2mat(nostim_well_out_ff{:,regi}));
    hfs5_tunnel_ff_degs=nonzeros(cell2mat(hfs5_well_out_ff{:,regi}));
    hfs40_tunnel_ff_degs=nonzeros(cell2mat(hfs40_well_out_ff{:,regi}));
    node_degrees_ff_tunnel=[node_degrees_ff_tunnel;...
        [mean(nostim_tunnel_ff_degs),mean(hfs5_tunnel_ff_degs),mean(hfs40_tunnel_ff_degs)]];
    node_degrees_ff_tunnel_se=[node_degrees_ff_tunnel_se;...
        [std(nostim_tunnel_ff_degs)/sqrt(length(nostim_tunnel_ff_degs)),...
        std(hfs5_tunnel_ff_degs)/sqrt(length(hfs5_tunnel_ff_degs)),...
        std(hfs40_tunnel_ff_degs)/sqrt(length(hfs40_tunnel_ff_degs))]];
    
    %fb in
    nostim_tunnel_fb_degs=nonzeros(cell2mat(nostim_well_out_fb{:,regi}));
    hfs5_tunnel_fb_degs=nonzeros(cell2mat(hfs5_well_out_fb{:,regi}));
    hfs40_tunnel_fb_degs=nonzeros(cell2mat(hfs40_well_out_fb{:,regi}));
    node_degrees_fb_tunnel=[node_degrees_fb_tunnel;...
        [mean(nostim_tunnel_fb_degs),mean(hfs5_tunnel_fb_degs),mean(hfs40_tunnel_fb_degs)]];
    node_degrees_fb_tunnel_se=[node_degrees_fb_tunnel_se;...
        [std(nostim_tunnel_fb_degs)/sqrt(length(nostim_tunnel_fb_degs)),...
        std(hfs5_tunnel_fb_degs)/sqrt(length(hfs5_tunnel_fb_degs)),...
        std(hfs40_tunnel_fb_degs)/sqrt(length(hfs40_tunnel_fb_degs))]];
    
    % ff in anova
    name_mat=[repmat("No Stim",length(nostim_well_ff_degs),1);repmat("5 HFS",length(hfs5_well_ff_degs),1);...
        repmat("40 HFS",length(hfs40_well_ff_degs),1)];
    [~,~,stats]=anova1([nostim_well_ff_degs;hfs5_well_ff_degs;hfs40_well_ff_degs],name_mat);
    [c,m]=multcompare(stats);
    ff_in_well_c{regi}=c;
    ff_in_well_m=[ff_in_well_m;m(:,1)'];
    ff_in_well_se=[ff_in_well_se;m(:,2)'];
    
    % fb in anova
    name_mat=[repmat("No Stim",length(nostim_well_fb_degs),1);repmat("5 HFS",length(hfs5_well_fb_degs),1);...
        repmat("40 HFS",length(hfs40_well_fb_degs),1)];
    [~,~,stats]=anova1([nostim_well_fb_degs;hfs5_well_fb_degs;hfs40_well_fb_degs],name_mat);
    [c,m]=multcompare(stats);
    fb_in_well_c{regi}=c;
    fb_in_well_m=[fb_in_well_m;m(:,1)'];
    fb_in_well_se=[fb_in_well_se;m(:,2)'];
    
    % ff out anova
    name_mat=[repmat("No Stim",length(nostim_tunnel_ff_degs),1);repmat("5 HFS",length(hfs5_tunnel_ff_degs),1);...
        repmat("40 HFS",length(hfs40_tunnel_ff_degs),1)];
    [~,~,stats]=anova1([nostim_tunnel_ff_degs;hfs5_tunnel_ff_degs;hfs40_tunnel_ff_degs],name_mat);
    [c,m]=multcompare(stats);
    ff_out_well_c{regi}=c;
    ff_out_well_m=[ff_out_well_m;m(:,1)'];
    ff_out_well_se=[ff_out_well_se;m(:,2)'];
    
    % fb in anova
    name_mat=[repmat("No Stim",length(nostim_tunnel_fb_degs),1);repmat("5 HFS",length(hfs5_tunnel_fb_degs),1);...
        repmat("40 HFS",length(hfs40_tunnel_fb_degs),1)];
    [~,~,stats]=anova1([nostim_tunnel_fb_degs;hfs5_tunnel_fb_degs;hfs40_tunnel_fb_degs],name_mat);
    [c,m]=multcompare(stats);
    fb_out_well_c{regi}=c;
    fb_out_well_m=[fb_out_well_m;m(:,1)'];
    fb_out_well_se=[fb_out_well_se;m(:,2)'];
    
    %in ff fb anova
    v=[nostim_well_ff_degs;nostim_well_fb_degs];
    g=[repmat("ff",length(nostim_well_ff_degs),1);repmat("fb",length(nostim_well_fb_degs),1)];
    [~,~,stats]=anova1(v,g);
    [c,~]=multcompare(stats);
    ff_fb_in_c{regi,1}=c;
    
    v=[hfs5_well_ff_degs;hfs5_well_fb_degs];
    g=[repmat("ff",length(hfs5_well_ff_degs),1);repmat("fb",length(hfs5_well_fb_degs),1)];
    [~,~,stats]=anova1(v,g);
    [c,~]=multcompare(stats);
    ff_fb_in_c{regi,2}=c;
    
    v=[hfs40_well_ff_degs;hfs40_well_fb_degs];
    g=[repmat("ff",length(hfs40_well_ff_degs),1);repmat("fb",length(hfs40_well_fb_degs),1)];
    [~,~,stats]=anova1(v,g);
    [c,~]=multcompare(stats);
    ff_fb_in_c{regi,3}=c;
    
    %out ff fb anova
    v=[nostim_tunnel_ff_degs;nostim_tunnel_fb_degs];
    g=[repmat("ff",length(nostim_tunnel_ff_degs),1);repmat("fb",length(nostim_tunnel_fb_degs),1)];
    [~,~,stats]=anova1(v,g);
    [c,~]=multcompare(stats);
    ff_fb_out_c{regi,1}=c;
    
    v=[hfs5_tunnel_ff_degs;hfs5_tunnel_fb_degs];
    g=[repmat("ff",length(hfs5_tunnel_ff_degs),1);repmat("fb",length(hfs5_tunnel_fb_degs),1)];
    [~,~,stats]=anova1(v,g);
    [c,~]=multcompare(stats);
    ff_fb_out_c{regi,2}=c;
    
    v=[hfs40_tunnel_ff_degs;hfs40_tunnel_fb_degs];
    g=[repmat("ff",length(hfs40_tunnel_ff_degs),1);repmat("fb",length(hfs40_tunnel_fb_degs),1)];
    [~,~,stats]=anova1(v,g);
    [c,~]=multcompare(stats);
    ff_fb_out_c{regi,3}=c;

    effect_in{regi}=ff_fb_meanEffect(nostim_well_ff_degs,nostim_well_fb_degs,hfs5_well_ff_degs,hfs5_well_fb_degs,hfs40_well_ff_degs,hfs40_well_fb_degs);
    effect_out{regi}=ff_fb_meanEffect(nostim_tunnel_ff_degs,nostim_tunnel_fb_degs,hfs5_tunnel_ff_degs,hfs5_tunnel_fb_degs,hfs40_tunnel_ff_degs,hfs40_tunnel_fb_degs);
end

ff_fb_c = reshape([ff_fb_in_c(:) ff_fb_out_c(:)]',2*size(ff_fb_in_c,1), []);
get_last=@(x) x(1,end);
ff_fb_p = cellfun(get_last,ff_fb_c);

% degs={node_degrees_ff_well,node_degrees_fb_well,node_degrees_ff_tunnel,node_degrees_fb_tunnel};
% degs_se={node_degrees_ff_well_se,node_degrees_fb_well_se,node_degrees_ff_tunnel_se,node_degrees_fb_tunnel_se};

degs={ff_in_well_m,fb_in_well_m,ff_out_well_m,fb_out_well_m};
degs_se={ff_in_well_se,fb_in_well_se,ff_out_well_se,fb_out_well_se};

figure('Position',[0,0,800,800])
t=tiledlayout(2,2,'TileSpacing','compact');
ylabel(t,"Centrality Per Node",'FontSize',24)
xlabel(t,"Subregion",'FontSize',24)
title(t,"Well Slope Centrality",'FontSize',24)
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
        title("Feedforward Out Centrality",'FontSize',20)
        b=bar(regs_cat,degs{i});
    elseif i==4
        title("Feedback Out Centrality",'FontSize',20)
        b=bar(regs_cat,degs{i});
    elseif i==1
        title("Feedforward In Centrality",'FontSize',20)
        b=bar(regs_cat,degs{i});
    elseif i==2
        title("Feedback In Centrality",'FontSize',20)
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
% saveas(gcf,"D:\Brewer lab data\HFS\Temporal Analysis\Graph Analysis\centrality well nonzero.png")
%% horzbar from anova
dat2plot_ff=[];
dat2plot_ff_se=[];
dat2plot_fb=[];
dat2plot_fb_se=[];
labels={'EC In','EC Out','DG In','DG Out','CA3 In','CA3 Out','CA1 In','CA1 Out'};
figure( 'Position', [100 100 1400 600])
for i=1:4
    dat2plot_ff=[dat2plot_ff;degs{1}(i,:);degs{3}(i,:)];
    dat2plot_ff_se=[dat2plot_ff_se;degs_se{1}(i,:);degs_se{3}(i,:)];
    dat2plot_fb=[dat2plot_fb;degs{2}(i,:);degs{4}(i,:)];
    dat2plot_fb_se=[dat2plot_fb_se;degs_se{2}(i,:);degs_se{4}(i,:)];
end

plot_horzbar_from_anova...
    (dat2plot_ff,dat2plot_ff_se,dat2plot_fb,dat2plot_fb_se,labels,ff_colors,fb_colors)
xlabel("FB Slope Centrality                  FF Slope Centrality")
xlim([-5,5])
xticks([-5:1:5])
ylim([0,9])
set(gca,'LineWidth',2,'TickLength',[0.02 0.02]);
set(gca,'FontSize',18)
% exportgraphics(gcf,".\three stim share figs\slope_well_centrality_comp.png",'Resolution',1500)
exportgraphics(gcf,".\Temporal Analysis\Graph Analysis\slope_well_centrality_comp.png",'Resolution',1500)

mean_diff=dat2plot_ff-dat2plot_fb;

%% FF/FB% slope centrality

ff_fb_perc=(dat2plot_ff./dat2plot_fb)*100;

% ff_fb_SD=dat2plot

ff_fb_SD=[];
ff_fb_lengths=[];
for regi=1:4
    
    %out
    %ff wells out
    nostim_well_ff_degs=nonzeros(cell2mat(nostim_well_in_ff{:,regi}));
    hfs5_well_ff_degs=nonzeros(cell2mat(hfs5_well_in_ff{:,regi}));
    hfs40_well_ff_degs=nonzeros(cell2mat(hfs40_well_in_ff{:,regi}));
    node_degrees_ff_well_m=[mean(nostim_well_ff_degs),mean(hfs5_well_ff_degs),mean(hfs40_well_ff_degs)];
    node_degrees_ff_well_sd=[std(nostim_well_ff_degs),std(hfs5_well_ff_degs),std(hfs40_well_ff_degs)];
    ff_out_lengths=[length(nostim_well_ff_degs),length(hfs5_well_ff_degs),length(hfs40_well_ff_degs)];

    %fb wells out
    nostim_well_fb_degs=nonzeros(cell2mat(nostim_well_in_fb{:,regi}));
    hfs5_well_fb_degs=nonzeros(cell2mat(hfs5_well_in_fb{:,regi}));
    hfs40_well_fb_degs=nonzeros(cell2mat(hfs40_well_in_fb{:,regi}));
    node_degrees_fb_well_m=[mean(nostim_well_fb_degs),mean(hfs5_well_fb_degs),mean(hfs40_well_fb_degs)];
    node_degrees_fb_well_sd=[std(nostim_well_fb_degs),std(hfs5_well_fb_degs),std(hfs40_well_fb_degs)];
    fb_out_lengths=[length(nostim_well_fb_degs),length(hfs5_well_fb_degs),length(hfs40_well_fb_degs)];
    
    %in
    %ff in
    nostim_tunnel_ff_degs=nonzeros(cell2mat(nostim_well_out_ff{:,regi}));
    hfs5_tunnel_ff_degs=nonzeros(cell2mat(hfs5_well_out_ff{:,regi}));
    hfs40_tunnel_ff_degs=nonzeros(cell2mat(hfs40_well_out_ff{:,regi}));
    node_degrees_ff_tunnel_m=[mean(nostim_tunnel_ff_degs),mean(hfs5_tunnel_ff_degs),mean(hfs40_tunnel_ff_degs)];
    node_degrees_ff_tunnel_sd=[std(nostim_tunnel_ff_degs),std(hfs5_tunnel_ff_degs), std(hfs40_tunnel_ff_degs)];
    ff_in_lengths=[length(nostim_tunnel_ff_degs),length(hfs5_tunnel_ff_degs),length(hfs40_tunnel_ff_degs)];

    %fb in
    nostim_tunnel_fb_degs=nonzeros(cell2mat(nostim_well_out_fb{:,regi}));
    hfs5_tunnel_fb_degs=nonzeros(cell2mat(hfs5_well_out_fb{:,regi}));
    hfs40_tunnel_fb_degs=nonzeros(cell2mat(hfs40_well_out_fb{:,regi}));
    node_degrees_fb_tunnel_m=[mean(nostim_tunnel_fb_degs),mean(hfs5_tunnel_fb_degs),mean(hfs40_tunnel_fb_degs)];
    node_degrees_fb_tunnel_sd=[std(nostim_tunnel_fb_degs),std(hfs5_tunnel_fb_degs),std(hfs40_tunnel_fb_degs)];
    fb_in_lengths=[length(nostim_tunnel_fb_degs),length(hfs5_tunnel_fb_degs),length(hfs40_tunnel_fb_degs)];

    ff_fb_SD=[ff_fb_SD;((node_degrees_ff_well_m./node_degrees_fb_well_m)*100).*...
        sqrt((node_degrees_ff_well_sd./node_degrees_fb_well_sd).^2+...
        (node_degrees_fb_well_sd./node_degrees_fb_well_sd).^2);...
        ((node_degrees_ff_tunnel_m./node_degrees_fb_tunnel_m)*100).*...
        sqrt((node_degrees_ff_tunnel_m./node_degrees_fb_tunnel_m).^2+...
        (node_degrees_fb_tunnel_sd./node_degrees_fb_tunnel_sd).^2)];
    ff_fb_lengths=[ff_fb_lengths;(ff_in_lengths+fb_in_lengths);(ff_out_lengths+fb_out_lengths)];
end

labels={'EC-DG In','EC-DG Out','DG-CA3 In','DG-CA3 Out','CA3-CA1 In','CA3-CA1 Out','CA1-EC In','CA1-EC Out'};
regs=categorical(labels);
regs=reordercats(regs,labels);

figure( 'Position', [100 100 700 600])
b=bar(regs,ff_fb_perc);
b(1).BaseValue=100;
b(1).BaseLine.LineStyle = "--";
b(1).BaseLine.LineWidth=2;

ff_fb_SE=ff_fb_SD./sqrt(ff_fb_lengths);
hold on
% grouped error bars
[ngroups,nbars]=size(ff_fb_perc);
x=[];
for j=1:nbars
    x(j,:)=b(j).XEndPoints;
end
% errorbar(x',ff_fb_perc,ff_fb_SE,'k','LineStyle','none')
hold off
ylabel("% FF/FB")
set(gca,'LineWidth',2,'TickLength',[0.02 0.02]);
set(gca,'FontSize',18)
exportgraphics(gcf,".\Temporal Analysis\Graph Analysis\slope_centrality_ff_fb_perc.png",'Resolution',1500)

%% Cohen's D Bar Graphs
CD_in=[];
CD_out=[];
CD_in_CE_neg=[];
CD_in_CE_pos=[];
CD_out_CE_neg=[];
CD_out_CE_pos=[];
labels={'EC In','EC Out','DG In','DG Out','CA3 In','CA3 Out','CA1 In','CA1 Out'};

for i=1:4
    CD_in=[CD_in;[effect_in{i}{1,1}{1}.Effect,effect_in{i}{1,2}{1}.Effect,effect_in{i}{1,3}{1}.Effect]];
    CD_in_CE_neg=[CD_in_CE_neg;[effect_in{i}{1,1}{1}.ConfidenceIntervals(1),effect_in{i}{1,2}{1}.ConfidenceIntervals(1),effect_in{i}{1,3}{1}.ConfidenceIntervals(1)]];
    CD_in_CE_pos=[CD_in_CE_pos;[effect_in{i}{1,1}{1}.ConfidenceIntervals(2),effect_in{i}{1,2}{1}.ConfidenceIntervals(2),effect_in{i}{1,3}{1}.ConfidenceIntervals(2)]];
    CD_out=[CD_out;[effect_out{i}{1,1}{1}.Effect,effect_out{i}{1,2}{1}.Effect,effect_out{i}{1,3}{1}.Effect]];
    CD_out_CE_neg=[CD_out_CE_neg;[effect_out{i}{1,1}{1}.ConfidenceIntervals(1),effect_out{i}{1,2}{1}.ConfidenceIntervals(1),effect_out{i}{1,3}{1}.ConfidenceIntervals(1)]];
    CD_out_CE_pos=[CD_out_CE_pos;[effect_out{i}{1,1}{1}.ConfidenceIntervals(2),effect_out{i}{1,2}{1}.ConfidenceIntervals(2),effect_out{i}{1,3}{1}.ConfidenceIntervals(2)]];
end

CD_m=interleave_rows(CD_in,CD_out);
CD_CE_neg=interleave_rows(CD_in_CE_neg,CD_out_CE_neg);
CD_CE_pos=interleave_rows(CD_in_CE_pos,CD_out_CE_pos);

figure( 'Position', [100 100 700 600])
plot_horzbar_CohenD...
    (CD_m,CD_CE_neg,CD_CE_pos,labels);
xlabel("Slope Centrality Directional Mean Effect Size")
xlim([-3.2,3.2])
ylim([0,9])
set(gca,'LineWidth',2,'TickLength',[0.02 0.02]);
set(gca,'FontSize',18)
exportgraphics(gcf,".\Temporal Analysis\Graph Analysis\slope_centrality_cohenD.png",'Resolution',1500)
%% tunnel in/out

% graphing degrees averaged
node_degrees_ff_well=[];
node_degrees_ff_well_se=[];
node_degrees_fb_well=[];
node_degrees_fb_well_se=[];
node_degrees_ff_tunnel=[];
node_degrees_ff_tunnel_se=[];
node_degrees_fb_tunnel=[];
node_degrees_fb_tunnel_se=[];

ff_in_well_m=[];
ff_in_well_se=[];
fb_in_well_m=[];
fb_in_well_se=[];
ff_out_well_m=[];
ff_out_well_se=[];
fb_out_well_m=[];
fb_out_well_se=[];

for regi=1:4
    
    %ff out tunnel
    nostim_well_ff_degs=nonzeros(cell2mat(nostim_tunnel_in_ff{:,regi}));
    hfs5_well_ff_degs=nonzeros(cell2mat(hfs5_tunnel_in_ff{:,regi}));
    hfs40_well_ff_degs=nonzeros(cell2mat(hfs40_tunnel_in_ff{:,regi}));
    node_degrees_ff_well=[node_degrees_ff_well;...
        [mean(nostim_well_ff_degs),mean(hfs5_well_ff_degs),mean(hfs40_well_ff_degs)]];
    node_degrees_ff_well_se=[node_degrees_ff_well_se;...
        [std(nostim_well_ff_degs)/sqrt(length(nostim_well_ff_degs)),...
        std(hfs5_well_ff_degs)/sqrt(length(hfs5_well_ff_degs)),...
        std(hfs40_well_ff_degs)/sqrt(length(hfs40_well_ff_degs))]];
    
    %fb out tunnel
    nostim_well_fb_degs=nonzeros(cell2mat(nostim_tunnel_in_fb{:,regi}));
    hfs5_well_fb_degs=nonzeros(cell2mat(hfs5_tunnel_in_fb{:,regi}));
    hfs40_well_fb_degs=nonzeros(cell2mat(hfs40_tunnel_in_fb{:,regi}));
    node_degrees_fb_well=[node_degrees_fb_well;...
        [mean(nostim_well_fb_degs),mean(hfs5_well_fb_degs),mean(hfs40_well_fb_degs)]];
    node_degrees_fb_well_se=[node_degrees_fb_well_se;...
        [std(nostim_well_fb_degs)/sqrt(length(nostim_well_fb_degs)),...
        std(hfs5_well_fb_degs)/sqrt(length(hfs5_well_fb_degs)),...
        std(hfs40_well_fb_degs)/sqrt(length(hfs40_well_fb_degs))]];
    
    %ff in tunnel
    nostim_tunnel_ff_degs=nonzeros(cell2mat(nostim_tunnel_out_ff{:,regi}));
    hfs5_tunnel_ff_degs=nonzeros(cell2mat(hfs5_tunnel_out_ff{:,regi}));
    hfs40_tunnel_ff_degs=nonzeros(cell2mat(hfs40_tunnel_out_ff{:,regi}));
    node_degrees_ff_tunnel=[node_degrees_ff_tunnel;...
        [mean(nostim_tunnel_ff_degs),mean(hfs5_tunnel_ff_degs),mean(hfs40_tunnel_ff_degs)]];
    node_degrees_ff_tunnel_se=[node_degrees_ff_tunnel_se;...
        [std(nostim_tunnel_ff_degs)/sqrt(length(nostim_tunnel_ff_degs)),...
        std(hfs5_tunnel_ff_degs)/sqrt(length(hfs5_tunnel_ff_degs)),...
        std(hfs40_tunnel_ff_degs)/sqrt(length(hfs40_tunnel_ff_degs))]];
    
    %fb in tunnel
    nostim_tunnel_fb_degs=nonzeros(cell2mat(nostim_tunnel_out_fb{:,regi}));
    hfs5_tunnel_fb_degs=nonzeros(cell2mat(hfs5_tunnel_out_fb{:,regi}));
    hfs40_tunnel_fb_degs=nonzeros(cell2mat(hfs40_tunnel_out_fb{:,regi}));
    node_degrees_fb_tunnel=[node_degrees_fb_tunnel;...
        [mean(nostim_tunnel_fb_degs),mean(hfs5_tunnel_fb_degs),mean(hfs40_tunnel_fb_degs)]];
    node_degrees_fb_tunnel_se=[node_degrees_fb_tunnel_se;...
        [std(nostim_tunnel_fb_degs)/sqrt(length(nostim_tunnel_fb_degs)),...
        std(hfs5_tunnel_fb_degs)/sqrt(length(hfs5_tunnel_fb_degs)),...
        std(hfs40_tunnel_fb_degs)/sqrt(length(hfs40_tunnel_fb_degs))]];
    
    % ff in anova
    name_mat=[repmat("No Stim",length(nostim_well_ff_degs),1);repmat("5 HFS",length(hfs5_well_ff_degs),1);...
        repmat("40 HFS",length(hfs40_well_ff_degs),1)];
    [~,~,stats]=anova1([nostim_well_ff_degs;hfs5_well_ff_degs;hfs40_well_ff_degs],name_mat);
    [c,m]=multcompare(stats);
    ff_in_well_c{regi}=c;
    ff_in_well_m=[ff_in_well_m;m(:,1)'];
    ff_in_well_se=[ff_in_well_se;m(:,2)'];
    
    % fb in anova
    name_mat=[repmat("No Stim",length(nostim_well_fb_degs),1);repmat("5 HFS",length(hfs5_well_fb_degs),1);...
        repmat("40 HFS",length(hfs40_well_fb_degs),1)];
    [~,~,stats]=anova1([nostim_well_fb_degs;hfs5_well_fb_degs;hfs40_well_fb_degs],name_mat);
    [c,m]=multcompare(stats);
    fb_in_well_c{regi}=c;
    fb_in_well_m=[fb_in_well_m;m(:,1)'];
    fb_in_well_se=[fb_in_well_se;m(:,2)'];
    
    % ff out anova
    name_mat=[repmat("No Stim",length(nostim_tunnel_ff_degs),1);repmat("5 HFS",length(hfs5_tunnel_ff_degs),1);...
        repmat("40 HFS",length(hfs40_tunnel_ff_degs),1)];
    [~,~,stats]=anova1([nostim_tunnel_ff_degs;hfs5_tunnel_ff_degs;hfs40_tunnel_ff_degs],name_mat);
    [c,m]=multcompare(stats);
    ff_out_well_c{regi}=c;
    ff_out_well_m=[ff_out_well_m;m(:,1)'];
    ff_out_well_se=[ff_out_well_se;m(:,2)'];
    
    % fb in anova
    name_mat=[repmat("No Stim",length(nostim_tunnel_fb_degs),1);repmat("5 HFS",length(hfs5_tunnel_fb_degs),1);...
        repmat("40 HFS",length(hfs40_tunnel_fb_degs),1)];
    [~,~,stats]=anova1([nostim_tunnel_fb_degs;hfs5_tunnel_fb_degs;hfs40_tunnel_fb_degs],name_mat);
    [c,m]=multcompare(stats);
    fb_out_well_c{regi}=c;
    fb_out_well_m=[fb_out_well_m;m(:,1)'];
    fb_out_well_se=[fb_out_well_se;m(:,2)'];
    
    %in ff fb anova
    v=[nostim_well_ff_degs;nostim_well_fb_degs];
    g=[repmat("ff",length(nostim_well_ff_degs),1);repmat("fb",length(nostim_well_fb_degs),1)];
    [~,~,stats]=anova1(v,g);
    [c,~]=multcompare(stats);
    ff_fb_in_c{regi,1}=c;
    
    v=[hfs5_well_ff_degs;hfs5_well_fb_degs];
    g=[repmat("ff",length(hfs5_well_ff_degs),1);repmat("fb",length(hfs5_well_fb_degs),1)];
    [~,~,stats]=anova1(v,g);
    [c,~]=multcompare(stats);
    ff_fb_in_c{regi,2}=c;
    
    v=[hfs40_well_ff_degs;hfs40_well_fb_degs];
    g=[repmat("ff",length(hfs40_well_ff_degs),1);repmat("fb",length(hfs40_well_fb_degs),1)];
    [~,~,stats]=anova1(v,g);
    [c,~]=multcompare(stats);
    ff_fb_in_c{regi,3}=c;
    
    %out ff fb anova
    v=[nostim_tunnel_ff_degs;nostim_tunnel_fb_degs];
    g=[repmat("ff",length(nostim_tunnel_ff_degs),1);repmat("fb",length(nostim_tunnel_fb_degs),1)];
    [~,~,stats]=anova1(v,g);
    [c,~]=multcompare(stats);
    ff_fb_out_c{regi,1}=c;
    
    v=[hfs5_tunnel_ff_degs;hfs5_tunnel_fb_degs];
    g=[repmat("ff",length(hfs5_tunnel_ff_degs),1);repmat("fb",length(hfs5_tunnel_fb_degs),1)];
    [~,~,stats]=anova1(v,g);
    [c,~]=multcompare(stats);
    ff_fb_out_c{regi,2}=c;
    
    v=[hfs40_tunnel_ff_degs;hfs40_tunnel_fb_degs];
    g=[repmat("ff",length(hfs40_tunnel_ff_degs),1);repmat("fb",length(hfs40_tunnel_fb_degs),1)];
    [~,~,stats]=anova1(v,g);
    [c,~]=multcompare(stats);
    ff_fb_out_c{regi,3}=c;

    effect_in{regi}=ff_fb_meanEffect(nostim_well_ff_degs,nostim_well_fb_degs,hfs5_well_ff_degs,hfs5_well_fb_degs,hfs40_well_ff_degs,hfs40_well_fb_degs);
    effect_out{regi}=ff_fb_meanEffect(nostim_tunnel_ff_degs,nostim_tunnel_fb_degs,hfs5_tunnel_ff_degs,hfs5_tunnel_fb_degs,hfs40_tunnel_ff_degs,hfs40_tunnel_fb_degs);
end

ff_fb_c = reshape([ff_fb_in_c(:) ff_fb_out_c(:)]',2*size(ff_fb_in_c,1), []);
get_last=@(x) x(1,end);
ff_fb_p = cellfun(get_last,ff_fb_c);

% degs={node_degrees_ff_well,node_degrees_fb_well,node_degrees_ff_tunnel,node_degrees_fb_tunnel};
% degs_se={node_degrees_ff_well_se,node_degrees_fb_well_se,node_degrees_ff_tunnel_se,node_degrees_fb_tunnel_se};

degs={ff_in_well_m,fb_in_well_m,ff_out_well_m,fb_out_well_m};
degs_se={ff_in_well_se,fb_in_well_se,ff_out_well_se,fb_out_well_se};

figure('Position',[0,0,800,800])
t=tiledlayout(2,2,'TileSpacing','compact');
ylabel(t,"Centrality Per Node",'FontSize',24)
xlabel(t,"Subregion",'FontSize',24)
title(t,"Tunnel Slope Centrality",'FontSize',24)
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
        title("Feedforward Out Centrality",'FontSize',20)
        b=bar(regs_tunnel_ff_cat,degs{i});
    elseif i==4
        title("Feedback Out Centrality",'FontSize',20)
        b=bar(regs_tunnel_fb_cat,degs{i});
    elseif i==1
        title("Feedforward In Centrality",'FontSize',20)
        b=bar(regs_tunnel_ff_cat,degs{i});
    elseif i==2
        title("Feedback In Centrality",'FontSize',20)
        b=bar(regs_tunnel_fb_cat,degs{i});
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
% saveas(gcf,"D:\Brewer lab data\HFS\Temporal Analysis\Graph Analysis\centrality tunnel nonzero.png")
%% horzbar from anova
dat2plot_ff=[];
dat2plot_ff_se=[];
dat2plot_fb=[];
dat2plot_fb_se=[];
labels={'EC-DG In','EC-DG Out','DG-CA3 In','DG-CA3 Out','CA3-CA1 In','CA3-CA1 Out','CA1-EC In','CA1-EC Out'};
figure( 'Position', [100 100 1400 600])
for i=1:4
    dat2plot_ff=[dat2plot_ff;degs{1}(i,:);degs{3}(i,:)];
    dat2plot_ff_se=[dat2plot_ff_se;degs_se{1}(i,:);degs_se{3}(i,:)];
    dat2plot_fb=[dat2plot_fb;degs{2}(i,:);degs{4}(i,:)];
    dat2plot_fb_se=[dat2plot_fb_se;degs_se{2}(i,:);degs_se{4}(i,:)];
end

plot_horzbar_from_anova...
    (dat2plot_ff,dat2plot_ff_se,dat2plot_fb,dat2plot_fb_se,labels,ff_colors,fb_colors)
xlabel("FB Edges                  FF Edges")
xlim([-18,18])
xticks([-18:3:18])
ylim([0,9])
set(gca,'LineWidth',2,'TickLength',[0.02 0.02]);
set(gca,'FontSize',18)
% exportgraphics(gcf,".\three stim share figs\slope_tunnel_centrality_comp.png",'Resolution',1500)
exportgraphics(gcf,".\Temporal Analysis\Graph Analysis\slope_tunnel_centrality_comp.png",'Resolution',1500)

mean_diff=dat2plot_ff-dat2plot_fb;
%% get centrality nonzero rsq
%out edges
[nostim_well_out_ff,nostim_well_out_fb,nostim_tunnel_out_ff,nostim_tunnel_out_fb]=...
    get_centrality_out(nostim_rsq_G,nostim_datatable,nostim_data_fldr,nostim_temporal_fldr);

[hfs5_well_out_ff,hfs5_well_out_fb,hfs5_tunnel_out_ff,hfs5_tunnel_out_fb]=...
    get_centrality_out(hfs5_rsq_G,hfs5_datatable,hfs5_data_fldr,hfs5_temporal_fldr);

[hfs40_well_out_ff,hfs40_well_out_fb,hfs40_tunnel_out_ff,hfs40_tunnel_out_fb]=...
    get_centrality_out(hfs40_rsq_G,hfs40_datatable,hfs40_data_fldr,hfs40_temporal_fldr);
% in edges
[nostim_well_in_ff,nostim_well_in_fb,nostim_tunnel_in_ff,nostim_tunnel_in_fb]=...
    get_centrality_in(nostim_rsq_G,nostim_datatable,nostim_data_fldr,nostim_temporal_fldr);

[hfs5_well_in_ff,hfs5_well_in_fb,hfs5_tunnel_in_ff,hfs5_tunnel_in_fb]=...
    get_centrality_in(hfs5_rsq_G,hfs5_datatable,hfs5_data_fldr,hfs5_temporal_fldr);

[hfs40_well_in_ff,hfs40_well_in_fb,hfs40_tunnel_in_ff,hfs40_tunnel_in_fb]=...
    get_centrality_in(hfs40_rsq_G,hfs40_datatable,hfs40_data_fldr,hfs40_temporal_fldr);

% graphing degrees averaged
node_degrees_ff_well=[];
node_degrees_ff_well_se=[];
node_degrees_fb_well=[];
node_degrees_fb_well_se=[];
node_degrees_ff_tunnel=[];
node_degrees_ff_tunnel_se=[];
node_degrees_fb_tunnel=[];
node_degrees_fb_tunnel_se=[];

%% well in/out
ff_in_well_m=[];
ff_in_well_se=[];
fb_in_well_m=[];
fb_in_well_se=[];
ff_out_well_m=[];
ff_out_well_se=[];
fb_out_well_m=[];
fb_out_well_se=[];

for regi=1:4
    
    %out
    %ff wells out
    nostim_well_ff_degs=nonzeros(cell2mat(nostim_well_in_ff{:,regi}));
    hfs5_well_ff_degs=nonzeros(cell2mat(hfs5_well_in_ff{:,regi}));
    hfs40_well_ff_degs=nonzeros(cell2mat(hfs40_well_in_ff{:,regi}));
    node_degrees_ff_well=[node_degrees_ff_well;...
        [mean(nostim_well_ff_degs),mean(hfs5_well_ff_degs),mean(hfs40_well_ff_degs)]];
    node_degrees_ff_well_se=[node_degrees_ff_well_se;...
        [std(nostim_well_ff_degs)/sqrt(length(nostim_well_ff_degs)),...
        std(hfs5_well_ff_degs)/sqrt(length(hfs5_well_ff_degs)),...
        std(hfs40_well_ff_degs)/sqrt(length(hfs40_well_ff_degs))]];
    
    %fb wells out
    nostim_well_fb_degs=nonzeros(cell2mat(nostim_well_in_fb{:,regi}));
    hfs5_well_fb_degs=nonzeros(cell2mat(hfs5_well_in_fb{:,regi}));
    hfs40_well_fb_degs=nonzeros(cell2mat(hfs40_well_in_fb{:,regi}));
    node_degrees_fb_well=[node_degrees_fb_well;...
        [mean(nostim_well_fb_degs),mean(hfs5_well_fb_degs),mean(hfs40_well_fb_degs)]];
    node_degrees_fb_well_se=[node_degrees_fb_well_se;...
        [std(nostim_well_fb_degs)/sqrt(length(nostim_well_fb_degs)),...
        std(hfs5_well_fb_degs)/sqrt(length(hfs5_well_fb_degs)),...
        std(hfs40_well_fb_degs)/sqrt(length(hfs40_well_fb_degs))]];
    
    %in
    %ff in
    nostim_tunnel_ff_degs=nonzeros(cell2mat(nostim_well_out_ff{:,regi}));
    hfs5_tunnel_ff_degs=nonzeros(cell2mat(hfs5_well_out_ff{:,regi}));
    hfs40_tunnel_ff_degs=nonzeros(cell2mat(hfs40_well_out_ff{:,regi}));
    node_degrees_ff_tunnel=[node_degrees_ff_tunnel;...
        [mean(nostim_tunnel_ff_degs),mean(hfs5_tunnel_ff_degs),mean(hfs40_tunnel_ff_degs)]];
    node_degrees_ff_tunnel_se=[node_degrees_ff_tunnel_se;...
        [std(nostim_tunnel_ff_degs)/sqrt(length(nostim_tunnel_ff_degs)),...
        std(hfs5_tunnel_ff_degs)/sqrt(length(hfs5_tunnel_ff_degs)),...
        std(hfs40_tunnel_ff_degs)/sqrt(length(hfs40_tunnel_ff_degs))]];
    
    %fb in
    nostim_tunnel_fb_degs=nonzeros(cell2mat(nostim_well_out_fb{:,regi}));
    hfs5_tunnel_fb_degs=nonzeros(cell2mat(hfs5_well_out_fb{:,regi}));
    hfs40_tunnel_fb_degs=nonzeros(cell2mat(hfs40_well_out_fb{:,regi}));
    node_degrees_fb_tunnel=[node_degrees_fb_tunnel;...
        [mean(nostim_tunnel_fb_degs),mean(hfs5_tunnel_fb_degs),mean(hfs40_tunnel_fb_degs)]];
    node_degrees_fb_tunnel_se=[node_degrees_fb_tunnel_se;...
        [std(nostim_tunnel_fb_degs)/sqrt(length(nostim_tunnel_fb_degs)),...
        std(hfs5_tunnel_fb_degs)/sqrt(length(hfs5_tunnel_fb_degs)),...
        std(hfs40_tunnel_fb_degs)/sqrt(length(hfs40_tunnel_fb_degs))]];
    
    % ff in anova
    name_mat=[repmat("No Stim",length(nostim_well_ff_degs),1);repmat("5 HFS",length(hfs5_well_ff_degs),1);...
        repmat("40 HFS",length(hfs40_well_ff_degs),1)];
    [~,~,stats]=anova1([nostim_well_ff_degs;hfs5_well_ff_degs;hfs40_well_ff_degs],name_mat);
    [c,m]=multcompare(stats);
    ff_in_well_c{regi}=c;
    ff_in_well_m=[ff_in_well_m;m(:,1)'];
    ff_in_well_se=[ff_in_well_se;m(:,2)'];
    
    % fb in anova
    name_mat=[repmat("No Stim",length(nostim_well_fb_degs),1);repmat("5 HFS",length(hfs5_well_fb_degs),1);...
        repmat("40 HFS",length(hfs40_well_fb_degs),1)];
    [~,~,stats]=anova1([nostim_well_fb_degs;hfs5_well_fb_degs;hfs40_well_fb_degs],name_mat);
    [c,m]=multcompare(stats);
    fb_in_well_c{regi}=c;
    fb_in_well_m=[fb_in_well_m;m(:,1)'];
    fb_in_well_se=[fb_in_well_se;m(:,2)'];
    
    % ff out anova
    name_mat=[repmat("No Stim",length(nostim_tunnel_ff_degs),1);repmat("5 HFS",length(hfs5_tunnel_ff_degs),1);...
        repmat("40 HFS",length(hfs40_tunnel_ff_degs),1)];
    [~,~,stats]=anova1([nostim_tunnel_ff_degs;hfs5_tunnel_ff_degs;hfs40_tunnel_ff_degs],name_mat);
    [c,m]=multcompare(stats);
    ff_out_well_c{regi}=c;
    ff_out_well_m=[ff_out_well_m;m(:,1)'];
    ff_out_well_se=[ff_out_well_se;m(:,2)'];
    
    % fb in anova
    name_mat=[repmat("No Stim",length(nostim_tunnel_fb_degs),1);repmat("5 HFS",length(hfs5_tunnel_fb_degs),1);...
        repmat("40 HFS",length(hfs40_tunnel_fb_degs),1)];
    [~,~,stats]=anova1([nostim_tunnel_fb_degs;hfs5_tunnel_fb_degs;hfs40_tunnel_fb_degs],name_mat);
    [c,m]=multcompare(stats);
    fb_out_well_c{regi}=c;
    fb_out_well_m=[fb_out_well_m;m(:,1)'];
    fb_out_well_se=[fb_out_well_se;m(:,2)'];
    
    %in ff fb anova
    v=[nostim_well_ff_degs;nostim_well_fb_degs];
    g=[repmat("ff",length(nostim_well_ff_degs),1);repmat("fb",length(nostim_well_fb_degs),1)];
    [~,~,stats]=anova1(v,g);
    [c,~]=multcompare(stats);
    ff_fb_in_c{regi,1}=c;
    
    v=[hfs5_well_ff_degs;hfs5_well_fb_degs];
    g=[repmat("ff",length(hfs5_well_ff_degs),1);repmat("fb",length(hfs5_well_fb_degs),1)];
    [~,~,stats]=anova1(v,g);
    [c,~]=multcompare(stats);
    ff_fb_in_c{regi,2}=c;
    
    v=[hfs40_well_ff_degs;hfs40_well_fb_degs];
    g=[repmat("ff",length(hfs40_well_ff_degs),1);repmat("fb",length(hfs40_well_fb_degs),1)];
    [~,~,stats]=anova1(v,g);
    [c,~]=multcompare(stats);
    ff_fb_in_c{regi,3}=c;
    
    %out ff fb anova
    v=[nostim_tunnel_ff_degs;nostim_tunnel_fb_degs];
    g=[repmat("ff",length(nostim_tunnel_ff_degs),1);repmat("fb",length(nostim_tunnel_fb_degs),1)];
    [~,~,stats]=anova1(v,g);
    [c,~]=multcompare(stats);
    ff_fb_out_c{regi,1}=c;
    
    v=[hfs5_tunnel_ff_degs;hfs5_tunnel_fb_degs];
    g=[repmat("ff",length(hfs5_tunnel_ff_degs),1);repmat("fb",length(hfs5_tunnel_fb_degs),1)];
    [~,~,stats]=anova1(v,g);
    [c,~]=multcompare(stats);
    ff_fb_out_c{regi,2}=c;
    
    v=[hfs40_tunnel_ff_degs;hfs40_tunnel_fb_degs];
    g=[repmat("ff",length(hfs40_tunnel_ff_degs),1);repmat("fb",length(hfs40_tunnel_fb_degs),1)];
    [~,~,stats]=anova1(v,g);
    [c,~]=multcompare(stats);
    ff_fb_out_c{regi,3}=c;

    effect_in{regi}=ff_fb_meanEffect(nostim_well_ff_degs,nostim_well_fb_degs,hfs5_well_ff_degs,hfs5_well_fb_degs,hfs40_well_ff_degs,hfs40_well_fb_degs);
    effect_out{regi}=ff_fb_meanEffect(nostim_tunnel_ff_degs,nostim_tunnel_fb_degs,hfs5_tunnel_ff_degs,hfs5_tunnel_fb_degs,hfs40_tunnel_ff_degs,hfs40_tunnel_fb_degs);
end

ff_fb_c = reshape([ff_fb_in_c(:) ff_fb_out_c(:)]',2*size(ff_fb_in_c,1), []);
get_last=@(x) x(1,end);
ff_fb_p = cellfun(get_last,ff_fb_c);

% degs={node_degrees_ff_well,node_degrees_fb_well,node_degrees_ff_tunnel,node_degrees_fb_tunnel};
% degs_se={node_degrees_ff_well_se,node_degrees_fb_well_se,node_degrees_ff_tunnel_se,node_degrees_fb_tunnel_se};

degs={ff_in_well_m,fb_in_well_m,ff_out_well_m,fb_out_well_m};
degs_se={ff_in_well_se,fb_in_well_se,ff_out_well_se,fb_out_well_se};

figure('Position',[0,0,800,800])
t=tiledlayout(2,2,'TileSpacing','compact');
ylabel(t,"Centrality Per Node",'FontSize',24)
xlabel(t,"Subregion",'FontSize',24)
title(t,"Well R^2 Centrality",'FontSize',24)
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
        title("Feedforward Out Centrality",'FontSize',20)
        b=bar(regs_cat,degs{i});
    elseif i==4
        title("Feedback Out Centrality",'FontSize',20)
        b=bar(regs_cat,degs{i});
    elseif i==1
        title("Feedforward In Centrality",'FontSize',20)
        b=bar(regs_cat,degs{i});
    elseif i==2
        title("Feedback In Centrality",'FontSize',20)
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
% saveas(gcf,"D:\Brewer lab data\HFS\Temporal Analysis\Graph Analysis\centrality well nonzero rsq.png")
%% horzbar from anova
dat2plot_ff=[];
dat2plot_ff_se=[];
dat2plot_fb=[];
dat2plot_fb_se=[];
labels={'EC In','EC Out','DG In','DG Out','CA3 In','CA3 Out','CA1 In','CA1 Out'};
figure( 'Position', [100 100 1400 600])
for i=1:4
    dat2plot_ff=[dat2plot_ff;degs{1}(i,:);degs{3}(i,:)];
    dat2plot_ff_se=[dat2plot_ff_se;degs_se{1}(i,:);degs_se{3}(i,:)];
    dat2plot_fb=[dat2plot_fb;degs{2}(i,:);degs{4}(i,:)];
    dat2plot_fb_se=[dat2plot_fb_se;degs_se{2}(i,:);degs_se{4}(i,:)];
end

plot_horzbar_from_anova...
    (dat2plot_ff,dat2plot_ff_se,dat2plot_fb,dat2plot_fb_se,labels,ff_colors,fb_colors)
xlabel("FB R^2 Centrality                  FF R^2 Centrality")
xlim([-1.5,1.5])
% xticks([-50:10:50])
ylim([0,9])
set(gca,'LineWidth',2,'TickLength',[0.02 0.02]);
set(gca,'FontSize',18)
% exportgraphics(gcf,".\three stim share figs\rsq_well_centrality_comp.png",'Resolution',1500)
exportgraphics(gcf,".\Temporal Analysis\Graph Analysis\rsq_well_centrality_comp.png",'Resolution',1500)

mean_diff=dat2plot_ff-dat2plot_fb;

%% FF/FB% degrees

ff_fb_perc=(dat2plot_ff./dat2plot_fb)*100;

% ff_fb_SD=dat2plot

ff_fb_SD=[];
ff_fb_lengths=[];
for regi=1:4
    
    %out
    %ff wells out
    nostim_well_ff_degs=nonzeros(cell2mat(nostim_well_in_ff{:,regi}));
    hfs5_well_ff_degs=nonzeros(cell2mat(hfs5_well_in_ff{:,regi}));
    hfs40_well_ff_degs=nonzeros(cell2mat(hfs40_well_in_ff{:,regi}));
    node_degrees_ff_well_m=[mean(nostim_well_ff_degs),mean(hfs5_well_ff_degs),mean(hfs40_well_ff_degs)];
    node_degrees_ff_well_sd=[std(nostim_well_ff_degs),std(hfs5_well_ff_degs),std(hfs40_well_ff_degs)];
    ff_out_lengths=[length(nostim_well_ff_degs),length(hfs5_well_ff_degs),length(hfs40_well_ff_degs)];

    %fb wells out
    nostim_well_fb_degs=nonzeros(cell2mat(nostim_well_in_fb{:,regi}));
    hfs5_well_fb_degs=nonzeros(cell2mat(hfs5_well_in_fb{:,regi}));
    hfs40_well_fb_degs=nonzeros(cell2mat(hfs40_well_in_fb{:,regi}));
    node_degrees_fb_well_m=[mean(nostim_well_fb_degs),mean(hfs5_well_fb_degs),mean(hfs40_well_fb_degs)];
    node_degrees_fb_well_sd=[std(nostim_well_fb_degs),std(hfs5_well_fb_degs),std(hfs40_well_fb_degs)];
    fb_out_lengths=[length(nostim_well_fb_degs),length(hfs5_well_fb_degs),length(hfs40_well_fb_degs)];
    
    %in
    %ff in
    nostim_tunnel_ff_degs=nonzeros(cell2mat(nostim_well_out_ff{:,regi}));
    hfs5_tunnel_ff_degs=nonzeros(cell2mat(hfs5_well_out_ff{:,regi}));
    hfs40_tunnel_ff_degs=nonzeros(cell2mat(hfs40_well_out_ff{:,regi}));
    node_degrees_ff_tunnel_m=[mean(nostim_tunnel_ff_degs),mean(hfs5_tunnel_ff_degs),mean(hfs40_tunnel_ff_degs)];
    node_degrees_ff_tunnel_sd=[std(nostim_tunnel_ff_degs),std(hfs5_tunnel_ff_degs), std(hfs40_tunnel_ff_degs)];
    ff_in_lengths=[length(nostim_tunnel_ff_degs),length(hfs5_tunnel_ff_degs),length(hfs40_tunnel_ff_degs)];

    %fb in
    nostim_tunnel_fb_degs=nonzeros(cell2mat(nostim_well_out_fb{:,regi}));
    hfs5_tunnel_fb_degs=nonzeros(cell2mat(hfs5_well_out_fb{:,regi}));
    hfs40_tunnel_fb_degs=nonzeros(cell2mat(hfs40_well_out_fb{:,regi}));
    node_degrees_fb_tunnel_m=[mean(nostim_tunnel_fb_degs),mean(hfs5_tunnel_fb_degs),mean(hfs40_tunnel_fb_degs)];
    node_degrees_fb_tunnel_sd=[std(nostim_tunnel_fb_degs),std(hfs5_tunnel_fb_degs),std(hfs40_tunnel_fb_degs)];
    fb_in_lengths=[length(nostim_tunnel_fb_degs),length(hfs5_tunnel_fb_degs),length(hfs40_tunnel_fb_degs)];

    ff_fb_SD=[ff_fb_SD;((node_degrees_ff_well_m./node_degrees_fb_well_m)*100).*...
        sqrt((node_degrees_ff_well_sd./node_degrees_fb_well_sd).^2+...
        (node_degrees_fb_well_sd./node_degrees_fb_well_sd).^2);...
        ((node_degrees_ff_tunnel_m./node_degrees_fb_tunnel_m)*100).*...
        sqrt((node_degrees_ff_tunnel_m./node_degrees_fb_tunnel_m).^2+...
        (node_degrees_fb_tunnel_sd./node_degrees_fb_tunnel_sd).^2)];
    ff_fb_lengths=[ff_fb_lengths;(ff_in_lengths+fb_in_lengths);(ff_out_lengths+fb_out_lengths)];
end

labels={'EC-DG In','EC-DG Out','DG-CA3 In','DG-CA3 Out','CA3-CA1 In','CA3-CA1 Out','CA1-EC In','CA1-EC Out'};
regs=categorical(labels);
regs=reordercats(regs,labels);

figure( 'Position', [100 100 700 600])
b=bar(regs,ff_fb_perc);
b(1).BaseValue=100;
b(1).BaseLine.LineStyle = "--";
b(1).BaseLine.LineWidth=2;

ff_fb_SE=ff_fb_SD./sqrt(ff_fb_lengths);
hold on
% grouped error bars
[ngroups,nbars]=size(ff_fb_perc);
x=[];
for j=1:nbars
    x(j,:)=b(j).XEndPoints;
end
% errorbar(x',ff_fb_perc,ff_fb_SE,'k','LineStyle','none')
hold off
ylabel("% FF/FB")
ylim([40,140])
set(gca,'LineWidth',2,'TickLength',[0.02 0.02]);
set(gca,'FontSize',18)
exportgraphics(gcf,".\Temporal Analysis\Graph Analysis\rsq_centrality_ff_fb_perc.png",'Resolution',1500)
%% Cohen's D Bar Graphs
CD_in=[];
CD_out=[];
CD_in_CE_neg=[];
CD_in_CE_pos=[];
CD_out_CE_neg=[];
CD_out_CE_pos=[];
labels={'EC In','EC Out','DG In','DG Out','CA3 In','CA3 Out','CA1 In','CA1 Out'};

for i=1:4
    CD_in=[CD_in;[effect_in{i}{1,1}{1}.Effect,effect_in{i}{1,2}{1}.Effect,effect_in{i}{1,3}{1}.Effect]];
    CD_in_CE_neg=[CD_in_CE_neg;[effect_in{i}{1,1}{1}.ConfidenceIntervals(1),effect_in{i}{1,2}{1}.ConfidenceIntervals(1),effect_in{i}{1,3}{1}.ConfidenceIntervals(1)]];
    CD_in_CE_pos=[CD_in_CE_pos;[effect_in{i}{1,1}{1}.ConfidenceIntervals(2),effect_in{i}{1,2}{1}.ConfidenceIntervals(2),effect_in{i}{1,3}{1}.ConfidenceIntervals(2)]];
    CD_out=[CD_out;[effect_out{i}{1,1}{1}.Effect,effect_out{i}{1,2}{1}.Effect,effect_out{i}{1,3}{1}.Effect]];
    CD_out_CE_neg=[CD_out_CE_neg;[effect_out{i}{1,1}{1}.ConfidenceIntervals(1),effect_out{i}{1,2}{1}.ConfidenceIntervals(1),effect_out{i}{1,3}{1}.ConfidenceIntervals(1)]];
    CD_out_CE_pos=[CD_out_CE_pos;[effect_out{i}{1,1}{1}.ConfidenceIntervals(2),effect_out{i}{1,2}{1}.ConfidenceIntervals(2),effect_out{i}{1,3}{1}.ConfidenceIntervals(2)]];
end

CD_m=interleave_rows(CD_in,CD_out);
CD_CE_neg=interleave_rows(CD_in_CE_neg,CD_out_CE_neg);
CD_CE_pos=interleave_rows(CD_in_CE_pos,CD_out_CE_pos);

figure( 'Position', [100 100 700 600])
plot_horzbar_CohenD...
    (CD_m,CD_CE_neg,CD_CE_pos,labels);
xlabel("R^2 Centrality Directional Mean Effect Size")
xlim([-3,3])
ylim([0,9])
set(gca,'LineWidth',2,'TickLength',[0.02 0.02]);
set(gca,'FontSize',18)
exportgraphics(gcf,".\Temporal Analysis\Graph Analysis\rsq_centrality_cohenD.png",'Resolution',1500)
%% tunnel in/out

% graphing degrees averaged
node_degrees_ff_well=[];
node_degrees_ff_well_se=[];
node_degrees_fb_well=[];
node_degrees_fb_well_se=[];
node_degrees_ff_tunnel=[];
node_degrees_ff_tunnel_se=[];
node_degrees_fb_tunnel=[];
node_degrees_fb_tunnel_se=[];

ff_in_well_m=[];
ff_in_well_se=[];
fb_in_well_m=[];
fb_in_well_se=[];
ff_out_well_m=[];
ff_out_well_se=[];
fb_out_well_m=[];
fb_out_well_se=[];

for regi=1:4
    
    %ff out tunnel
    nostim_well_ff_degs=nonzeros(cell2mat(nostim_tunnel_in_ff{:,regi}));
    hfs5_well_ff_degs=nonzeros(cell2mat(hfs5_tunnel_in_ff{:,regi}));
    hfs40_well_ff_degs=nonzeros(cell2mat(hfs40_tunnel_in_ff{:,regi}));
    node_degrees_ff_well=[node_degrees_ff_well;...
        [mean(nostim_well_ff_degs),mean(hfs5_well_ff_degs),mean(hfs40_well_ff_degs)]];
    node_degrees_ff_well_se=[node_degrees_ff_well_se;...
        [std(nostim_well_ff_degs)/sqrt(length(nostim_well_ff_degs)),...
        std(hfs5_well_ff_degs)/sqrt(length(hfs5_well_ff_degs)),...
        std(hfs40_well_ff_degs)/sqrt(length(hfs40_well_ff_degs))]];
    
    %fb out tunnel
    nostim_well_fb_degs=nonzeros(cell2mat(nostim_tunnel_in_fb{:,regi}));
    hfs5_well_fb_degs=nonzeros(cell2mat(hfs5_tunnel_in_fb{:,regi}));
    hfs40_well_fb_degs=nonzeros(cell2mat(hfs40_tunnel_in_fb{:,regi}));
    node_degrees_fb_well=[node_degrees_fb_well;...
        [mean(nostim_well_fb_degs),mean(hfs5_well_fb_degs),mean(hfs40_well_fb_degs)]];
    node_degrees_fb_well_se=[node_degrees_fb_well_se;...
        [std(nostim_well_fb_degs)/sqrt(length(nostim_well_fb_degs)),...
        std(hfs5_well_fb_degs)/sqrt(length(hfs5_well_fb_degs)),...
        std(hfs40_well_fb_degs)/sqrt(length(hfs40_well_fb_degs))]];
    
    %ff in tunnel
    nostim_tunnel_ff_degs=nonzeros(cell2mat(nostim_tunnel_out_ff{:,regi}));
    hfs5_tunnel_ff_degs=nonzeros(cell2mat(hfs5_tunnel_out_ff{:,regi}));
    hfs40_tunnel_ff_degs=nonzeros(cell2mat(hfs40_tunnel_out_ff{:,regi}));
    node_degrees_ff_tunnel=[node_degrees_ff_tunnel;...
        [mean(nostim_tunnel_ff_degs),mean(hfs5_tunnel_ff_degs),mean(hfs40_tunnel_ff_degs)]];
    node_degrees_ff_tunnel_se=[node_degrees_ff_tunnel_se;...
        [std(nostim_tunnel_ff_degs)/sqrt(length(nostim_tunnel_ff_degs)),...
        std(hfs5_tunnel_ff_degs)/sqrt(length(hfs5_tunnel_ff_degs)),...
        std(hfs40_tunnel_ff_degs)/sqrt(length(hfs40_tunnel_ff_degs))]];
    
    %fb in tunnel
    nostim_tunnel_fb_degs=nonzeros(cell2mat(nostim_tunnel_out_fb{:,regi}));
    hfs5_tunnel_fb_degs=nonzeros(cell2mat(hfs5_tunnel_out_fb{:,regi}));
    hfs40_tunnel_fb_degs=nonzeros(cell2mat(hfs40_tunnel_out_fb{:,regi}));
    node_degrees_fb_tunnel=[node_degrees_fb_tunnel;...
        [mean(nostim_tunnel_fb_degs),mean(hfs5_tunnel_fb_degs),mean(hfs40_tunnel_fb_degs)]];
    node_degrees_fb_tunnel_se=[node_degrees_fb_tunnel_se;...
        [std(nostim_tunnel_fb_degs)/sqrt(length(nostim_tunnel_fb_degs)),...
        std(hfs5_tunnel_fb_degs)/sqrt(length(hfs5_tunnel_fb_degs)),...
        std(hfs40_tunnel_fb_degs)/sqrt(length(hfs40_tunnel_fb_degs))]];
    
    % ff in anova
    name_mat=[repmat("No Stim",length(nostim_well_ff_degs),1);repmat("5 HFS",length(hfs5_well_ff_degs),1);...
        repmat("40 HFS",length(hfs40_well_ff_degs),1)];
    [~,~,stats]=anova1([nostim_well_ff_degs;hfs5_well_ff_degs;hfs40_well_ff_degs],name_mat);
    [c,m]=multcompare(stats);
    ff_in_well_c{regi}=c;
    ff_in_well_m=[ff_in_well_m;m(:,1)'];
    ff_in_well_se=[ff_in_well_se;m(:,2)'];
    
    % fb in anova
    name_mat=[repmat("No Stim",length(nostim_well_fb_degs),1);repmat("5 HFS",length(hfs5_well_fb_degs),1);...
        repmat("40 HFS",length(hfs40_well_fb_degs),1)];
    [~,~,stats]=anova1([nostim_well_fb_degs;hfs5_well_fb_degs;hfs40_well_fb_degs],name_mat);
    [c,m]=multcompare(stats);
    fb_in_well_c{regi}=c;
    fb_in_well_m=[fb_in_well_m;m(:,1)'];
    fb_in_well_se=[fb_in_well_se;m(:,2)'];
    
    % ff out anova
    name_mat=[repmat("No Stim",length(nostim_tunnel_ff_degs),1);repmat("5 HFS",length(hfs5_tunnel_ff_degs),1);...
        repmat("40 HFS",length(hfs40_tunnel_ff_degs),1)];
    [~,~,stats]=anova1([nostim_tunnel_ff_degs;hfs5_tunnel_ff_degs;hfs40_tunnel_ff_degs],name_mat);
    [c,m]=multcompare(stats);
    ff_out_well_c{regi}=c;
    ff_out_well_m=[ff_out_well_m;m(:,1)'];
    ff_out_well_se=[ff_out_well_se;m(:,2)'];
    
    % fb in anova
    name_mat=[repmat("No Stim",length(nostim_tunnel_fb_degs),1);repmat("5 HFS",length(hfs5_tunnel_fb_degs),1);...
        repmat("40 HFS",length(hfs40_tunnel_fb_degs),1)];
    [~,~,stats]=anova1([nostim_tunnel_fb_degs;hfs5_tunnel_fb_degs;hfs40_tunnel_fb_degs],name_mat);
    [c,m]=multcompare(stats);
    fb_out_well_c{regi}=c;
    fb_out_well_m=[fb_out_well_m;m(:,1)'];
    fb_out_well_se=[fb_out_well_se;m(:,2)'];
    
    %in ff fb anova
    v=[nostim_well_ff_degs;nostim_well_fb_degs];
    g=[repmat("ff",length(nostim_well_ff_degs),1);repmat("fb",length(nostim_well_fb_degs),1)];
    [~,~,stats]=anova1(v,g);
    [c,~]=multcompare(stats);
    ff_fb_in_c{regi,1}=c;
    
    v=[hfs5_well_ff_degs;hfs5_well_fb_degs];
    g=[repmat("ff",length(hfs5_well_ff_degs),1);repmat("fb",length(hfs5_well_fb_degs),1)];
    [~,~,stats]=anova1(v,g);
    [c,~]=multcompare(stats);
    ff_fb_in_c{regi,2}=c;
    
    v=[hfs40_well_ff_degs;hfs40_well_fb_degs];
    g=[repmat("ff",length(hfs40_well_ff_degs),1);repmat("fb",length(hfs40_well_fb_degs),1)];
    [~,~,stats]=anova1(v,g);
    [c,~]=multcompare(stats);
    ff_fb_in_c{regi,3}=c;
    
    %out ff fb anova
    v=[nostim_tunnel_ff_degs;nostim_tunnel_fb_degs];
    g=[repmat("ff",length(nostim_tunnel_ff_degs),1);repmat("fb",length(nostim_tunnel_fb_degs),1)];
    [~,~,stats]=anova1(v,g);
    [c,~]=multcompare(stats);
    ff_fb_out_c{regi,1}=c;
    
    v=[hfs5_tunnel_ff_degs;hfs5_tunnel_fb_degs];
    g=[repmat("ff",length(hfs5_tunnel_ff_degs),1);repmat("fb",length(hfs5_tunnel_fb_degs),1)];
    [~,~,stats]=anova1(v,g);
    [c,~]=multcompare(stats);
    ff_fb_out_c{regi,2}=c;
    
    v=[hfs40_tunnel_ff_degs;hfs40_tunnel_fb_degs];
    g=[repmat("ff",length(hfs40_tunnel_ff_degs),1);repmat("fb",length(hfs40_tunnel_fb_degs),1)];
    [~,~,stats]=anova1(v,g);
    [c,~]=multcompare(stats);
    ff_fb_out_c{regi,3}=c;

    effect_in{regi}=ff_fb_meanEffect(nostim_well_ff_degs,nostim_well_fb_degs,hfs5_well_ff_degs,hfs5_well_fb_degs,hfs40_well_ff_degs,hfs40_well_fb_degs);
    effect_out{regi}=ff_fb_meanEffect(nostim_tunnel_ff_degs,nostim_tunnel_fb_degs,hfs5_tunnel_ff_degs,hfs5_tunnel_fb_degs,hfs40_tunnel_ff_degs,hfs40_tunnel_fb_degs);
end

ff_fb_c = reshape([ff_fb_in_c(:) ff_fb_out_c(:)]',2*size(ff_fb_in_c,1), []);
get_last=@(x) x(1,end);
ff_fb_p = cellfun(get_last,ff_fb_c);

% degs={node_degrees_ff_well,node_degrees_fb_well,node_degrees_ff_tunnel,node_degrees_fb_tunnel};
% degs_se={node_degrees_ff_well_se,node_degrees_fb_well_se,node_degrees_ff_tunnel_se,node_degrees_fb_tunnel_se};

degs={ff_in_well_m,fb_in_well_m,ff_out_well_m,fb_out_well_m};
degs_se={ff_in_well_se,fb_in_well_se,ff_out_well_se,fb_out_well_se};

figure('Position',[0,0,800,800])
t=tiledlayout(2,2,'TileSpacing','compact');
ylabel(t,"Centrality Per Node",'FontSize',24)
xlabel(t,"Subregion",'FontSize',24)
title(t,"Tunnel R^2 Centrality",'FontSize',24)
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
        title("Feedforward Out Centrality",'FontSize',20)
        b=bar(regs_tunnel_ff_cat,degs{i});
    elseif i==4
        title("Feedback Out Centrality",'FontSize',20)
        b=bar(regs_tunnel_fb_cat,degs{i});
    elseif i==1
        title("Feedforward In Centrality",'FontSize',20)
        b=bar(regs_tunnel_ff_cat,degs{i});
    elseif i==2
        title("Feedback In Centrality",'FontSize',20)
        b=bar(regs_tunnel_fb_cat,degs{i});
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
% saveas(gcf,"D:\Brewer lab data\HFS\Temporal Analysis\Graph Analysis\centrality tunnel nonzero rsq.png")

%% horzbar from anova
dat2plot_ff=[];
dat2plot_ff_se=[];
dat2plot_fb=[];
dat2plot_fb_se=[];
labels={'EC-DG In','EC-DG Out','DG-CA3 In','DG-CA3 Out','CA3-CA1 In','CA3-CA1 Out','CA1-EC In','CA1-EC Out'};
figure( 'Position', [100 100 1400 600])
for i=1:4
    dat2plot_ff=[dat2plot_ff;degs{1}(i,:);degs{3}(i,:)];
    dat2plot_ff_se=[dat2plot_ff_se;degs_se{1}(i,:);degs_se{3}(i,:)];
    dat2plot_fb=[dat2plot_fb;degs{2}(i,:);degs{4}(i,:)];
    dat2plot_fb_se=[dat2plot_fb_se;degs_se{2}(i,:);degs_se{4}(i,:)];
end

plot_horzbar_from_anova...
    (dat2plot_ff,dat2plot_ff_se,dat2plot_fb,dat2plot_fb_se,labels,ff_colors,fb_colors)
xlabel("FB Edges                  FF Edges")
xlim([-6,6])
% xticks([-50:10:50])
ylim([0,9])
set(gca,'LineWidth',2,'TickLength',[0.02 0.02]);
set(gca,'FontSize',18)
% exportgraphics(gcf,".\three stim share figs\rsq_tunnel_centrality_comp.png",'Resolution',1500)
exportgraphics(gcf,".\Temporal Analysis\Graph Analysis\rsq_tunnel_centrality_comp.png",'Resolution',1500)

mean_diff=dat2plot_ff-dat2plot_fb;
%% edges per array

[nostim_well_out_ff,nostim_well_out_fb,nostim_tunnel_out_ff,nostim_tunnel_out_fb]=...
    get_edgecount_out(nostim_slope_G,nostim_datatable,nostim_data_fldr,nostim_temporal_fldr);

[hfs5_well_out_ff,hfs5_well_out_fb,hfs5_tunnel_out_ff,hfs5_tunnel_out_fb]=...
    get_edgecount_out(hfs5_slope_G,hfs5_datatable,hfs5_data_fldr,hfs5_temporal_fldr);

[hfs40_well_out_ff,hfs40_well_out_fb,hfs40_tunnel_out_ff,hfs40_tunnel_out_fb]=...
    get_edgecount_out(hfs40_slope_G,hfs40_datatable,hfs40_data_fldr,hfs40_temporal_fldr);

sum_table=@(x)cellfun(@sum,x);

no_1=varfun(sum_table,nostim_well_out_ff(:,[1:4]));
no_2=varfun(sum_table,nostim_well_out_fb(:,[1:4]));
no_3=varfun(sum_table,nostim_tunnel_out_ff(:,[1:4]));
no_4=varfun(sum_table,nostim_tunnel_out_fb(:,[1:4]));

for i=1:4
    no_all(:,i)=sum([no_1{:,i},no_2{:,i},no_3{:,i},no_4{:,i}],2);
end

hfs5_1=varfun(sum_table,hfs5_well_out_ff(:,[1:4]));
hfs5_2=varfun(sum_table,hfs5_well_out_fb(:,[1:4]));
hfs5_3=varfun(sum_table,hfs5_tunnel_out_ff(:,[1:4]));
hfs5_4=varfun(sum_table,hfs5_tunnel_out_fb(:,[1:4]));

for i=1:4
    hfs5_all(:,i)=sum([hfs5_1{:,i},hfs5_2{:,i},hfs5_3{:,i},hfs5_4{:,i}],2);
end

hfs40_1=varfun(sum_table,hfs40_well_out_ff(:,[1:4]));
hfs40_2=varfun(sum_table,hfs40_well_out_fb(:,[1:4]));
hfs40_3=varfun(sum_table,hfs40_tunnel_out_ff(:,[1:4]));
hfs40_4=varfun(sum_table,hfs40_tunnel_out_fb(:,[1:4]));

for i=1:4
    hfs40_all(:,i)=sum([hfs40_1{:,i},hfs40_2{:,i},hfs40_3{:,i},hfs40_4{:,i}],2);
end

w_t_ff={no_1,hfs5_1,hfs40_1};
w_t_fb={no_2,hfs5_2,hfs40_2};
t_w_ff={no_3,hfs5_3,hfs40_3};
t_w_fb={no_4,hfs5_4,hfs40_4};

wells=["EC","DG","CA3","CA1"];
wells_cat=categorical(wells);
well_cat=reordercats(wells_cat,wells);

ff_tunnels=["EC-DG","DG-CA3","CA3-CA1","CA1-EC"];
ff_tunnels_cat=categorical(ff_tunnels);
ff_tunnels_cat=reordercats(ff_tunnels_cat,ff_tunnels);

fb_tunnels=["DG-EC","CA3-DG","CA1-CA3","EC-CA1"];
fb_tunnels_cat=categorical(fb_tunnels);
fb_tunnels_cat=reordercats(fb_tunnels_cat,fb_tunnels);

%all arrays
edges_per_subregion(w_t_ff,well_cat)
title("Well to Tunnel FF",'FontSize',20)
ylabel("Edges per Subregion",'FontSize',18)
% exportgraphics(gcf,'.\three stim share figs\w_t_ff_edges_per_array.png','Resolution',1500)
exportgraphics(gcf,'.\Temporal Analysis\Graph Analysis\w_t_ff_edges_per_array.png','Resolution',1500)

edges_per_subregion(w_t_fb,well_cat)
title("Well to Tunnel FB",'FontSize',20)
ylabel("Edges per Subregion",'FontSize',18)
% exportgraphics(gcf,'.\three stim share figs\w_t_fb_edges_per_array.png','Resolution',1500)
exportgraphics(gcf,'.\Temporal Analysis\Graph Analysis\w_t_fb_edges_per_array.png','Resolution',1500)

edges_per_subregion(t_w_ff,ff_tunnels_cat)
title("Tunnel to Well FF",'FontSize',20)
ylabel("Edges per Subregion",'FontSize',18)
% exportgraphics(gcf,'.\three stim share figs\t_w_ff_edges_per_array.png','Resolution',1500)
exportgraphics(gcf,'.\Temporal Analysis\Graph Analysis\t_w_ff_edges_per_array.png','Resolution',1500)

edges_per_subregion(t_w_fb,fb_tunnels_cat)
title("Tunnel to Well FB",'FontSize',20)
ylabel("Edges per Subregion",'FontSize',18)
% exportgraphics(gcf,'.\three stim share figs\t_w_fb_edges_per_array.png','Resolution',1500)
exportgraphics(gcf,'.\Temporal Analysis\Graph Analysis\t_w_fb_edges_per_array.png','Resolution',1500)

%% array 2
edges_per_subregion_array(w_t_ff,wells_cat,2)
title("Well to Tunnel FF",'FontSize',20)
ylabel("Edges per Subregion",'FontSize',18)
% exportgraphics(gcf,'.\three stim share figs\w_t_ff_edges_per_array_no2.png','Resolution',1500)
exportgraphics(gcf,'.\Temporal Analysis\Graph Analysis\w_t_ff_edges_per_array_no2.png','Resolution',1500)

edges_per_subregion_array(w_t_fb,wells_cat,2)
title("Well to Tunnel FB",'FontSize',20)
ylabel("Edges per Subregion",'FontSize',18)
% exportgraphics(gcf,'.\three stim share figs\w_t_fb_edges_per_array_no2.png','Resolution',1500)
exportgraphics(gcf,'.\Temporal Analysis\Graph Analysis\w_t_fb_edges_per_array_no2.png','Resolution',1500)

edges_per_subregion_array(t_w_ff,ff_tunnels_cat,2)
title("Tunnel to Well FF",'FontSize',20)
ylabel("Edges per Subregion",'FontSize',18)
% exportgraphics(gcf,'.\three stim share figs\t_w_ff_edges_per_array_no2.png','Resolution',1500)
exportgraphics(gcf,'.\Temporal Analysis\Graph Analysis\t_w_ff_edges_per_array_no2.png','Resolution',1500)

edges_per_subregion_array(t_w_fb,fb_tunnels_cat,2)
title("Tunnel to Well FB",'FontSize',20)
ylabel("Edges per Subregion",'FontSize',18)
% exportgraphics(gcf,'.\three stim share figs\t_w_fb_edges_per_array_no2.png','Resolution',1500)
exportgraphics(gcf,'.\Temporal Analysis\Graph Analysis\t_w_fb_edges_per_array_no2.png','Resolution',1500)

%% Horz Bar
labels={'EC/EC-DG','EC-DG/DG','DG/DG-CA3','DG-CA3/CA3','CA3/CA3-CA1','CA3-CA1/CA1','CA1/CA1-EC','CA1-EC/EC'};

%exclude data
datatables={nostim_datatable,hfs5_datatable,hfs40_datatable};
for i=1:3
    wt_ff_exclude{i}=exclude_FID(w_t_ff{i},datatables{i},10);
end
for i=1:3
    wt_fb_exclude{i}=exclude_FID(w_t_fb{i},datatables{i},10);
end
for i=1:3
    tw_ff_exclude{i}=exclude_FID(t_w_ff{i},datatables{i},10);
end
for i=1:3
    tw_fb_exclude{i}=exclude_FID(t_w_fb{i},datatables{i},10);
end

% plot_horzbar(w_t_ff,w_t_fb,t_w_ff,t_w_fb,labels)
plot_horzbar(wt_ff_exclude,wt_fb_exclude,tw_ff_exclude,tw_fb_exclude,labels)
xlabel("FB Edges                  FF Edges")
xlim([-50,50])
xticks([-50:10:50])
% exportgraphics(gcf,".\three stim share figs\edges_per_subregion_comp.png",'Resolution',1500)
exportgraphics(gcf,".\Temporal Analysis\Graph Analysis\edges_per_subregion_comp.png",'Resolution',1500)

%% subregion anova edges
wt_ff_m=[];
wt_ff_se=[];
for i=1:4
    vector=[wt_ff_exclude{1}{:,i};wt_ff_exclude{2}{:,i};wt_ff_exclude{3}{:,i}];
    g=[repmat("No Stim",length(wt_ff_exclude{1}{:,i}),1);repmat("5 HFS",length(wt_ff_exclude{2}{:,i}),1);repmat("40 HFS",length(wt_ff_exclude{3}{:,i}),1)];
    [~,~,stats]=anova1(vector,g);
    [c,m]=multcompare(stats);
    wt_ff_c{i}=c;
    wt_ff_m=[wt_ff_m;m(:,1)'];
    wt_ff_se=[wt_ff_se;m(:,2)'];
end

wt_fb_m=[];
wt_fb_se=[];
for i=1:4
    vector=[wt_fb_exclude{1}{:,i};wt_fb_exclude{2}{:,i};wt_fb_exclude{3}{:,i}];
    g=[repmat("No Stim",length(wt_fb_exclude{1}{:,i}),1);repmat("5 HFS",length(wt_fb_exclude{2}{:,i}),1);repmat("40 HFS",length(wt_fb_exclude{3}{:,i}),1)];
    [~,~,stats]=anova1(vector,g);
    [c,m]=multcompare(stats);
    wt_fb_c{i}=c;
    wt_fb_m=[wt_fb_m;m(:,1)'];
    wt_fb_se=[wt_fb_se;m(:,2)'];
end

tw_ff_m=[];
tw_ff_se=[];
for i=1:4
    vector=[tw_ff_exclude{1}{:,i};tw_ff_exclude{2}{:,i};tw_ff_exclude{3}{:,i}];
    g=[repmat("No Stim",length(tw_ff_exclude{1}{:,i}),1);repmat("5 HFS",length(tw_ff_exclude{2}{:,i}),1);repmat("40 HFS",length(tw_ff_exclude{3}{:,i}),1)];
    [~,~,stats]=anova1(vector,g);
    [c,m]=multcompare(stats);
    tw_ff_c{i}=c;
    tw_ff_m=[tw_ff_m;m(:,1)'];
    tw_ff_se=[tw_ff_se;m(:,2)'];
end

tw_fb_m=[];
tw_fb_se=[];
for i=1:4
    vector=[tw_fb_exclude{1}{:,i};tw_fb_exclude{2}{:,i};tw_fb_exclude{3}{:,i}];
    g=[repmat("No Stim",length(tw_fb_exclude{1}{:,i}),1);repmat("5 HFS",length(tw_fb_exclude{2}{:,i}),1);repmat("40 HFS",length(tw_fb_exclude{3}{:,i}),1)];
    [~,~,stats]=anova1(vector,g);
    [c,m]=multcompare(stats);
    tw_fb_c{i}=c;
    tw_fb_m=[tw_fb_m;m(:,1)'];
    tw_fb_se=[tw_fb_se;m(:,2)'];
end

%% FF FB ANOVA
wt_ff_fb_c=[];
tw_ff_fb_c=[];
wt_ff_fb_m=[];
tw_ff_fb_m=[];
wt_ff_fb_se=[];
tw_ff_fb_se=[];

for subregion=1:4
    for stims=1:3
        vec_wt=[wt_ff_exclude{stims}{:,subregion},wt_fb_exclude{stims}{:,subregion}];
        vec_tw=[tw_ff_exclude{stims}{:,subregion},tw_fb_exclude{stims}{:,subregion}];
        [~,~,stats]=anova1(vec_wt);
        [wt_c,wt_m]=multcompare(stats);
        wt_ff_fb_c{stims}{subregion}=wt_c;
        wt_ff_fb_m{stims}{subregion}=wt_m(:,1);
        wt_ff_fb_se{stims}{subregion}=wt_c(:,2);
        
        [~,~,stats]=anova1(vec_tw);
        [tw_c,tw_m]=multcompare(stats);
        tw_ff_fb_c{stims}{subregion}=tw_c;
        tw_ff_fb_m{stims}{subregion}=tw_m(:,1);
        tw_ff_fb_se{stims}{subregion}=tw_c(:,2);
    end
end

%% horzbar from anova
dat2plot_ff=[];
dat2plot_ff_se=[];
dat2plot_fb=[];
dat2plot_fb_se=[];
labels={'EC/EC-DG','EC-DG/DG','DG/DG-CA3','DG-CA3/CA3','CA3/CA3-CA1','CA3-CA1/CA1','CA1/CA1-EC','CA1-EC/EC'};
figure( 'Position', [100 100 1400 600])
for i=1:4
    dat2plot_ff=[dat2plot_ff;wt_ff_m(i,:);tw_ff_m(i,:)];
    dat2plot_ff_se=[dat2plot_ff_se;wt_ff_se(i,:);tw_ff_se(i,:)];
    dat2plot_fb=[dat2plot_fb;wt_fb_m(i,:);tw_fb_m(i,:)];
    dat2plot_fb_se=[dat2plot_fb_se;wt_fb_se(i,:);tw_fb_se(i,:)];
end

plot_horzbar_from_anova...
    (dat2plot_ff,dat2plot_ff_se,dat2plot_fb,dat2plot_fb_se,labels,ff_colors,fb_colors)
xlabel("FB Edges                  FF Edges")
xlim([-50,50])
xticks([-50:10:50])
ylim([0,9])
set(gca,'LineWidth',2,'TickLength',[0.02 0.02]);
set(gca,'FontSize',18)
% exportgraphics(gcf,".\three stim share figs\edges_per_subregion_comp_multcompare.png",'Resolution',1500)
exportgraphics(gcf,".\Temporal Analysis\Graph Analysis\edges_per_subregion_comp_multcompare.png",'Resolution',1500)

%% Get N

%out edges
[nostim_well_out_ff,nostim_well_out_fb,nostim_tunnel_out_ff,nostim_tunnel_out_fb]=...
    get_out_degrees(nostim_slope_G,nostim_datatable,nostim_data_fldr,nostim_temporal_fldr);

[hfs5_well_out_ff,hfs5_well_out_fb,hfs5_tunnel_out_ff,hfs5_tunnel_out_fb]=...
    get_out_degrees(hfs5_slope_G,hfs5_datatable,hfs5_data_fldr,hfs5_temporal_fldr);

[hfs40_well_out_ff,hfs40_well_out_fb,hfs40_tunnel_out_ff,hfs40_tunnel_out_fb]=...
    get_out_degrees(hfs40_slope_G,hfs40_datatable,hfs40_data_fldr,hfs40_temporal_fldr);
% in edges
[nostim_well_in_ff,nostim_well_in_fb,nostim_tunnel_in_ff,nostim_tunnel_in_fb]=...
    get_in_degrees(nostim_slope_G,nostim_datatable,nostim_data_fldr,nostim_temporal_fldr);

[hfs5_well_in_ff,hfs5_well_in_fb,hfs5_tunnel_in_ff,hfs5_tunnel_in_fb]=...
    get_in_degrees(hfs5_slope_G,hfs5_datatable,hfs5_data_fldr,hfs5_temporal_fldr);

[hfs40_well_in_ff,hfs40_well_in_fb,hfs40_tunnel_in_ff,hfs40_tunnel_in_fb]=...
    get_in_degrees(hfs40_slope_G,hfs40_datatable,hfs40_data_fldr,hfs40_temporal_fldr);

%nostim well
for i=1:4
    nostim_well_in_ff_n(i)=numel(nonzeros(cell2mat(nostim_well_in_ff{:,i})));
end

for i=1:4
    nostim_well_in_fb_n(i)=numel(nonzeros(cell2mat(nostim_well_in_fb{:,i})));
end

for i=1:4
    nostim_well_out_ff_n(i)=numel(nonzeros(cell2mat(nostim_well_out_ff{:,i})));
end

for i=1:4
    nostim_well_out_fb_n(i)=numel(nonzeros(cell2mat(nostim_well_out_fb{:,i})));
end

%nostim tunnel
for i=1:4
    nostim_tunnel_in_ff_n(i)=numel(nonzeros(cell2mat(nostim_tunnel_in_ff{:,i})));
end

for i=1:4
    nostim_tunnel_in_fb_n(i)=numel(nonzeros(cell2mat(nostim_tunnel_in_fb{:,i})));
end

for i=1:4
    nostim_tunnel_out_ff_n(i)=numel(nonzeros(cell2mat(nostim_tunnel_out_ff{:,i})));
end

for i=1:4
    nostim_tunnel_out_fb_n(i)=numel(nonzeros(cell2mat(nostim_tunnel_out_fb{:,i})));
end

%hfs5 well
for i=1:4
    hfs5_well_in_ff_n(i)=numel(nonzeros(cell2mat(hfs5_well_in_ff{:,i})));
end

for i=1:4
    hfs5_well_in_fb_n(i)=numel(nonzeros(cell2mat(hfs5_well_in_fb{:,i})));
end

for i=1:4
    hfs5_well_out_ff_n(i)=numel(nonzeros(cell2mat(hfs5_well_out_ff{:,i})));
end

for i=1:4
    hfs5_well_out_fb_n(i)=numel(nonzeros(cell2mat(hfs5_well_out_fb{:,i})));
end

%hfs40 tunnel
for i=1:4
    hfs5_tunnel_in_ff_n(i)=numel(nonzeros(cell2mat(hfs5_tunnel_in_ff{:,i})));
end

for i=1:4
    hfs5_tunnel_in_fb_n(i)=numel(nonzeros(cell2mat(hfs5_tunnel_in_fb{:,i})));
end

for i=1:4
    hfs5_tunnel_out_ff_n(i)=numel(nonzeros(cell2mat(hfs5_tunnel_out_ff{:,i})));
end

for i=1:4
    hfs5_tunnel_out_fb_n(i)=numel(nonzeros(cell2mat(hfs5_tunnel_out_fb{:,i})));
end

%hfs40 well
for i=1:4
    hfs40_well_in_ff_n(i)=numel(nonzeros(cell2mat(hfs40_well_in_ff{:,i})));
end

for i=1:4
    hfs40_well_in_fb_n(i)=numel(nonzeros(cell2mat(hfs40_well_in_fb{:,i})));
end

for i=1:4
    hfs40_well_out_ff_n(i)=numel(nonzeros(cell2mat(hfs40_well_out_ff{:,i})));
end

for i=1:4
    hfs40_well_out_fb_n(i)=numel(nonzeros(cell2mat(hfs40_well_out_fb{:,i})));
end

%hfs40 tunnel
for i=1:4
    hfs40_tunnel_in_ff_n(i)=numel(nonzeros(cell2mat(hfs40_tunnel_in_ff{:,i})));
end

for i=1:4
    hfs40_tunnel_in_fb_n(i)=numel(nonzeros(cell2mat(hfs40_tunnel_in_fb{:,i})));
end

for i=1:4
    hfs40_tunnel_out_ff_n(i)=numel(nonzeros(cell2mat(hfs40_tunnel_out_ff{:,i})));
end

for i=1:4
    hfs40_tunnel_out_fb_n(i)=numel(nonzeros(cell2mat(hfs40_tunnel_out_fb{:,i})));
end

% No stim
nostim_well_out_ff_count=cellfun(@sum,table2cell(nostim_well_out_ff(:,[1:4])));
nostim_well_out_fb_count=cellfun(@sum,table2cell(nostim_well_out_fb(:,[1:4])));
nostim_well_in_ff_count=cellfun(@sum,table2cell(nostim_well_in_ff(:,[1:4])));
nostim_well_in_fb_count=cellfun(@sum,table2cell(nostim_well_in_fb(:,[1:4])));

% 5 HFS
hfs5_well_out_ff_count=cellfun(@sum,table2cell(hfs5_well_out_ff(:,[1:4])));
hfs5_well_out_fb_count=cellfun(@sum,table2cell(hfs5_well_out_fb(:,[1:4])));
hfs5_well_in_ff_count=cellfun(@sum,table2cell(hfs5_well_in_ff(:,[1:4])));
hfs5_well_in_fb_count=cellfun(@sum,table2cell(hfs5_well_in_fb(:,[1:4])));

% 40 HFS
hfs40_well_out_ff_count=cellfun(@sum,table2cell(hfs40_well_out_ff(:,[1:4])));
hfs40_well_out_fb_count=cellfun(@sum,table2cell(hfs40_well_out_fb(:,[1:4])));
hfs40_well_in_ff_count=cellfun(@sum,table2cell(hfs40_well_in_ff(:,[1:4])));
hfs40_well_in_fb_count=cellfun(@sum,table2cell(hfs40_well_in_fb(:,[1:4])));

% plot N
dat2plot_ff=[];
dat2plot_fb=[];
ff_SE=[];
fb_SE=[];
labels={'EC In','EC Out','DG In','DG Out','CA3 In','CA3 Out','CA1 In','CA1 Out'};
figure( 'Position', [100 100 1400 600])
for i=1:4
%     dat2plot_ff=[dat2plot_ff;[nostim_well_in_ff_n(i)/9,hfs5_well_in_ff_n(i)/6,hfs40_well_in_ff_n(i)/6];...
%         [nostim_well_out_ff_n(i)/9,hfs5_well_out_ff_n(i)/6,hfs40_well_out_ff_n(i)/6]];
%     dat2plot_fb=[dat2plot_fb;[nostim_well_in_fb_n(i)/9,hfs5_well_in_fb_n(i)/6,hfs40_well_in_fb_n(i)/6];...
%         [nostim_well_out_fb_n(i)/9,hfs5_well_out_fb_n(i)/6,hfs40_well_out_fb_n(i)/6]];

%     dat2plot_ff=[dat2plot_ff;[mean(nostim_well_in_ff_count(:,i)),mean(hfs5_well_in_ff_count(:,i)),mean(hfs40_well_in_ff_count(:,i))];...
%         [mean(nostim_well_out_ff_count(:,i)),mean(hfs5_well_out_ff_count(:,i)),mean(hfs40_well_out_ff_count(:,i))]];
%     dat2plot_fb=[dat2plot_fb;[mean(nostim_well_in_fb_count(:,i)),mean(hfs5_well_in_fb_count(:,i)),mean(hfs40_well_in_fb_count(:,i))];...
%         [mean(nostim_well_out_fb_count(:,i)),mean(hfs5_well_out_fb_count(:,i)),mean(hfs40_well_out_fb_count(:,i))]];
%     ff_SE=[ff_SE;[std(nostim_well_in_ff_count(:,i))/sqrt(length(nostim_well_in_ff_count(:,i))),...
%         std(hfs5_well_in_ff_count(:,i))/sqrt(length(hfs5_well_in_ff_count(:,i))),...
%         std(hfs40_well_in_ff_count(:,i))/sqrt(length(hfs40_well_in_ff_count(:,i)))];...
%         [std(nostim_well_out_ff_count(:,i))/sqrt(length(nostim_well_out_ff_count(:,i))),...
%         std(hfs5_well_out_ff_count(:,i))/sqrt(length(hfs5_well_out_ff_count(:,i))),...
%         std(hfs40_well_out_ff_count(:,i))/sqrt(length(hfs40_well_out_ff_count(:,i)))]];
%     fb_SE=[fb_SE;[std(nostim_well_in_fb_count(:,i))/sqrt(length(nostim_well_in_fb_count(:,i))),...
%         std(hfs5_well_in_fb_count(:,i))/sqrt(length(hfs5_well_in_fb_count(:,i))),...
%         std(hfs40_well_in_fb_count(:,i))/sqrt(length(hfs40_well_in_fb_count(:,i)))];...
%         [std(nostim_well_out_fb_count(:,i))/sqrt(length(nostim_well_out_fb_count(:,i))),...
%         std(hfs5_well_out_fb_count(:,i))/sqrt(length(hfs5_well_out_fb_count(:,i))),...
%         std(hfs40_well_out_fb_count(:,i))/sqrt(length(hfs40_well_out_fb_count(:,i)))]];
    % Anova
    [p,~,stats]=anova1([(nostim_well_in_ff_count(:,i));(hfs5_well_in_ff_count(:,i));(hfs40_well_in_ff_count(:,i))],...
        [repmat("nostim",length(nostim_well_in_ff_count(:,i)),1);repmat("hfs5",length(hfs5_well_in_ff_count(:,i)),1);...
        repmat("hfs40",length(hfs40_well_in_ff_count(:,i)),1)]);
    [ff_in_c,ff_in_m]=multcompare(stats);
    [~,~,stats]=anova1([(nostim_well_in_fb_count(:,i));(hfs5_well_in_fb_count(:,i));(hfs40_well_in_fb_count(:,i))],...
        [repmat("nostim",length(nostim_well_in_fb_count(:,i)),1);repmat("hfs5",length(hfs5_well_in_fb_count(:,i)),1);...
        repmat("hfs40",length(hfs40_well_in_fb_count(:,i)),1)]);
    [fb_in_c,fb_in_m]=multcompare(stats);
    [~,~,stats]=anova1([(nostim_well_out_ff_count(:,i));(hfs5_well_out_ff_count(:,i));(hfs40_well_out_ff_count(:,i))],...
        [repmat("nostim",length(nostim_well_out_ff_count(:,i)),1);repmat("hfs5",length(hfs5_well_out_ff_count(:,i)),1);...
        repmat("hfs40",length(hfs40_well_out_ff_count(:,i)),1)]);
    [ff_out_c,ff_out_m]=multcompare(stats);
    [~,~,stats]=anova1([(nostim_well_out_fb_count(:,i));(hfs5_well_out_fb_count(:,i));(hfs40_well_out_fb_count(:,i))],...
        [repmat("nostim",length(nostim_well_out_fb_count(:,i)),1);repmat("hfs5",length(hfs5_well_out_fb_count(:,i)),1);...
        repmat("hfs40",length(hfs40_well_out_fb_count(:,i)),1)]);
    [fb_out_c,fb_out_m]=multcompare(stats);

    dat2plot_ff=[dat2plot_ff;ff_in_m(:,1)';ff_out_m(:,1)'];
    dat2plot_fb=[dat2plot_fb;fb_in_m(:,1)';fb_out_m(:,1)'];
    ff_SE=[ff_SE;ff_in_m(:,1)';ff_out_m(:,1)'];
    fb_SE=[fb_SE;fb_in_m(:,2)';fb_out_m(:,2)'];
end

% plot_horzbar_from_anova...
%     (dat2plot_ff,zeros(size(dat2plot_ff)),dat2plot_fb,zeros(size(dat2plot_fb)),labels,ff_colors,fb_colors)
plot_horzbar_from_anova...
    (dat2plot_ff,ff_SE,dat2plot_fb,fb_SE,labels,ff_colors,fb_colors)
xlabel("FB Edges                  FF Edges")
% xlim([-15,15])
% xticks([-15:3:15])
xlim([-60,60])
xticks([-60:10:60])
ylim([0,9])
set(gca,'LineWidth',2,'TickLength',[0.02 0.02]);
set(gca,'FontSize',18)
exportgraphics(gcf,".\Temporal Analysis\Graph Analysis\normalized_number_edges_per_array.png",'Resolution',1500)

%% Average weight slopes
%get weights, well is out edges for wells, tunnel is out edges for
%tunnels. Therefore well_ff is out, well_fb is in, tunnel_ff is in,
%tunnel_fb is in.
[nostim_well_ff,nostim_well_fb,nostim_tunnel_ff,nostim_tunnel_fb]=...
    get_weights(nostim_slope_G,nostim_datatable,nostim_data_fldr,nostim_temporal_fldr);
[hfs5_well_ff,hfs5_well_fb,hfs5_tunnel_ff,hfs5_tunnel_fb]=...
    get_weights(hfs5_slope_G,hfs5_datatable,hfs5_data_fldr,hfs5_temporal_fldr);
[hfs40_well_ff,hfs40_well_fb,hfs40_tunnel_ff,hfs40_tunnel_fb]=...
    get_weights(hfs40_slope_G,hfs40_datatable,hfs40_data_fldr,hfs40_temporal_fldr);

wt_ff_m=[];
wt_ff_se=[];
wt_fb_m=[];
wt_fb_se=[];
tw_ff_m=[];
tw_ff_se=[];
tw_fb_m=[];
tw_fb_se=[];
% anova
for i=1:4
    [~,~,stats]=anova1(...
        [vertcat(nostim_well_ff{:,i}{:});vertcat(hfs5_well_ff{:,i}{:});vertcat(hfs40_well_ff{:,i}{:})],...
        [repmat("nostim",length(vertcat(nostim_well_ff{:,i}{:})),1); ...
        repmat("hfs5",length(vertcat(hfs5_well_ff{:,i}{:})),1);...
        repmat("hfs40",length(vertcat(hfs40_well_ff{:,i}{:})),1)]);
    [c,m]=multcompare(stats);
    
    wt_ff_c{i}=c;
    wt_ff_m=[wt_ff_m;m(:,1)'];
    wt_ff_se=[wt_ff_se;m(:,2)'];

    [~,~,stats]=anova1(...
        [vertcat(nostim_well_fb{:,i}{:});vertcat(hfs5_well_fb{:,i}{:});vertcat(hfs40_well_fb{:,i}{:})],...
        [repmat("nostim",length(vertcat(nostim_well_fb{:,i}{:})),1); ...
        repmat("hfs5",length(vertcat(hfs5_well_fb{:,i}{:})),1);...
        repmat("hfs40",length(vertcat(hfs40_well_fb{:,i}{:})),1)]);
    [c,m]=multcompare(stats);

    wt_fb_c{i}=c;
    wt_fb_m=[wt_fb_m;m(:,1)'];
    wt_fb_se=[wt_fb_se;m(:,2)'];

    [~,~,stats]=anova1(...
        [vertcat(nostim_tunnel_ff{:,i}{:});vertcat(hfs5_tunnel_ff{:,i}{:});vertcat(hfs40_tunnel_ff{:,i}{:})],...
        [repmat("nostim",length(vertcat(nostim_tunnel_ff{:,i}{:})),1); ...
        repmat("hfs5",length(vertcat(hfs5_tunnel_ff{:,i}{:})),1);...
        repmat("hfs40",length(vertcat(hfs40_tunnel_ff{:,i}{:})),1)]);
    [c,m]=multcompare(stats);
    
    tw_ff_c{i}=c;
    tw_ff_m=[tw_ff_m;m(:,1)'];
    tw_ff_se=[tw_ff_se;m(:,2)'];

    [~,~,stats]=anova1(...
        [vertcat(nostim_tunnel_fb{:,i}{:});vertcat(hfs5_tunnel_fb{:,i}{:});vertcat(hfs40_tunnel_fb{:,i}{:})],...
        [repmat("nostim",length(vertcat(nostim_tunnel_fb{:,i}{:})),1); ...
        repmat("hfs5",length(vertcat(hfs5_tunnel_fb{:,i}{:})),1);...
        repmat("hfs40",length(vertcat(hfs40_tunnel_fb{:,i}{:})),1)]);
    [c,m]=multcompare(stats);

    tw_fb_c{i}=c;
    tw_fb_m=[tw_fb_m;m(:,1)'];
    tw_fb_se=[tw_fb_se;m(:,2)'];

    %in ff fb anova
    v=[vertcat(nostim_well_ff{:,i}{:});vertcat(nostim_well_fb{:,i}{:})];
    g=[repmat("ff",length(vertcat(nostim_well_ff{:,i}{:})),1);...
        repmat("fb",length(vertcat(nostim_well_fb{:,i}{:})),1)];
    [~,~,stats]=anova1(v,g);
    [c,~]=multcompare(stats);
    ff_fb_in_c{i,1}=c;
    
    v=[vertcat(hfs5_well_ff{:,i}{:});vertcat(hfs5_well_fb{:,i}{:})];
    g=[repmat("ff",length(vertcat(hfs5_well_ff{:,i}{:})),1);...
        repmat("fb",length(vertcat(hfs5_well_fb{:,i}{:})),1)];
    [~,~,stats]=anova1(v,g);
    [c,~]=multcompare(stats);
    ff_fb_in_c{i,2}=c;
    
    v=[vertcat(hfs40_well_ff{:,i}{:});vertcat(hfs40_well_fb{:,i}{:})];
    g=[repmat("ff",length(vertcat(hfs40_well_ff{:,i}{:})),1);...
        repmat("fb",length(vertcat(hfs40_well_fb{:,i}{:})),1)];
    [~,~,stats]=anova1(v,g);
    [c,~]=multcompare(stats);
    ff_fb_in_c{i,3}=c;
    
    %out ff fb anova
    v=[vertcat(nostim_tunnel_ff{:,i}{:});vertcat(nostim_tunnel_fb{:,i}{:})];
    g=[repmat("ff",length(vertcat(nostim_tunnel_ff{:,i}{:})),1);...
        repmat("fb",length(vertcat(nostim_tunnel_fb{:,i}{:})),1)];
    [~,~,stats]=anova1(v,g);
    [c,~]=multcompare(stats);
    ff_fb_out_c{i,1}=c;
    
    v=[vertcat(hfs5_tunnel_ff{:,i}{:});vertcat(hfs5_tunnel_fb{:,i}{:})];
    g=[repmat("ff",length(vertcat(hfs5_tunnel_ff{:,i}{:})),1);...
        repmat("fb",length(vertcat(hfs5_tunnel_fb{:,i}{:})),1)];
    [~,~,stats]=anova1(v,g);
    [c,~]=multcompare(stats);
    ff_fb_out_c{i,2}=c;
    
    v=[vertcat(hfs40_tunnel_ff{:,i}{:});vertcat(hfs40_tunnel_fb{:,i}{:})];
    g=[repmat("ff",length(vertcat(hfs40_tunnel_ff{:,i}{:})),1);...
        repmat("fb",length(vertcat(hfs40_tunnel_fb{:,i}{:})),1)];
    [~,~,stats]=anova1(v,g);
    [c,~]=multcompare(stats);
    ff_fb_out_c{i,3}=c;

    effect_out{i}=ff_fb_meanEffect(vertcat(nostim_well_ff{:,i}{:}),vertcat(nostim_well_fb{:,i}{:}),...
        vertcat(hfs5_well_ff{:,i}{:}),vertcat(hfs5_well_fb{:,i}{:}),...
        vertcat(hfs40_well_ff{:,i}{:}),vertcat(hfs40_well_fb{:,i}{:}));
    effect_in{i}=ff_fb_meanEffect(vertcat(nostim_tunnel_ff{:,i}{:}),vertcat(nostim_tunnel_fb{:,i}{:}),...
        vertcat(hfs5_tunnel_ff{:,i}{:}),vertcat(hfs5_tunnel_fb{:,i}{:}),...
        vertcat(hfs40_tunnel_ff{:,i}{:}),vertcat(hfs40_tunnel_fb{:,i}{:}));    
end

% horzbar from anova slope
dat2plot_ff=[];
dat2plot_ff_se=[];
dat2plot_fb=[];
dat2plot_fb_se=[];
labels={'EC In','EC Out','DG In','DG Out','CA3 In','CA3 Out','CA1 In','CA1 Out'};
figure( 'Position', [100 100 1400 600])
for i=1:4
    dat2plot_ff=[dat2plot_ff;tw_ff_m(i,:);wt_ff_m(i,:)];
    dat2plot_ff_se=[dat2plot_ff_se;tw_ff_se(i,:);wt_ff_se(i,:)];
    dat2plot_fb=[dat2plot_fb;tw_fb_m(i,:);wt_fb_m(i,:)];
    dat2plot_fb_se=[dat2plot_fb_se;tw_fb_se(i,:);wt_fb_se(i,:)];
end

plot_horzbar_from_anova...
    (dat2plot_ff,dat2plot_ff_se,dat2plot_fb,dat2plot_fb_se,labels,ff_colors,fb_colors)
xlabel("FB Slope Weight                  FF Slope Weight")
xlim([-2,2])
xticks([-5:0.25:5])
ylim([0,9])
set(gca,'LineWidth',2,'TickLength',[0.02 0.02]);
set(gca,'FontSize',18)
% exportgraphics(gcf,".\three stim share figs\edges_per_subregion_comp_multcompare.png",'Resolution',1500)
exportgraphics(gcf,".\Temporal Analysis\Graph Analysis\slope_weights.png",'Resolution',1500)
%% FF/FB% degrees

ff_fb_perc=(dat2plot_ff./dat2plot_fb)*100;

% ff_fb_SD=dat2plot

ff_fb_SD=[];
ff_fb_lengths=[];
for regi=1:4
    
    %in
    %ff wells in
    nostim_well_ff_degs=nonzeros(cell2mat(nostim_tunnel_ff{:,regi}));
    hfs5_well_ff_degs=nonzeros(cell2mat(hfs5_tunnel_ff{:,regi}));
    hfs40_well_ff_degs=nonzeros(cell2mat(hfs40_tunnel_ff{:,regi}));
    node_degrees_ff_well_m=[mean(nostim_well_ff_degs),mean(hfs5_well_ff_degs),mean(hfs40_well_ff_degs)];
    node_degrees_ff_well_sd=[std(nostim_well_ff_degs),std(hfs5_well_ff_degs),std(hfs40_well_ff_degs)];
    ff_out_lengths=[length(nostim_well_ff_degs),length(hfs5_well_ff_degs),length(hfs40_well_ff_degs)];

    %fb wells in
    nostim_well_fb_degs=nonzeros(cell2mat(nostim_tunnel_fb{:,regi}));
    hfs5_well_fb_degs=nonzeros(cell2mat(hfs5_tunnel_fb{:,regi}));
    hfs40_well_fb_degs=nonzeros(cell2mat(hfs40_tunnel_fb{:,regi}));
    node_degrees_fb_well_m=[mean(nostim_well_fb_degs),mean(hfs5_well_fb_degs),mean(hfs40_well_fb_degs)];
    node_degrees_fb_well_sd=[std(nostim_well_fb_degs),std(hfs5_well_fb_degs),std(hfs40_well_fb_degs)];
    fb_out_lengths=[length(nostim_well_fb_degs),length(hfs5_well_fb_degs),length(hfs40_well_fb_degs)];
    
    %out
    %ff out
    nostim_tunnel_ff_degs=nonzeros(cell2mat(nostim_well_ff{:,regi}));
    hfs5_tunnel_ff_degs=nonzeros(cell2mat(hfs5_well_ff{:,regi}));
    hfs40_tunnel_ff_degs=nonzeros(cell2mat(hfs40_well_ff{:,regi}));
    node_degrees_ff_tunnel_m=[mean(nostim_tunnel_ff_degs),mean(hfs5_tunnel_ff_degs),mean(hfs40_tunnel_ff_degs)];
    node_degrees_ff_tunnel_sd=[std(nostim_tunnel_ff_degs),std(hfs5_tunnel_ff_degs), std(hfs40_tunnel_ff_degs)];
    ff_in_lengths=[length(nostim_tunnel_ff_degs),length(hfs5_tunnel_ff_degs),length(hfs40_tunnel_ff_degs)];

    %fb out
    nostim_tunnel_fb_degs=nonzeros(cell2mat(nostim_well_fb{:,regi}));
    hfs5_tunnel_fb_degs=nonzeros(cell2mat(hfs5_well_fb{:,regi}));
    hfs40_tunnel_fb_degs=nonzeros(cell2mat(hfs40_well_fb{:,regi}));
    node_degrees_fb_tunnel_m=[mean(nostim_tunnel_fb_degs),mean(hfs5_tunnel_fb_degs),mean(hfs40_tunnel_fb_degs)];
    node_degrees_fb_tunnel_sd=[std(nostim_tunnel_fb_degs),std(hfs5_tunnel_fb_degs),std(hfs40_tunnel_fb_degs)];
    fb_in_lengths=[length(nostim_tunnel_fb_degs),length(hfs5_tunnel_fb_degs),length(hfs40_tunnel_fb_degs)];

    ff_fb_SD=[ff_fb_SD;((node_degrees_ff_well_m./node_degrees_fb_well_m)*100).*...
        sqrt((node_degrees_ff_well_sd./node_degrees_fb_well_sd).^2+...
        (node_degrees_fb_well_sd./node_degrees_fb_well_sd).^2);...
        ((node_degrees_ff_tunnel_m./node_degrees_fb_tunnel_m)*100).*...
        sqrt((node_degrees_ff_tunnel_m./node_degrees_fb_tunnel_m).^2+...
        (node_degrees_fb_tunnel_sd./node_degrees_fb_tunnel_sd).^2)];
    ff_fb_lengths=[ff_fb_lengths;(ff_in_lengths+fb_in_lengths);(ff_out_lengths+fb_out_lengths)];
end

labels={'EC-DG In','EC-DG Out','DG-CA3 In','DG-CA3 Out','CA3-CA1 In','CA3-CA1 Out','CA1-EC In','CA1-EC Out'};
regs=categorical(labels);
regs=reordercats(regs,labels);

figure( 'Position', [100 100 700 600])
b=bar(regs,ff_fb_perc);
b(1).BaseValue=100;
b(1).BaseLine.LineStyle = "--";
b(1).BaseLine.LineWidth=2;

ff_fb_SE=ff_fb_SD./sqrt(ff_fb_lengths);
hold on
% grouped error bars
[ngroups,nbars]=size(ff_fb_perc);
x=[];
for j=1:nbars
    x(j,:)=b(j).XEndPoints;
end
% errorbar(x',ff_fb_perc,ff_fb_SE,'k','LineStyle','none')
hold off
ylabel("% FF/FB")
ylim([50,180])
set(gca,'LineWidth',2,'TickLength',[0.02 0.02]);
set(gca,'FontSize',18)
exportgraphics(gcf,".\Temporal Analysis\Graph Analysis\slope_ff_fb_perc.png",'Resolution',1500)
%% Cohen's D Bar Graphs
CD_in=[];
CD_out=[];
CD_in_CE_neg=[];
CD_in_CE_pos=[];
CD_out_CE_neg=[];
CD_out_CE_pos=[];
labels={'EC In','EC Out','DG In','DG Out','CA3 In','CA3 Out','CA1 In','CA1 Out'};

for i=1:4
    CD_in=[CD_in;[effect_in{i}{1,1}{1}.Effect,effect_in{i}{1,2}{1}.Effect,effect_in{i}{1,3}{1}.Effect]];
    CD_in_CE_neg=[CD_in_CE_neg;[effect_in{i}{1,1}{1}.ConfidenceIntervals(1),effect_in{i}{1,2}{1}.ConfidenceIntervals(1),effect_in{i}{1,3}{1}.ConfidenceIntervals(1)]];
    CD_in_CE_pos=[CD_in_CE_pos;[effect_in{i}{1,1}{1}.ConfidenceIntervals(2),effect_in{i}{1,2}{1}.ConfidenceIntervals(2),effect_in{i}{1,3}{1}.ConfidenceIntervals(2)]];
    CD_out=[CD_out;[effect_out{i}{1,1}{1}.Effect,effect_out{i}{1,2}{1}.Effect,effect_out{i}{1,3}{1}.Effect]];
    CD_out_CE_neg=[CD_out_CE_neg;[effect_out{i}{1,1}{1}.ConfidenceIntervals(1),effect_out{i}{1,2}{1}.ConfidenceIntervals(1),effect_out{i}{1,3}{1}.ConfidenceIntervals(1)]];
    CD_out_CE_pos=[CD_out_CE_pos;[effect_out{i}{1,1}{1}.ConfidenceIntervals(2),effect_out{i}{1,2}{1}.ConfidenceIntervals(2),effect_out{i}{1,3}{1}.ConfidenceIntervals(2)]];
end

CD_m=interleave_rows(CD_in,CD_out);
CD_CE_neg=interleave_rows(CD_in_CE_neg,CD_out_CE_neg);
CD_CE_pos=interleave_rows(CD_in_CE_pos,CD_out_CE_pos);

figure( 'Position', [100 100 700 600])
plot_horzbar_CohenD...
    (CD_m,CD_CE_neg,CD_CE_pos,labels);
xlabel("Slope Directional Mean Effect Size")
xlim([-3,3])
ylim([0,9])
set(gca,'LineWidth',2,'TickLength',[0.02 0.02]);
set(gca,'FontSize',18)
exportgraphics(gcf,".\Temporal Analysis\Graph Analysis\slope_cohenD.png",'Resolution',1500)

%% Average weight rsq
%get weights
[nostim_well_ff,nostim_well_fb,nostim_tunnel_ff,nostim_tunnel_fb]=...
    get_weights(nostim_rsq_G,nostim_datatable,nostim_data_fldr,nostim_temporal_fldr);
[hfs5_well_ff,hfs5_well_fb,hfs5_tunnel_ff,hfs5_tunnel_fb]=...
    get_weights(hfs5_rsq_G,hfs5_datatable,hfs5_data_fldr,hfs5_temporal_fldr);
[hfs40_well_ff,hfs40_well_fb,hfs40_tunnel_ff,hfs40_tunnel_fb]=...
    get_weights(hfs40_rsq_G,hfs40_datatable,hfs40_data_fldr,hfs40_temporal_fldr);

wt_ff_m=[];
wt_ff_se=[];
wt_fb_m=[];
wt_fb_se=[];
tw_ff_m=[];
tw_ff_se=[];
tw_fb_m=[];
tw_fb_se=[];
% anova
for i=1:4
    [~,~,stats]=anova1(...
        [vertcat(nostim_well_ff{:,i}{:});vertcat(hfs5_well_ff{:,i}{:});vertcat(hfs40_well_ff{:,i}{:})],...
        [repmat("nostim",length(vertcat(nostim_well_ff{:,i}{:})),1); ...
        repmat("hfs5",length(vertcat(hfs5_well_ff{:,i}{:})),1);...
        repmat("hfs40",length(vertcat(hfs40_well_ff{:,i}{:})),1)]);
    [c,m]=multcompare(stats);
    
    wt_ff_c{i}=c;
    wt_ff_m=[wt_ff_m;m(:,1)'];
    wt_ff_se=[wt_ff_se;m(:,2)'];

    [~,~,stats]=anova1(...
        [vertcat(nostim_well_fb{:,i}{:});vertcat(hfs5_well_fb{:,i}{:});vertcat(hfs40_well_fb{:,i}{:})],...
        [repmat("nostim",length(vertcat(nostim_well_fb{:,i}{:})),1); ...
        repmat("hfs5",length(vertcat(hfs5_well_fb{:,i}{:})),1);...
        repmat("hfs40",length(vertcat(hfs40_well_fb{:,i}{:})),1)]);
    [c,m]=multcompare(stats);

    wt_fb_c{i}=c;
    wt_fb_m=[wt_fb_m;m(:,1)'];
    wt_fb_se=[wt_fb_se;m(:,2)'];

    [~,~,stats]=anova1(...
        [vertcat(nostim_tunnel_ff{:,i}{:});vertcat(hfs5_tunnel_ff{:,i}{:});vertcat(hfs40_tunnel_ff{:,i}{:})],...
        [repmat("nostim",length(vertcat(nostim_tunnel_ff{:,i}{:})),1); ...
        repmat("hfs5",length(vertcat(hfs5_tunnel_ff{:,i}{:})),1);...
        repmat("hfs40",length(vertcat(hfs40_tunnel_ff{:,i}{:})),1)]);
    [c,m]=multcompare(stats);
    
    tw_ff_c{i}=c;
    tw_ff_m=[tw_ff_m;m(:,1)'];
    tw_ff_se=[tw_ff_se;m(:,2)'];

    [~,~,stats]=anova1(...
        [vertcat(nostim_tunnel_fb{:,i}{:});vertcat(hfs5_tunnel_fb{:,i}{:});vertcat(hfs40_tunnel_fb{:,i}{:})],...
        [repmat("nostim",length(vertcat(nostim_tunnel_fb{:,i}{:})),1); ...
        repmat("hfs5",length(vertcat(hfs5_tunnel_fb{:,i}{:})),1);...
        repmat("hfs40",length(vertcat(hfs40_tunnel_fb{:,i}{:})),1)]);
    [c,m]=multcompare(stats);

    tw_fb_c{i}=c;
    tw_fb_m=[tw_fb_m;m(:,1)'];
    tw_fb_se=[tw_fb_se;m(:,2)'];

    %in ff fb anova
    v=[vertcat(nostim_well_ff{:,i}{:});vertcat(nostim_well_fb{:,i}{:})];
    g=[repmat("ff",length(vertcat(nostim_well_ff{:,i}{:})),1);...
        repmat("fb",length(vertcat(nostim_well_fb{:,i}{:})),1)];
    [~,~,stats]=anova1(v,g);
    [c,~]=multcompare(stats);
    ff_fb_in_c{i,1}=c;
    
    v=[vertcat(hfs5_well_ff{:,i}{:});vertcat(hfs5_well_fb{:,i}{:})];
    g=[repmat("ff",length(vertcat(hfs5_well_ff{:,i}{:})),1);...
        repmat("fb",length(vertcat(hfs5_well_fb{:,i}{:})),1)];
    [~,~,stats]=anova1(v,g);
    [c,~]=multcompare(stats);
    ff_fb_in_c{i,2}=c;
    
    v=[vertcat(hfs40_well_ff{:,i}{:});vertcat(hfs40_well_fb{:,i}{:})];
    g=[repmat("ff",length(vertcat(hfs40_well_ff{:,i}{:})),1);...
        repmat("fb",length(vertcat(hfs40_well_fb{:,i}{:})),1)];
    [~,~,stats]=anova1(v,g);
    [c,~]=multcompare(stats);
    ff_fb_in_c{i,3}=c;
    
    %out ff fb anova
    v=[vertcat(nostim_tunnel_ff{:,i}{:});vertcat(nostim_tunnel_fb{:,i}{:})];
    g=[repmat("ff",length(vertcat(nostim_tunnel_ff{:,i}{:})),1);...
        repmat("fb",length(vertcat(nostim_tunnel_fb{:,i}{:})),1)];
    [~,~,stats]=anova1(v,g);
    [c,~]=multcompare(stats);
    ff_fb_out_c{i,1}=c;
    
    v=[vertcat(hfs5_tunnel_ff{:,i}{:});vertcat(hfs5_tunnel_fb{:,i}{:})];
    g=[repmat("ff",length(vertcat(hfs5_tunnel_ff{:,i}{:})),1);...
        repmat("fb",length(vertcat(hfs5_tunnel_fb{:,i}{:})),1)];
    [~,~,stats]=anova1(v,g);
    [c,~]=multcompare(stats);
    ff_fb_out_c{i,2}=c;
    
    v=[vertcat(hfs40_tunnel_ff{:,i}{:});vertcat(hfs40_tunnel_fb{:,i}{:})];
    g=[repmat("ff",length(vertcat(hfs40_tunnel_ff{:,i}{:})),1);...
        repmat("fb",length(vertcat(hfs40_tunnel_fb{:,i}{:})),1)];
    [~,~,stats]=anova1(v,g);
    [c,~]=multcompare(stats);
    ff_fb_out_c{i,3}=c;

    effect_in{i}=ff_fb_meanEffect(vertcat(nostim_well_ff{:,i}{:}),vertcat(nostim_well_fb{:,i}{:}),...
        vertcat(hfs5_well_ff{:,i}{:}),vertcat(hfs5_well_fb{:,i}{:}),...
        vertcat(hfs40_well_ff{:,i}{:}),vertcat(hfs40_well_fb{:,i}{:}));
    effect_out{i}=ff_fb_meanEffect(vertcat(nostim_tunnel_ff{:,i}{:}),vertcat(nostim_tunnel_fb{:,i}{:}),...
        vertcat(hfs5_tunnel_ff{:,i}{:}),vertcat(hfs5_tunnel_fb{:,i}{:}),...
        vertcat(hfs40_tunnel_ff{:,i}{:}),vertcat(hfs40_tunnel_fb{:,i}{:}));      
end

% horzbar from anova slope
dat2plot_ff=[];
dat2plot_ff_se=[];
dat2plot_fb=[];
dat2plot_fb_se=[];
labels={'EC In','EC Out','DG In','DG Out','CA3 In','CA3 Out','CA1 In','CA1 Out'};
figure( 'Position', [100 100 1400 600])
for i=1:4
    dat2plot_ff=[dat2plot_ff;tw_ff_m(i,:);wt_ff_m(i,:)];
    dat2plot_ff_se=[dat2plot_ff_se;tw_ff_se(i,:);wt_ff_se(i,:)];
    dat2plot_fb=[dat2plot_fb;tw_fb_m(i,:);wt_fb_m(i,:)];
    dat2plot_fb_se=[dat2plot_fb_se;tw_fb_se(i,:);wt_fb_se(i,:)];
end

plot_horzbar_from_anova...
    (dat2plot_ff,dat2plot_ff_se,dat2plot_fb,dat2plot_fb_se,labels,ff_colors,fb_colors)
xlabel("FB R^2 Weight                  FF R^2 Weight")
xlim([-0.7,0.7])
xticks([-1:.1:1])
ylim([0,9])
set(gca,'LineWidth',2,'TickLength',[0.02 0.02]);
set(gca,'FontSize',18)
% exportgraphics(gcf,".\three stim share figs\edges_per_subregion_comp_multcompare.png",'Resolution',1500)
exportgraphics(gcf,".\Temporal Analysis\Graph Analysis\rsq_weights.png",'Resolution',1500)
%% FF/FB% degrees

ff_fb_perc=(dat2plot_ff./dat2plot_fb)*100;

% ff_fb_SD=dat2plot

ff_fb_SD=[];
ff_fb_lengths=[];
for regi=1:4
    
    %in
    %ff wells in
    nostim_well_ff_degs=nonzeros(cell2mat(nostim_tunnel_ff{:,regi}));
    hfs5_well_ff_degs=nonzeros(cell2mat(hfs5_tunnel_ff{:,regi}));
    hfs40_well_ff_degs=nonzeros(cell2mat(hfs40_tunnel_ff{:,regi}));
    node_degrees_ff_well_m=[mean(nostim_well_ff_degs),mean(hfs5_well_ff_degs),mean(hfs40_well_ff_degs)];
    node_degrees_ff_well_sd=[std(nostim_well_ff_degs),std(hfs5_well_ff_degs),std(hfs40_well_ff_degs)];
    ff_out_lengths=[length(nostim_well_ff_degs),length(hfs5_well_ff_degs),length(hfs40_well_ff_degs)];

    %fb wells in
    nostim_well_fb_degs=nonzeros(cell2mat(nostim_tunnel_fb{:,regi}));
    hfs5_well_fb_degs=nonzeros(cell2mat(hfs5_tunnel_fb{:,regi}));
    hfs40_well_fb_degs=nonzeros(cell2mat(hfs40_tunnel_fb{:,regi}));
    node_degrees_fb_well_m=[mean(nostim_well_fb_degs),mean(hfs5_well_fb_degs),mean(hfs40_well_fb_degs)];
    node_degrees_fb_well_sd=[std(nostim_well_fb_degs),std(hfs5_well_fb_degs),std(hfs40_well_fb_degs)];
    fb_out_lengths=[length(nostim_well_fb_degs),length(hfs5_well_fb_degs),length(hfs40_well_fb_degs)];
    
    %out
    %ff out
    nostim_tunnel_ff_degs=nonzeros(cell2mat(nostim_well_ff{:,regi}));
    hfs5_tunnel_ff_degs=nonzeros(cell2mat(hfs5_well_ff{:,regi}));
    hfs40_tunnel_ff_degs=nonzeros(cell2mat(hfs40_well_ff{:,regi}));
    node_degrees_ff_tunnel_m=[mean(nostim_tunnel_ff_degs),mean(hfs5_tunnel_ff_degs),mean(hfs40_tunnel_ff_degs)];
    node_degrees_ff_tunnel_sd=[std(nostim_tunnel_ff_degs),std(hfs5_tunnel_ff_degs), std(hfs40_tunnel_ff_degs)];
    ff_in_lengths=[length(nostim_tunnel_ff_degs),length(hfs5_tunnel_ff_degs),length(hfs40_tunnel_ff_degs)];

    %fb out
    nostim_tunnel_fb_degs=nonzeros(cell2mat(nostim_well_fb{:,regi}));
    hfs5_tunnel_fb_degs=nonzeros(cell2mat(hfs5_well_fb{:,regi}));
    hfs40_tunnel_fb_degs=nonzeros(cell2mat(hfs40_well_fb{:,regi}));
    node_degrees_fb_tunnel_m=[mean(nostim_tunnel_fb_degs),mean(hfs5_tunnel_fb_degs),mean(hfs40_tunnel_fb_degs)];
    node_degrees_fb_tunnel_sd=[std(nostim_tunnel_fb_degs),std(hfs5_tunnel_fb_degs),std(hfs40_tunnel_fb_degs)];
    fb_in_lengths=[length(nostim_tunnel_fb_degs),length(hfs5_tunnel_fb_degs),length(hfs40_tunnel_fb_degs)];

    ff_fb_SD=[ff_fb_SD;((node_degrees_ff_well_m./node_degrees_fb_well_m)*100).*...
        sqrt((node_degrees_ff_well_sd./node_degrees_fb_well_sd).^2+...
        (node_degrees_fb_well_sd./node_degrees_fb_well_sd).^2);...
        ((node_degrees_ff_tunnel_m./node_degrees_fb_tunnel_m)*100).*...
        sqrt((node_degrees_ff_tunnel_m./node_degrees_fb_tunnel_m).^2+...
        (node_degrees_fb_tunnel_sd./node_degrees_fb_tunnel_sd).^2)];
    ff_fb_lengths=[ff_fb_lengths;(ff_in_lengths+fb_in_lengths);(ff_out_lengths+fb_out_lengths)];
end

labels={'EC-DG In','EC-DG Out','DG-CA3 In','DG-CA3 Out','CA3-CA1 In','CA3-CA1 Out','CA1-EC In','CA1-EC Out'};
regs=categorical(labels);
regs=reordercats(regs,labels);

figure( 'Position', [100 100 700 600])
b=bar(regs,ff_fb_perc);
b(1).BaseValue=100;
b(1).BaseLine.LineStyle = "--";
b(1).BaseLine.LineWidth=2;

ff_fb_SE=ff_fb_SD./sqrt(ff_fb_lengths);
hold on
% grouped error bars
[ngroups,nbars]=size(ff_fb_perc);
x=[];
for j=1:nbars
    x(j,:)=b(j).XEndPoints;
end
% errorbar(x',ff_fb_perc,ff_fb_SE,'k','LineStyle','none')
hold off
ylabel("% FF/FB")
set(gca,'LineWidth',2,'TickLength',[0.02 0.02]);
set(gca,'FontSize',18)
exportgraphics(gcf,".\Temporal Analysis\Graph Analysis\rsq_ff_fb_perc.png",'Resolution',1500)
%% Cohen's D Bar Graphs
CD_in=[];
CD_out=[];
CD_in_CE_neg=[];
CD_in_CE_pos=[];
CD_out_CE_neg=[];
CD_out_CE_pos=[];
labels={'EC In','EC Out','DG In','DG Out','CA3 In','CA3 Out','CA1 In','CA1 Out'};

for i=1:4
    CD_in=[CD_in;[effect_in{i}{1,1}{1}.Effect,effect_in{i}{1,2}{1}.Effect,effect_in{i}{1,3}{1}.Effect]];
    CD_in_CE_neg=[CD_in_CE_neg;[effect_in{i}{1,1}{1}.ConfidenceIntervals(1),effect_in{i}{1,2}{1}.ConfidenceIntervals(1),effect_in{i}{1,3}{1}.ConfidenceIntervals(1)]];
    CD_in_CE_pos=[CD_in_CE_pos;[effect_in{i}{1,1}{1}.ConfidenceIntervals(2),effect_in{i}{1,2}{1}.ConfidenceIntervals(2),effect_in{i}{1,3}{1}.ConfidenceIntervals(2)]];
    CD_out=[CD_out;[effect_out{i}{1,1}{1}.Effect,effect_out{i}{1,2}{1}.Effect,effect_out{i}{1,3}{1}.Effect]];
    CD_out_CE_neg=[CD_out_CE_neg;[effect_out{i}{1,1}{1}.ConfidenceIntervals(1),effect_out{i}{1,2}{1}.ConfidenceIntervals(1),effect_out{i}{1,3}{1}.ConfidenceIntervals(1)]];
    CD_out_CE_pos=[CD_out_CE_pos;[effect_out{i}{1,1}{1}.ConfidenceIntervals(2),effect_out{i}{1,2}{1}.ConfidenceIntervals(2),effect_out{i}{1,3}{1}.ConfidenceIntervals(2)]];
end

CD_m=interleave_rows(CD_in,CD_out);
CD_CE_neg=interleave_rows(CD_in_CE_neg,CD_out_CE_neg);
CD_CE_pos=interleave_rows(CD_in_CE_pos,CD_out_CE_pos);

figure( 'Position', [100 100 700 600])
plot_horzbar_CohenD...
    (CD_m,CD_CE_neg,CD_CE_pos,labels);
xlabel("R^2 Directional Mean Effect Size")
xlim([-1.5,1.5])
ylim([0,9])
set(gca,'LineWidth',2,'TickLength',[0.02 0.02]);
set(gca,'FontSize',18)
exportgraphics(gcf,".\Temporal Analysis\Graph Analysis\rsq_cohenD.png",'Resolution',1500)
%% Average increase in RSQ

for i=1:4
    wt_ff_diff(i,1)=(diff(wt_ff_m(i,[1,2]))/wt_ff_m(i,1))*100;
    wt_ff_diff(i,2)=(diff(wt_ff_m(i,[1,3]))/wt_ff_m(i,1))*100;
    wt_fb_diff(i,1)=(diff(wt_fb_m(i,[1,2]))/wt_fb_m(i,1))*100;
    wt_fb_diff(i,2)=(diff(wt_fb_m(i,[1,3]))/wt_fb_m(i,1))*100;
    tw_ff_diff(i,1)=(diff(tw_ff_m(i,[1,2]))/tw_ff_m(i,1))*100;
    tw_ff_diff(i,2)=(diff(tw_ff_m(i,[1,3]))/tw_ff_m(i,1))*100;
    tw_fb_diff(i,1)=(diff(tw_fb_m(i,[1,2]))/tw_fb_m(i,1))*100;
    tw_fb_diff(i,2)=(diff(tw_fb_m(i,[1,3]))/tw_fb_m(i,1))*100;
end
pos_wtff=wt_ff_diff(wt_ff_diff>0);
pos_wtfb=wt_fb_diff(wt_fb_diff>0);
pos_twff=tw_ff_diff(tw_ff_diff>0);
pos_twfb=tw_fb_diff(tw_fb_diff>0);

average_rsq_increase=mean([pos_wtff;pos_wtfb;pos_twff;pos_twfb]);
average_rsq_increase_ff=mean([pos_wtff;pos_twff]);
average_rsq_increase_fb=mean([pos_wtfb;pos_twfb]);
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

function [G_cat_slopes,G_cat_rsq]=cat_edges(allregion_unit_matched,source,target,NodeTable,dataInfo)
% does not save tunnel-tunnel interactions in G for slopes or rsq
G_cat_slopes=[];
G_cat_rsq=[];

time_vec = [0.05:0.1:299.95];
no_ti = length(time_vec);
for graph_type = ["slope","rsq"]
    for fi=1:length(allregion_unit_matched)
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
    close all
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
%% Discontinued functions
% function [ff_source,ff_target,fb_source,fb_target] = split_dir_nodes(G,source,target, data_fldr, dataInfo)
% 
% ff_source=[];
% ff_target=[];
% fb_source=[];
% fb_target=[];
% 
% %% get source electrodes
% for fi=1:length(G)
%     sub_source=source.AFR_table(source.AFR_table.fi==fi,:);
%     source_wells=sub_source.well_names;
%     source_tunnels=sub_source.tunnel_names;
%     source_tunnels=split(source_tunnels,"-");
%     source_tunnels=source_tunnels(:,1);
%     
%     [source_electrodes,source_index]=unique([source_wells,source_tunnels],'rows','stable');
%     
%     source_reg=sub_source.well_reg(source_index);
% %     regions=unique(source_reg,'stable');
%     regions=["EC","DG","CA3","CA1"];
%     
%     for reg=1:length(unique(source_reg))
%         mea=load(data_fldr+'\mea_'+regions(reg)+'.mat');
%         reg_idx=(source_reg==regions(reg));
%         reg_elec=source_electrodes(reg_idx,:);
%         
%         % 0.6 is FB, 0.1 is FF
%         edges2use=and(ismember(cellstr(G{fi}.Edges.EndNodes(:,1)),...
%             reg_elec(:,1)),ismember(cellstr(G{fi}.Edges.EndNodes(:,2)),reg_elec(:,2)));
%         ff_g_edge=G{fi}.Edges(and(edges2use,G{fi}.Edges.Code==0.1),:);
%         fb_g_edge=G{fi}.Edges(and(edges2use,G{fi}.Edges.Code==0.6),:);
%         
%         ff_g_nodes=string(ff_g_edge.EndNodes);
%         ff_g_nodes=ff_g_nodes(:);
%         ff_g_nodes=G{fi}.Nodes(contains(G{fi}.Nodes.Name,ff_g_nodes),:);
% %         ff_g_wells=string(ff_g_edge.EndNodes);
% %         ff_g_wells=ff_g_wells(:,1);
%         
%         fb_g_nodes=string(fb_g_edge.EndNodes);
%         fb_g_nodes=fb_g_nodes(:);
%         fb_g_nodes=G{fi}.Nodes(contains(G{fi}.Nodes.Name,fb_g_nodes),:);
% %         fb_g_wells=string(fb_g_edge.EndNodes);
% %         fb_g_wells=fb_g_wells(:,1);
%         
%         if dataInfo.orientation=="CW"
%             wells=mea.cw.channel_names;
%         else
%             wells=mea.ccw.channel_names;
%         end
%         
%         ff_source{fi}{reg,1}=ff_g_edge;
%         ff_source{fi}{reg,2}=ff_g_nodes;
% %         ff_source{fi}{reg,3}=unique(ff_g_wells,'stable');
%         ff_source{fi}{reg,3}=wells;
%         ff_source{fi}{reg,4}=regions(reg);
%         if dataInfo.orientation=="CW"
%             ff_source{fi}{reg,5}="CW";
%         else
%             ff_source{fi}{reg,5}="CCW";
%         end
%         
%         fb_source{fi}{reg,1}=fb_g_edge;
%         fb_source{fi}{reg,2}=fb_g_nodes;
% %         fb_source{fi}{reg,3}=unique(fb_g_wells,'stable');
%         fb_source{fi}{reg,3}=wells;
%         fb_source{fi}{reg,4}=regions(reg);
%         if dataInfo.orientation=="CW"
%             fb_source{fi}{reg,5}="CW";
%         else
%             fb_source{fi}{reg,5}="CCW";
%         end
%         
%     end
%     ff_source{fi}=cell2table(ff_source{fi},"VariableNames",["Edges","Nodes","Electrodes","Subregion","Orientation"]);
%     fb_source{fi}=cell2table(fb_source{fi},"VariableNames",["Edges","Nodes","Electrodes","Subregion","Orientation"]);
% 
% end
% 
% %% Get Target Electrodes
% % for fi=1:length(G)
% %     sub_target=target.AFR_table(target.AFR_table.fi==fi,:);
% %     target_wells=sub_target.well_names;
% %     target_tunnels=sub_target.tunnel_names;
% %     target_tunnels=split(target_tunnels,"-");
% %     target_tunnels=target_tunnels(:,1);
% %     
% %     [target_electrodes,target_index]=unique([target_tunnels,target_wells],'rows','stable');
% %     
% %     target_reg=sub_target.tunnel_reg(target_index);
% % %     regions=unique(target_reg,'stable');
% %     regions=["EC-DG","DG-CA3","CA3-CA1","CA1-EC"];
% %     
% %     for reg=1:length(unique(target_reg))
% % %         mea=load(data_fldr+'\mea_'+regions(reg)+'.mat');
% %         reg_idx=(target_reg==regions(reg));
% %         reg_elec=target_electrodes(reg_idx,:);
% %         
% %         % 0.6 is FB, 0.1 is FF
% %         edges2use=and(ismember(cellstr(G{fi}.Edges.EndNodes(:,1)),...
% %             reg_elec(:,1)),ismember(cellstr(G{fi}.Edges.EndNodes(:,2)),reg_elec(:,2)));
% %         ff_g_edge=G{fi}.Edges(and(edges2use,G{fi}.Edges.Code==0.1),:);
% %         fb_g_edge=G{fi}.Edges(and(edges2use,G{fi}.Edges.Code==0.6),:);
% %         
% %         ff_g_nodes=string(ff_g_edge.EndNodes);
% %         ff_g_nodes=ff_g_nodes(:);
% %         ff_g_nodes=G{fi}.Nodes(contains(G{fi}.Nodes.Name,ff_g_nodes),:);
% % %         ff_g_wells=string(ff_g_edge.EndNodes);
% % %         ff_g_wells=ff_g_wells(:,1);
% %         
% %         fb_g_nodes=string(fb_g_edge.EndNodes);
% %         fb_g_nodes=fb_g_nodes(:);
% %         fb_g_nodes=G{fi}.Nodes(contains(G{fi}.Nodes.Name,fb_g_nodes),:);
% % %         fb_g_wells=string(fb_g_edge.EndNodes);
% % %         fb_g_wells=fb_g_wells(:,1);
% %         
% %         if dataInfo.orientation(fi)=="CW"
% %             tunnels=load(data_fldr+'\cw_vars_tunnel.mat');
% %             tunnel_elec=tunnels.channel_names(string(tunnels.subregions)==regions(reg));
% %         else
% %             tunnels=load(data_fldr+'\ccw_vars_tunnel.mat');
% %             tunnel_elec=tunnels.channel_names(string(tunnels.subregions)==regions(reg));
% %         end
% %         
% %         ff_target{fi}{reg,1}=ff_g_edge;
% %         ff_target{fi}{reg,2}=ff_g_nodes;
% % %         ff_target{fi}{reg,3}=unique(ff_g_wells,'stable');
% %         ff_target{fi}{reg,3}=tunnel_elec;
% %         ff_target{fi}{reg,4}=regions(reg);
% %         if dataInfo.orientation=="CW"
% %             ff_target{fi}{reg,5}="CW";
% %         else
% %             ff_target{fi}{reg,5}="CCW";
% %         end
% %         
% %         fb_target{fi}{reg,1}=fb_g_edge;
% %         fb_target{fi}{reg,2}=fb_g_nodes;
% % %         fb_target{fi}{reg,3}=unique(fb_g_wells,'stable');
% %         fb_target{fi}{reg,3}=tunnel_elec;
% %         fb_target{fi}{reg,4}=regions(reg);
% %         if dataInfo.orientation=="CW"
% %             fb_target{fi}{reg,5}="CW";
% %         else
% %             fb_target{fi}{reg,5}="CCW";
% %         end
% %         
% %     end
% %     ff_target{fi}=cell2table(ff_target{fi},"VariableNames",["Edges","Nodes","Electrodes","Subregion","Orientation"]);
% %     fb_target{fi}=cell2table(fb_target{fi},"VariableNames",["Edges","Nodes","Electrodes","Subregion","Orientation"]);
% % 
% % end
% end
% 
% function [source_out,target_out]=regional_out_nodes(G,data_fldr,dataInfo)
% 
% source_out=[];
% target_out=[];
% 
% %source out
% regions=["EC","DG","CA3","CA1"];
% for fi=1:length(G)
%     for regi=1:4
%         mea=load(data_fldr+'\mea_'+regions(reg)+'.mat');
%         if dataInfo{fi}.orientation=="CW"
%             
%         else
%             
%         end
%     end
% end
% 
% end