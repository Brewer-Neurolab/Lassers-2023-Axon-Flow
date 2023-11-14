%% IBI comparison bar graphs and anova
clear 
clc
ff_colors=["#3869EE","#FF2303","#EE9B00"];
fb_colors=["#4fb7ff","#FF9B9B","#FFDA49"];
%% Feedforward lower

nostim=load('IBI_no_stim_ff_ancova_vars_1.mat');
theta=load('IBI_theta_stim_ff_ancova_vars_1.mat');
hfs=load('IBI_HFS_ff_ancova_vars_1.mat');

subregions_ff=["EC-DG","DG-CA3","CA3-CA1","CA1-EC"];
ttest_struct=[];
perc_change=[];
for i=1:length(subregions_ff)
    nostim_idx=contains(string(nostim.regLabel_post_ff_1),subregions_ff(i));
    theta_idx=contains(string(theta.regLabel_post_ff_1),subregions_ff(i));
    hfs_idx=contains(string(hfs.regLabel_post_ff_1),subregions_ff(i));
    sp_feature=[nostim.sp_feature_post_ff_1(nostim_idx);theta.sp_feature_post_ff_1(theta_idx);hfs.sp_feature_post_ff_1(hfs_idx)];
    prob_vec=[nostim.prob_vec_post_ff_1(nostim_idx);theta.prob_vec_post_ff_1(theta_idx);hfs.prob_vec_post_ff_1(hfs_idx)];
    regLabel=[nostim.regLabel_post_ff_1(nostim_idx);theta.regLabel_post_ff_1(theta_idx);hfs.regLabel_post_ff_1(hfs_idx)];
    
    [~,~,~,stats] = aoctool(sp_feature,prob_vec,regLabel,0.05,'','','','off');
    [ff_c{i},ff_means{i}]=multcompare(stats,0.05,'off','','s');
    
    perc_change{i}{1}(1,:)=[1,2];
    perc_change{i}{1}(2,:)=[1,3];
    perc_change{i}{1}(3,:)=[2,3];
    
    perc_change{i}{2}(1,:)=-(ff_means{i}(2,1)-ff_means{i}(1,1))/(ff_means{i}(1,1))*100;
    perc_change{i}{2}(2,:)=-(ff_means{i}(3,1)-ff_means{i}(1,1))/(ff_means{i}(1,1))*100;
    perc_change{i}{2}(3,:)=-(ff_means{i}(3,1)-ff_means{i}(2,1))/(ff_means{i}(2,1))*100;
    
    ttest_struct{i}.groups(1,:)=[1,2];
    ttest_struct{i}.groups(2,:)=[1,3];
    ttest_struct{i}.groups(3,:)=[2,3];
    [~,pval]=ttest2(nostim.prob_vec_post_ff_1(nostim_idx),theta.prob_vec_post_ff_1(theta_idx));
    ttest_struct{i}.pval(1)=pval;
    [~,pval]=ttest2(nostim.prob_vec_post_ff_1(nostim_idx),hfs.prob_vec_post_ff_1(hfs_idx));
    ttest_struct{i}.pval(2)=pval;
    [~,pval]=ttest2(theta.prob_vec_post_ff_1(theta_idx),hfs.prob_vec_post_ff_1(hfs_idx));
    ttest_struct{i}.pval(3)=pval;
end

% bar graphs
ydata=[];
errordata=[];
for i=1:length(subregions_ff)
    ydata=[ydata;ff_means{i}(:,1)'];
    errordata=[errordata;ff_means{i}(:,2)'];
end

figure( 'Position', [100 100 700 600])
subreg_lab=categorical(subregions_ff);
subreg_lab=reordercats(subreg_lab,subregions_ff);
b=bar(subreg_lab,ydata);
hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(ydata);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
x=x';
errorbar(x,ydata,errordata,'k','linestyle','none');
ax=gca;
ax.FontSize=24;
%legend('No Stimulation','Theta Stimulation','HFS Stimulation')
% title("Feedforward IBI")
ylabel 'Slope of IBI (s^{-1})'
ylim([-1,0])
hold off
ax.TickLength=[.02,.02];
ax.LineWidth=2;
exportgraphics(gcf,'.\three stim share figs\ff_ibi_bar_1.png','Resolution',1500)


for i=1:length(subregions_ff)
    ax.XAxis.Limits=categorical({char(subregions_ff(i)),char(subregions_ff(i))});
    exportgraphics(gcf,strcat('.\three stim share figs\ff_ibi_bar_1',subregions_ff(i),'.png'),'Resolution',1500)
end

%% Feedback lower

nostim=load('IBI_no_stim_fb_ancova_vars_1.mat');
theta=load('IBI_theta_stim_fb_ancova_vars_1.mat');
hfs=load('IBI_HFS_fb_ancova_vars_1.mat');

subregions_fb=["EC-DG","DG-CA3","CA3-CA1","CA1-EC"];
ttest_struct=[];
perc_change=[];
for i=1:length(subregions_fb)
    nostim_idx=contains(string(nostim.regLabel_post_fb_1),subregions_fb(i));
    theta_idx=contains(string(theta.regLabel_post_fb_1),subregions_fb(i));
    hfs_idx=contains(string(hfs.regLabel_post_fb_1),subregions_fb(i));
    sp_feature=[nostim.sp_feature_post_fb_1(nostim_idx);theta.sp_feature_post_fb_1(theta_idx);hfs.sp_feature_post_fb_1(hfs_idx)];
    prob_vec=[nostim.prob_vec_post_fb_1(nostim_idx);theta.prob_vec_post_fb_1(theta_idx);hfs.prob_vec_post_fb_1(hfs_idx)];
    regLabel=[nostim.regLabel_post_fb_1(nostim_idx);theta.regLabel_post_fb_1(theta_idx);hfs.regLabel_post_fb_1(hfs_idx)];
    
    [~,~,~,stats] = aoctool(sp_feature,prob_vec,regLabel,0.05,'','','','off');
    [fb_c{i},fb_means{i}]=multcompare(stats,0.05,'off','','s');
    
    perc_change{i}{1}(1,:)=[1,2];
    perc_change{i}{1}(2,:)=[1,3];
    perc_change{i}{1}(3,:)=[2,3];
    
    perc_change{i}{2}(1,:)=-(fb_means{i}(2,1)-fb_means{i}(1,1))/(fb_means{i}(1,1))*100;
    perc_change{i}{2}(2,:)=-(fb_means{i}(3,1)-fb_means{i}(1,1))/(fb_means{i}(1,1))*100;
    perc_change{i}{2}(3,:)=-(fb_means{i}(3,1)-fb_means{i}(2,1))/(fb_means{i}(2,1))*100;
    
    ttest_struct{i}.groups(1,:)=[1,2];
    ttest_struct{i}.groups(2,:)=[1,3];
    ttest_struct{i}.groups(3,:)=[2,3];
    [~,pval]=ttest2(nostim.prob_vec_post_fb_1(nostim_idx),theta.prob_vec_post_fb_1(theta_idx));
    ttest_struct{i}.pval(1)=pval;
    [~,pval]=ttest2(nostim.prob_vec_post_fb_1(nostim_idx),hfs.prob_vec_post_fb_1(hfs_idx));
    ttest_struct{i}.pval(2)=pval;
    [~,pval]=ttest2(theta.prob_vec_post_fb_1(theta_idx),hfs.prob_vec_post_fb_1(hfs_idx));
    ttest_struct{i}.pval(3)=pval;
end

% bar graphs
ydata=[];
errordata=[];
for i=1:length(subregions_fb)
    ydata=[ydata;fb_means{i}(:,1)'];
    errordata=[errordata;fb_means{i}(:,2)'];
end

figure( 'Position', [100 100 700 600])
subregions_fb=["DG-EC","CA3-DG","CA1-CA3","EC-CA1"];
subreg_lab=categorical(subregions_fb);
subreg_lab=reordercats(subreg_lab,subregions_fb);
b=bar(subreg_lab,ydata);
hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(ydata);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
x=x';
errorbar(x,ydata,errordata,'k','linestyle','none');
ax=gca;
ax.FontSize=24;
%legend('No Stimulation','Theta Stimulation','HFS Stimulation')
% title("Feedback IBI")
ylabel 'Slope of IBI (s^{-1})'
ylim([-1,0])
hold off
ax.TickLength=[.02,.02];
ax.LineWidth=2;
exportgraphics(gcf,'.\three stim share figs\fb_ibi_bar_1.png','Resolution',1500)


for i=1:length(subregions_fb)
    ax.XAxis.Limits=categorical({char(subregions_fb(i)),char(subregions_fb(i))});
    exportgraphics(gcf,strcat('.\three stim share figs\fb_ibi_bar_1',subregions_fb(i),'.png'),'Resolution',1500)
end
%% Horzbar Slopes lower

dat2plot_ff=[];
dat2plot_ff_se=[];
dat2plot_fb=[];
dat2plot_fb_se=[];
labels={'EC-DG','DG-CA3','CA3-CA1','CA1-EC'};
figure( 'Position', [100 100 1400 600])
for i=1:4
    dat2plot_ff=[dat2plot_ff,ff_means{i}(:,1)];
    dat2plot_ff_se=[dat2plot_ff_se,ff_means{i}(:,2)];
    dat2plot_fb=[dat2plot_fb,fb_means{i}(:,1)];
    dat2plot_fb_se=[dat2plot_fb_se,fb_means{i}(:,2)];
end

plot_horzbar_from_anova...
    (-dat2plot_ff',dat2plot_ff_se',-dat2plot_fb',dat2plot_fb_se',labels,ff_colors,fb_colors);
xlabel("FB Slope of IBI (s^{-1})                                      FF Slope of IBI (s^{-1})")
xlim([-1,1])
xticks([-1.5:.1:1.5])
% ylim([0,5])
set(gca,'LineWidth',2,'TickLength',[0.02 0.02]);
set(gca,'FontSize',12)

exportgraphics(gcf,strcat('.\three stim share figs\ff_fb_ibi_1_horzbar.png'),'Resolution',1500)
%%
labels=flip({'EC-DG','DG-CA3','CA3-CA1','CA1-EC'});
for i=4:-1:1
%     dat2plot_ff=[];
%     dat2plot_ff_se=[];
%     dat2plot_fb=[];
%     dat2plot_fb_se=[];
% 
%     figure( 'Position', [100 100 1400 600])
%     dat2plot_ff=ff_means{i}(:,1);
%     dat2plot_ff_se=ff_means{i}(:,2);
%     dat2plot_fb=fb_means{i}(:,1);
%     dat2plot_fb_se=fb_means{i}(:,2);
% 
%     plot_horzbar_from_anova...
%         (-dat2plot_ff',dat2plot_ff_se',-dat2plot_fb',dat2plot_fb_se',labels{i},ff_colors(i),fb_colors(i));
%     xlabel("FB Slope                                      FF Slope")
%     xlim([-1.2,1.2])
%     xticks([-1.5:.1:1.5])
%     ylim([0,5])
%     set(gca,'LineWidth',2,'TickLength',[0.02 0.02]);
%     set(gca,'FontSize',12)
    ylim([i-.5,i+.5])
    exportgraphics(gcf,strcat('.\three stim share figs\ff_fb_ibi_1_horzbar',labels{i},'.png'),'Resolution',1500)
end
%% FF/FB% 

ff_fb_perc=((dat2plot_ff-dat2plot_fb)./(dat2plot_ff+dat2plot_fb))*100;

ff_lower=dat2plot_ff;
ff_lower_se=dat2plot_ff_se;
fb_lower=dat2plot_fb;
fb_lower_se=dat2plot_fb_se;

% ff_fb_SD=dat2plot

ff_fb_SD=[];
ff_fb_lengths=[];

labels={'EC-DG','DG-CA3','CA3-CA1','CA1-EC'};
regs=categorical(labels);
regs=reordercats(regs,labels);

figure( 'Position', [100 100 700 600])
% b=barh(regs,ff_fb_perc);

b=plot_horzbar_CohenD(ff_fb_perc',zeros(size(ff_fb_perc')),zeros(size(ff_fb_perc')),labels);
b(1).BaseValue=0;
b(1).BaseLine.LineStyle = "--";
b(1).BaseLine.LineWidth=2;
xlim([-70,70])
xticks([-70:10:70])
xticklabels([-70:10:70])

xlabel("<- Excess FB            Excess FF ->")
set(gca,'LineWidth',2,'TickLength',[0.02 0.02]);
set(gca,'FontSize',18)
exportgraphics(gcf,'.\three stim share figs\ff_fb_IBI_low_perc.png','Resolution',1500)
%% Well Lower

nostim=load('IBI_no_stim_well_ancova_vars_1_5SD_500maxSD.mat');
theta=load('IBI_Theta_well_ancova_vars_1_5SD_500maxSD.mat');
hfs=load('IBI_HFS_well_ancova_vars_1_5SD_500maxSD.mat');

subregions_well=["EC","DG","CA3","CA1"];
ttest_struct=[];
perc_change=[];
for i=1:length(subregions_well)
    nostim_idx=contains(string(nostim.regLabel_post_well_1),subregions_well(i));
    theta_idx=contains(string(theta.regLabel_post_well_1),subregions_well(i));
    hfs_idx=contains(string(hfs.regLabel_post_well_1),subregions_well(i));
    sp_feature=[nostim.sp_feature_post_well_1(nostim_idx);theta.sp_feature_post_well_1(theta_idx);hfs.sp_feature_post_well_1(hfs_idx)];
    prob_vec=[nostim.prob_vec_post_well_1(nostim_idx);theta.prob_vec_post_well_1(theta_idx);hfs.prob_vec_post_well_1(hfs_idx)];
    regLabel=[nostim.regLabel_post_well_1(nostim_idx);theta.regLabel_post_well_1(theta_idx);hfs.regLabel_post_well_1(hfs_idx)];
    
    [~,~,~,stats] = aoctool(sp_feature,prob_vec,regLabel,0.05,'','','','off');
    [ff_c{i},ff_means{i}]=multcompare(stats,0.05,'off','','s');
    
    perc_change{i}{1}(1,:)=[1,2];
    perc_change{i}{1}(2,:)=[1,3];
    perc_change{i}{1}(3,:)=[2,3];
    
    perc_change{i}{2}(1,:)=-(ff_means{i}(2,1)-ff_means{i}(1,1))/(ff_means{i}(1,1))*100;
    perc_change{i}{2}(2,:)=-(ff_means{i}(3,1)-ff_means{i}(1,1))/(ff_means{i}(1,1))*100;
    perc_change{i}{2}(3,:)=-(ff_means{i}(3,1)-ff_means{i}(2,1))/(ff_means{i}(2,1))*100;
    
    ttest_struct{i}.groups(1,:)=[1,2];
    ttest_struct{i}.groups(2,:)=[1,3];
    ttest_struct{i}.groups(3,:)=[2,3];
    [~,pval]=ttest2(nostim.prob_vec_post_well_1(nostim_idx),theta.prob_vec_post_well_1(theta_idx));
    ttest_struct{i}.pval(1)=pval;
    [~,pval]=ttest2(nostim.prob_vec_post_well_1(nostim_idx),hfs.prob_vec_post_well_1(hfs_idx));
    ttest_struct{i}.pval(2)=pval;
    [~,pval]=ttest2(theta.prob_vec_post_well_1(theta_idx),hfs.prob_vec_post_well_1(hfs_idx));
    ttest_struct{i}.pval(3)=pval;
end

% bar graphs
ydata=[];
errordata=[];
for i=1:length(subregions_well)
    ydata=[ydata;ff_means{i}(:,1)'];
    errordata=[errordata;ff_means{i}(:,2)'];
end

figure( 'Position', [100 100 700 600])
subreg_lab=categorical(subregions_well);
subreg_lab=reordercats(subreg_lab,subregions_well);
b=bar(subreg_lab,ydata);

well_lower=ydata';

hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(ydata);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
x=x';
errorbar(x,ydata,errordata,'k','linestyle','none');
ax=gca;
ax.FontSize=24;
%legend('No Stimulation','Theta Stimulation','HFS Stimulation')
% title("Lower Well IBI")
ylabel 'Slope of IBI (s^{-1})'
ylim([-1,0])
hold off
ax.TickLength=[.02,.02];
ax.LineWidth=2;
exportgraphics(gcf,'.\three stim share figs\well_5SD_500maxSD_ibi_bar_1.png','Resolution',1500)


for i=1:length(subregions_well)
    ax.XAxis.Limits=categorical({char(subregions_well(i)),char(subregions_well(i))});
    exportgraphics(gcf,strcat('.\three stim share figs\well_5SD_500maxSD_ibi_bar_1',subregions_well(i),'.png'),'Resolution',1500)
end
%% Feedforward upper

nostim=load('IBI_no_stim_ff_ancova_vars_2.mat');
theta=load('IBI_theta_stim_ff_ancova_vars_2.mat');
hfs=load('IBI_HFS_ff_ancova_vars_2.mat');

subregions_ff=["EC-DG","DG-CA3","CA3-CA1","CA1-EC"];
ttest_struct=[];
perc_change=[];
for i=1:length(subregions_ff)
    nostim_idx=contains(string(nostim.regLabel_post_ff_2),subregions_ff(i));
    theta_idx=contains(string(theta.regLabel_post_ff_2),subregions_ff(i));
    hfs_idx=contains(string(hfs.regLabel_post_ff_2),subregions_ff(i));
    sp_feature=[nostim.sp_feature_post_ff_2(nostim_idx);theta.sp_feature_post_ff_2(theta_idx);hfs.sp_feature_post_ff_2(hfs_idx)];
    prob_vec=[nostim.prob_vec_post_ff_2(nostim_idx);theta.prob_vec_post_ff_2(theta_idx);hfs.prob_vec_post_ff_2(hfs_idx)];
    regLabel=[nostim.regLabel_post_ff_2(nostim_idx);theta.regLabel_post_ff_2(theta_idx);hfs.regLabel_post_ff_2(hfs_idx)];
    
    [~,~,~,stats] = aoctool(sp_feature,prob_vec,regLabel,0.05,'','','','off');
    [ff_c{i},ff_means{i}]=multcompare(stats,0.05,'off','','s');
    
    perc_change{i}{1}(1,:)=[1,2];
    perc_change{i}{1}(2,:)=[1,3];
    perc_change{i}{1}(3,:)=[2,3];
    
    perc_change{i}{2}(1,:)=-(ff_means{i}(2,1)-ff_means{i}(1,1))/(ff_means{i}(1,1))*100;
    perc_change{i}{2}(2,:)=-(ff_means{i}(3,1)-ff_means{i}(1,1))/(ff_means{i}(1,1))*100;
    perc_change{i}{2}(3,:)=-(ff_means{i}(3,1)-ff_means{i}(2,1))/(ff_means{i}(2,1))*100;
    
    ttest_struct{i}.groups(1,:)=[1,2];
    ttest_struct{i}.groups(2,:)=[1,3];
    ttest_struct{i}.groups(3,:)=[2,3];
    [~,pval]=ttest2(nostim.prob_vec_post_ff_2(nostim_idx),theta.prob_vec_post_ff_2(theta_idx));
    ttest_struct{i}.pval(1)=pval;
    [~,pval]=ttest2(nostim.prob_vec_post_ff_2(nostim_idx),hfs.prob_vec_post_ff_2(hfs_idx));
    ttest_struct{i}.pval(2)=pval;
    [~,pval]=ttest2(theta.prob_vec_post_ff_2(theta_idx),hfs.prob_vec_post_ff_2(hfs_idx));
    ttest_struct{i}.pval(3)=pval;
end

% bar graphs
ydata=[];
errordata=[];
for i=1:length(subregions_ff)
    ydata=[ydata;ff_means{i}(:,1)'];
    errordata=[errordata;ff_means{i}(:,2)'];
end

figure( 'Position', [100 100 700 600])
subreg_lab=categorical(subregions_ff);
subreg_lab=reordercats(subreg_lab,subregions_ff);
b=bar(subreg_lab,ydata);

well_upper=ydata';

hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(ydata);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
x=x';
errorbar(x,ydata,errordata,'k','linestyle','none');
ax=gca;
ax.FontSize=24;
%legend('No Stimulation','Theta Stimulation','HFS Stimulation')
% title("Feedforward IBI")
ylabel 'Slope of IBI (s^{-1})'
ylim([-2.5,0])
hold off
ax.TickLength=[.02,.02];
ax.LineWidth=2;
exportgraphics(gcf,'.\three stim share figs\ff_ibi_bar_2.png','Resolution',1500)


for i=1:length(subregions_ff)
    ax.XAxis.Limits=categorical({char(subregions_ff(i)),char(subregions_ff(i))});
    exportgraphics(gcf,strcat('.\three stim share figs\ff_ibi_bar_2',subregions_ff(i),'.png'),'Resolution',1500)
end

%% Feedback upper

nostim=load('IBI_no_stim_fb_ancova_vars_2.mat');
theta=load('IBI_theta_stim_fb_ancova_vars_2.mat');
hfs=load('IBI_HFS_fb_ancova_vars_2.mat');

subregions_fb=["EC-DG","DG-CA3","CA3-CA1","CA1-EC"];
ttest_struct=[];
perc_change=[];
fb_means=[];
fb_c=[];
for i=1:length(subregions_fb)
    nostim_idx=contains(string(nostim.regLabel_post_fb_2),subregions_fb(i));
    theta_idx=contains(string(theta.regLabel_post_fb_2),subregions_fb(i));
    hfs_idx=contains(string(hfs.regLabel_post_fb_2),subregions_fb(i));
    sp_feature=[nostim.sp_feature_post_fb_2(nostim_idx);theta.sp_feature_post_fb_2(theta_idx);hfs.sp_feature_post_fb_2(hfs_idx)];
    prob_vec=[nostim.prob_vec_post_fb_2(nostim_idx);theta.prob_vec_post_fb_2(theta_idx);hfs.prob_vec_post_fb_2(hfs_idx)];
    regLabel=[nostim.regLabel_post_fb_2(nostim_idx);theta.regLabel_post_fb_2(theta_idx);hfs.regLabel_post_fb_2(hfs_idx)];
    
    %test
%     sp_feature=[theta.sp_feature_post_fb_2(theta_idx)];
%     prob_vec=[theta.prob_vec_post_fb_2(theta_idx)];
%     regLabel=[theta.regLabel_post_fb_2(theta_idx)];
    
    [~,~,~,stats] = aoctool(sp_feature,prob_vec,regLabel,0.05,'','','','off');
    [fb_c{i},fb_means{i}]=multcompare(stats,0.05,'off','','s');
    
    perc_change{i}{1}(1,:)=[1,2];
    perc_change{i}{1}(2,:)=[1,3];
    perc_change{i}{1}(3,:)=[2,3];
    
    perc_change{i}{2}(1,:)=-(fb_means{i}(2,1)-fb_means{i}(1,1))/(fb_means{i}(1,1))*100;
    perc_change{i}{2}(2,:)=-(fb_means{i}(3,1)-fb_means{i}(1,1))/(fb_means{i}(1,1))*100;
    perc_change{i}{2}(3,:)=-(fb_means{i}(3,1)-fb_means{i}(2,1))/(fb_means{i}(2,1))*100;
    
    ttest_struct{i}.groups(1,:)=[1,2];
    ttest_struct{i}.groups(2,:)=[1,3];
    ttest_struct{i}.groups(3,:)=[2,3];
    [~,pval]=ttest2(nostim.prob_vec_post_fb_2(nostim_idx),theta.prob_vec_post_fb_2(theta_idx));
    ttest_struct{i}.pval(1)=pval;
    [~,pval]=ttest2(nostim.prob_vec_post_fb_2(nostim_idx),hfs.prob_vec_post_fb_2(hfs_idx));
    ttest_struct{i}.pval(2)=pval;
    [~,pval]=ttest2(theta.prob_vec_post_fb_2(theta_idx),hfs.prob_vec_post_fb_2(hfs_idx));
    ttest_struct{i}.pval(3)=pval;
end

% bar graphs
ydata=[];
errordata=[];
for i=1:length(subregions_fb)
    ydata=[ydata;fb_means{i}(:,1)'];
    errordata=[errordata;fb_means{i}(:,2)'];
end

figure( 'Position', [100 100 700 600])
subregions_fb=["DG-EC","CA3-DG","CA1-CA3","EC-CA1"];
subreg_lab=categorical(subregions_fb);
subreg_lab=reordercats(subreg_lab,subregions_fb);
b=bar(subreg_lab,ydata);
hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(ydata);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
x=x';
errorbar(x,ydata,errordata,'k','linestyle','none');
ax=gca;
ax.FontSize=24;
%legend('No Stimulation','Theta Stimulation','HFS Stimulation')
% title("Feedback IBI")
ylabel 'Slope of IBI (s^{-1})'
ylim([-2.5,0])
hold off
ax.TickLength=[.02,.02];
ax.LineWidth=2;
exportgraphics(gcf,'.\three stim share figs\fb_ibi_bar_2.png','Resolution',1500)


for i=1:length(subregions_fb)
    ax.XAxis.Limits=categorical({char(subregions_fb(i)),char(subregions_fb(i))});
    exportgraphics(gcf,strcat('.\three stim share figs\fb_ibi_bar_2',subregions_fb(i),'.png'),'Resolution',1500)
end
%% Horzbar Slopes upper

dat2plot_ff=[];
dat2plot_ff_se=[];
dat2plot_fb=[];
dat2plot_fb_se=[];
labels={'EC-DG','DG-CA3','CA3-CA1','CA1-EC'};
figure( 'Position', [100 100 1400 600])
for i=1:4
    dat2plot_ff=[dat2plot_ff,ff_means{i}(:,1)];
    dat2plot_ff_se=[dat2plot_ff_se,ff_means{i}(:,2)];
    dat2plot_fb=[dat2plot_fb,fb_means{i}(:,1)];
    dat2plot_fb_se=[dat2plot_fb_se,fb_means{i}(:,2)];
end

plot_horzbar_from_anova...
    (-dat2plot_ff',dat2plot_ff_se',-dat2plot_fb',dat2plot_fb_se',labels,ff_colors,fb_colors);
xlabel("FB Slope of IBI (s^{-1})                                      FF Slope of IBI (s^{-1})")
xlim([-2.5,2.5])
xticks([-2.5:.1:2.5])
% ylim([0,5])
set(gca,'LineWidth',2,'TickLength',[0.02 0.02]);
set(gca,'FontSize',12)

exportgraphics(gcf,strcat('.\three stim share figs\ff_fb_ibi_2_horzbar.png'),'Resolution',1500)
%%
labels=flip({'EC-DG','DG-CA3','CA3-CA1','CA1-EC'});
for i=4:-1:1
%     dat2plot_ff=[];
%     dat2plot_ff_se=[];
%     dat2plot_fb=[];
%     dat2plot_fb_se=[];
% 
%     figure( 'Position', [100 100 1400 600])
%     dat2plot_ff=ff_means{i}(:,1);
%     dat2plot_ff_se=ff_means{i}(:,2);
%     dat2plot_fb=fb_means{i}(:,1);
%     dat2plot_fb_se=fb_means{i}(:,2);
% 
%     plot_horzbar_from_anova...
%         (-dat2plot_ff',dat2plot_ff_se',-dat2plot_fb',dat2plot_fb_se',labels{i},ff_colors(i),fb_colors(i));
%     xlabel("FB Slope                                      FF Slope")
%     xlim([-1.2,1.2])
%     xticks([-1.5:.1:1.5])
%     ylim([0,5])
%     set(gca,'LineWidth',2,'TickLength',[0.02 0.02]);
%     set(gca,'FontSize',12)
    ylim([i-.5,i+.5])
    exportgraphics(gcf,strcat('.\three stim share figs\ff_fb_ibi_2_horzbar',labels{i},'.png'),'Resolution',1500)
end
%% FF/FB% 

ff_fb_perc=((dat2plot_ff-dat2plot_fb)./(dat2plot_ff+dat2plot_fb))*100;

% ff_fb_SD=dat2plot

ff_upper=dat2plot_ff;
ff_upper_se=dat2plot_ff_se;
fb_upper=dat2plot_fb;
fb_upper_se=dat2plot_fb_se;

ff_fb_SD=[];
ff_fb_lengths=[];

labels={'EC-DG','DG-CA3','CA3-CA1','CA1-EC'};
regs=categorical(labels);
regs=reordercats(regs,labels);

figure( 'Position', [100 100 700 600])
% b=barh(regs,ff_fb_perc);

b=plot_horzbar_CohenD(ff_fb_perc',zeros(size(ff_fb_perc')),zeros(size(ff_fb_perc')),labels);
b(1).BaseValue=0;
b(1).BaseLine.LineStyle = "--";
b(1).BaseLine.LineWidth=2;
xlim([-50,50])
xticks([-50:10:50])
xticklabels([-50:10:50])

xlabel("<- Excess FB            Excess FF ->")
set(gca,'LineWidth',2,'TickLength',[0.02 0.02]);
set(gca,'FontSize',18)
exportgraphics(gcf,'.\three stim share figs\ff_fb_IBI_high_perc.png','Resolution',1500)
%% Well Upper

nostim=load('IBI_no_stim_well_ancova_vars_2_5SD_500maxSD.mat');
theta=load('IBI_theta_well_ancova_vars_2_5SD_500maxSD.mat');
hfs=load('IBI_HFS_well_ancova_vars_2_5SD_500maxSD.mat');

subregions_well=["EC","DG","CA3","CA1"];
ttest_struct=[];
perc_change=[];
for i=1:length(subregions_well)
    nostim_idx=contains(string(nostim.regLabel_post_well_2),subregions_well(i));
    theta_idx=contains(string(theta.regLabel_post_well_2),subregions_well(i));
    hfs_idx=contains(string(hfs.regLabel_post_well_2),subregions_well(i));
    sp_feature=[nostim.sp_feature_post_well_2(nostim_idx);theta.sp_feature_post_well_2(theta_idx);hfs.sp_feature_post_well_2(hfs_idx)];
    prob_vec=[nostim.prob_vec_post_well_2(nostim_idx);theta.prob_vec_post_well_2(theta_idx);hfs.prob_vec_post_well_2(hfs_idx)];
    regLabel=[nostim.regLabel_post_well_2(nostim_idx);theta.regLabel_post_well_2(theta_idx);hfs.regLabel_post_well_2(hfs_idx)];
    
    [~,~,~,stats] = aoctool(sp_feature,prob_vec,regLabel,0.05,'','','','off');
    [ff_c{i},ff_means{i}]=multcompare(stats,0.05,'off','','s');
    
    perc_change{i}{1}(1,:)=[1,2];
    perc_change{i}{1}(2,:)=[1,3];
    perc_change{i}{1}(3,:)=[2,3];
    
    perc_change{i}{2}(1,:)=-(ff_means{i}(2,1)-ff_means{i}(1,1))/(ff_means{i}(1,1))*100;
    perc_change{i}{2}(2,:)=-(ff_means{i}(3,1)-ff_means{i}(1,1))/(ff_means{i}(1,1))*100;
    perc_change{i}{2}(3,:)=-(ff_means{i}(3,1)-ff_means{i}(2,1))/(ff_means{i}(2,1))*100;
    
    ttest_struct{i}.groups(1,:)=[1,2];
    ttest_struct{i}.groups(2,:)=[1,3];
    ttest_struct{i}.groups(3,:)=[2,3];
    [~,pval]=ttest2(nostim.prob_vec_post_well_2(nostim_idx),theta.prob_vec_post_well_2(theta_idx));
    ttest_struct{i}.pval(1)=pval;
    [~,pval]=ttest2(nostim.prob_vec_post_well_2(nostim_idx),hfs.prob_vec_post_well_2(hfs_idx));
    ttest_struct{i}.pval(2)=pval;
    [~,pval]=ttest2(theta.prob_vec_post_well_2(theta_idx),hfs.prob_vec_post_well_2(hfs_idx));
    ttest_struct{i}.pval(3)=pval;
end

% bar graphs
ydata=[];
errordata=[];
for i=1:length(subregions_well)
    ydata=[ydata;ff_means{i}(:,1)'];
    errordata=[errordata;ff_means{i}(:,2)'];
end

figure( 'Position', [100 100 700 600])
subreg_lab=categorical(subregions_well);
subreg_lab=reordercats(subreg_lab,subregions_well);
b=bar(subreg_lab,ydata);
hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(ydata);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
x=x';
errorbar(x,ydata,errordata,'k','linestyle','none');
ax=gca;
ax.FontSize=24;
%legend('No Stimulation','Theta Stimulation','HFS Stimulation')
% title("Higher Well IBI")
ylabel 'Slope of IBI (s^{-1})'
ylim([-2.5,0])
hold off
ax.TickLength=[.02,.02];
ax.LineWidth=2;
exportgraphics(gcf,'.\three stim share figs\well_5SD_500maxSD_ibi_bar_2.png','Resolution',1500)


for i=1:length(subregions_well)
    ax.XAxis.Limits=categorical({char(subregions_well(i)),char(subregions_well(i))});
    exportgraphics(gcf,strcat('.\three stim share figs\well_5SD_500maxSD_ibi_bar_2',subregions_well(i),'.png'),'Resolution',1500)
end
%% Bar graphs by subregion ff 

nostim_plot=load('ff_IBI_points_plotted_nostim.mat');
nostim_fit_1=load('ff_IBI_fit_plotted_nostim_1.mat');
nostim_fit_2=load('ff_IBI_fit_plotted_nostim_2.mat');
theta_plot=load('ff_IBI_points_plotted_thetastim.mat');
theta_fit_1=load('ff_IBI_fit_plotted_thetastim_1.mat');
theta_fit_2=load('ff_IBI_fit_plotted_thetastim_2.mat');
hfs_plot=load('ff_IBI_points_plotted_HFS.mat');
hfs_fit_1=load('ff_IBI_fit_plotted_HFS_1.mat');
hfs_fit_2=load('ff_IBI_fit_plotted_HFS_2.mat');

subregions_ff=["EC-DG","DG-CA3","CA3-CA1","CA1-EC"];
subregions_fb=["DG-EC","CA3-DG","CA1-CA3","EC-CA1"];

prob_fifty_ISI_1=[];
prob_fifty_ISI_2=[];
for i=1:length(subregions_ff)
    figure( 'Position', [100 100 700 600])
    hold on
    plot(nostim_plot.points_plotted{i}{1},nostim_plot.points_plotted{i}{2},'Color',ff_colors(1),'LineWidth',4)
    plot(nostim_fit_1.plotted_fit_1{i}{1},nostim_fit_1.plotted_fit_1{i}{2},'--k','LineWidth',4)
    plot(nostim_fit_2.plotted_fit_2{i}{1},nostim_fit_2.plotted_fit_2{i}{2},':','LineWidth',4,'Color',[0.5 0.5 0.5])
    plot(theta_plot.points_plotted{i}{1},theta_plot.points_plotted{i}{2},'Color',ff_colors(2),'LineWidth',4)
    plot(theta_fit_1.plotted_fit_1{i}{1},theta_fit_1.plotted_fit_1{i}{2},'--k','LineWidth',4)
    plot(theta_fit_2.plotted_fit_2{i}{1},theta_fit_2.plotted_fit_2{i}{2},':','LineWidth',4,'Color',[0.5 0.5 0.5])
    plot(hfs_plot.points_plotted{i}{1},hfs_plot.points_plotted{i}{2},'Color',ff_colors(3),'LineWidth',4)
    plot(hfs_fit_1.plotted_fit_1{i}{1},hfs_fit_1.plotted_fit_1{i}{2},'--k','LineWidth',4)
    plot(hfs_fit_2.plotted_fit_2{i}{1},hfs_fit_2.plotted_fit_2{i}{2},':','LineWidth',4,'Color',[0.5 0.5 0.5])
    set(gca,'XScale','log')
    ax = gca;
    ax.FontSize = 24;
    ax.YScale = 'log';
    ax.XLim = [-inf 5];
    ax.YLim = [1e-2,1];
    set(gca,'linewidth',3)
    ax.TickLength=[0.05,0.05];
    axis square
    grid on
    %exportgraphics(gcf, strcat('./three stim share figs/all_ff-isi-z1-',subregions_ff(i),'.png'),'Resolution',1500)
    xlabel('Interburst Interval (s)')
    ylabel('Pr(X \geq x)')
    xlim([1e-1 25]);
    ylim([1e-2 1])
    xticks([0.1,0.25,1,2.5,5,10,20,100])
    yticks([1e-2 1e-1, 0.25, 0.5, 1])%cdf
    set(gca,'YScale','log')
    set(gca,'XScale','log')
%     title(strcat(subregions_ff(i)," IBI Feedforward"))
    hold off
    exportgraphics(gcf, strcat('./three stim share figs/all_ff-ibi-',subregions_ff(i),'.png'),'Resolution',1500)
    
    prob_fifty_ISI_1{i}{1}(1,:)=[1];
    prob_fifty_ISI_1{i}{1}(2,:)=[2];
    prob_fifty_ISI_1{i}{1}(3,:)=[3];
    
    prob=nostim_plot.points_plotted{i}{2}(abs(nostim_plot.points_plotted{i}{1}-0.5)==min(abs(nostim_plot.points_plotted{i}{1}-0.5)));
    prob_fifty_ISI_1{i}{2}(1)=prob(1);
    prob=theta_plot.points_plotted{i}{2}(abs(theta_plot.points_plotted{i}{1}-0.5)==min(abs(theta_plot.points_plotted{i}{1}-0.5)));
    prob_fifty_ISI_1{i}{2}(2)=prob(1);
    prob=hfs_plot.points_plotted{i}{2}(abs(hfs_plot.points_plotted{i}{1}-0.5)==min(abs(hfs_plot.points_plotted{i}{1}-0.5)));
    prob_fifty_ISI_1{i}{2}(3)=prob(1);
    
    prob_fifty_ISI_2{i}{1}(1,:)=[1];
    prob_fifty_ISI_2{i}{1}(2,:)=[2];
    prob_fifty_ISI_2{i}{1}(3,:)=[3];
    
    prob=nostim_plot.points_plotted{i}{2}(abs(nostim_plot.points_plotted{i}{1}-0.5)==min(abs(nostim_plot.points_plotted{i}{1}-0.5)));
    prob_fifty_ISI_2{i}{2}(1)=prob(1);
    prob=theta_plot.points_plotted{i}{2}(abs(theta_plot.points_plotted{i}{1}-0.5)==min(abs(theta_plot.points_plotted{i}{1}-0.5)));
    prob_fifty_ISI_2{i}{2}(2)=prob(1);
    prob=hfs_plot.points_plotted{i}{2}(abs(hfs_plot.points_plotted{i}{1}-0.5)==min(abs(hfs_plot.points_plotted{i}{1}-0.5)));
    prob_fifty_ISI_2{i}{2}(3)=prob(1);
end

%% fb Bar graphs by subregion lower

nostim_plot=load('fb_IBI_points_plotted_nostim.mat');
nostim_fit_1=load('fb_IBI_fit_plotted_nostim_1.mat');
nostim_fit_2=load('fb_IBI_fit_plotted_nostim_2.mat');
theta_plot=load('fb_IBI_points_plotted_thetastim.mat');
theta_fit_1=load('fb_IBI_fit_plotted_thetastim_1.mat');
theta_fit_2=load('fb_IBI_fit_plotted_thetastim_2.mat');
hfs_plot=load('fb_IBI_points_plotted_HFS.mat');
hfs_fit_1=load('fb_IBI_fit_plotted_HFS_1.mat');
hfs_fit_2=load('fb_IBI_fit_plotted_HFS_2.mat');

subregions_ff=["EC-DG","DG-CA3","CA3-CA1","CA1-EC"];
subregions_fb=["DG-EC","CA3-DG","CA1-CA3","EC-CA1"];

prob_fifty_ISI_1=[];
prob_fifty_ISI_2=[];
for i=1:length(subregions_fb)
    figure( 'Position', [100 100 700 600])
    hold on
    plot(nostim_plot.points_plotted{i}{1},nostim_plot.points_plotted{i}{2},'Color',fb_colors(1),'LineWidth',2)
    plot(nostim_fit_1.plotted_fit_1{i}{1},nostim_fit_1.plotted_fit_1{i}{2},'--k','LineWidth',2)
    plot(nostim_fit_2.plotted_fit_2{i}{1},nostim_fit_2.plotted_fit_2{i}{2},':','LineWidth',3,'Color',[0.5 0.5 0.5])
    plot(theta_plot.points_plotted{i}{1},theta_plot.points_plotted{i}{2},'Color',fb_colors(2),'LineWidth',2)
    plot(theta_fit_1.plotted_fit_1{i}{1},theta_fit_1.plotted_fit_1{i}{2},'--k','LineWidth',2)
    plot(theta_fit_2.plotted_fit_2{i}{1},theta_fit_2.plotted_fit_2{i}{2},':','LineWidth',3,'Color',[0.5 0.5 0.5])
    plot(hfs_plot.points_plotted{i}{1},hfs_plot.points_plotted{i}{2},'Color',fb_colors(3),'LineWidth',2)
    plot(hfs_fit_1.plotted_fit_1{i}{1},hfs_fit_1.plotted_fit_1{i}{2},'--k','LineWidth',2)
    plot(hfs_fit_2.plotted_fit_2{i}{1},hfs_fit_2.plotted_fit_2{i}{2},':','LineWidth',3,'Color',[0.5 0.5 0.5])
    set(gca,'XScale','log')
    ax = gca;
    ax.FontSize = 24;
    ax.YScale = 'log';
    ax.XLim = [-inf 5];
    ax.YLim = [1e-2,1];
    set(gca,'linewidth',3)
    ax.TickLength=[0.05,0.05];
    axis square
    grid on
    %exportgraphics(gcf, strcat('./three stim share figs/all_fb-isi-z1-',subregions_fb(i),'.png'),'Resolution',1500)
    xlabel('Interburst Interval (s)')
    ylabel('Pr(X \geq x)')
    xlim([1e-1 25]);
    ylim([1e-2 1])
    xticks([0.1,0.25,1,2.5,5,10,20,100])
    yticks([1e-2 1e-1, 0.25, 0.5, 1])%cdf
    set(gca,'YScale','log')
    set(gca,'XScale','log')
%     title(strcat(subregions_fb(i)," IBI Feedback"))
    hold off
    exportgraphics(gcf, strcat('./three stim share figs/all_fb-ibi-',subregions_fb(i),'.png'),'Resolution',1500)
    
    prob_fifty_ISI_1{i}{1}(1,:)=[1,2];
    prob_fifty_ISI_1{i}{1}(2,:)=[1,3];
    prob_fifty_ISI_1{i}{1}(3,:)=[2,3];
    
    prob=nostim_plot.points_plotted{i}{2}(abs(nostim_plot.points_plotted{i}{1}-0.5)==min(abs(nostim_plot.points_plotted{i}{1}-0.5)));
    prob_fifty_ISI_1{i}{2}(1)=prob(1);
    prob=theta_plot.points_plotted{i}{2}(abs(theta_plot.points_plotted{i}{1}-0.5)==min(abs(theta_plot.points_plotted{i}{1}-0.5)));
    prob_fifty_ISI_1{i}{2}(2)=prob(1);
    prob=hfs_plot.points_plotted{i}{2}(abs(hfs_plot.points_plotted{i}{1}-0.5)==min(abs(hfs_plot.points_plotted{i}{1}-0.5)));
    prob_fifty_ISI_1{i}{2}(3)=prob(1);
    
    prob_fifty_ISI_2{i}{1}(1,:)=[1,2];
    prob_fifty_ISI_2{i}{1}(2,:)=[1,3];
    prob_fifty_ISI_2{i}{1}(3,:)=[2,3];
    
    prob=nostim_plot.points_plotted{i}{2}(abs(nostim_plot.points_plotted{i}{1}-0.5)==min(abs(nostim_plot.points_plotted{i}{1}-0.5)));
    prob_fifty_ISI_2{i}{2}(1)=prob(1);
    prob=theta_plot.points_plotted{i}{2}(abs(theta_plot.points_plotted{i}{1}-0.5)==min(abs(theta_plot.points_plotted{i}{1}-0.5)));
    prob_fifty_ISI_2{i}{2}(2)=prob(1);
    prob=hfs_plot.points_plotted{i}{2}(abs(hfs_plot.points_plotted{i}{1}-0.5)==min(abs(hfs_plot.points_plotted{i}{1}-0.5)));
    prob_fifty_ISI_2{i}{2}(3)=prob(1);
end

%% Well Bar graph by subregion

nostim_plot=load('well_IBI_points_plotted_nostim_5SD_500maxSD.mat');
nostim_fit_1=load('well_IBI_fit_plotted_nostim_1_5SD_500maxSD.mat');
nostim_fit_2=load('well_IBI_fit_plotted_nostim_2_5SD_500maxSD.mat');
theta_plot=load('well_IBI_points_plotted_theta_5SD_500maxSD.mat');
theta_fit_1=load('well_IBI_fit_plotted_theta_1_5SD_500maxSD.mat');
theta_fit_2=load('well_IBI_fit_plotted_theta_2_5SD_500maxSD.mat');
hfs_plot=load('well_IBI_points_plotted_HFS_5SD_500maxSD.mat');
hfs_fit_1=load('well_IBI_fit_plotted_HFS_1_5SD_500maxSD.mat');
hfs_fit_2=load('well_IBI_fit_plotted_HFS_2_5SD_500maxSD.mat');

subregions_well=["EC","DG","CA3","CA1"];

prob_fifty_ISI_1=[];
prob_fifty_ISI_2=[];
for i=1:length(subregions_well)
    figure( 'Position', [100 100 700 600])
    hold on
    plot(nostim_plot.points_plotted{i}{1},nostim_plot.points_plotted{i}{2},'Color',[0 0.4470 0.7410],'LineWidth',2)
    plot(nostim_fit_1.plotted_fit_1{i}{1},nostim_fit_1.plotted_fit_1{i}{2},'--k','LineWidth',2)
    plot(nostim_fit_2.plotted_fit_2{i}{1},nostim_fit_2.plotted_fit_2{i}{2},':','LineWidth',3,'Color',[0.5 0.5 0.5])
    plot(theta_plot.points_plotted{i}{1},theta_plot.points_plotted{i}{2},'Color',[0.8500 0.3250 0.0980],'LineWidth',2)
    plot(theta_fit_1.plotted_fit_1{i}{1},theta_fit_1.plotted_fit_1{i}{2},'--k','LineWidth',2)
    plot(theta_fit_2.plotted_fit_2{i}{1},theta_fit_2.plotted_fit_2{i}{2},':','LineWidth',3,'Color',[0.5 0.5 0.5])
    plot(hfs_plot.points_plotted{i}{1},hfs_plot.points_plotted{i}{2},'Color',[0.9290 0.6940 0.1250],'LineWidth',2)
    plot(hfs_fit_1.plotted_fit_1{i}{1},hfs_fit_1.plotted_fit_1{i}{2},'--k','LineWidth',2)
    plot(hfs_fit_2.plotted_fit_2{i}{1},hfs_fit_2.plotted_fit_2{i}{2},':','LineWidth',3,'Color',[0.5 0.5 0.5])
    set(gca,'XScale','log')
    ax = gca;
    ax.FontSize = 24;
    ax.YScale = 'log';
    ax.XLim = [-inf 5];
    ax.YLim = [1e-2,1];
    set(gca,'linewidth',3)
    ax.TickLength=[0.05,0.05];
    axis square
    grid on
    %exportgraphics(gcf, strcat('./three stim share figs/all_ff-isi-z1-',subregions_ff(i),'.png'),'Resolution',1500)
    xlabel('Interburst Interval (s)')
    ylabel('Pr(X \geq x)')
    xlim([1e-1 25]);
    ylim([1e-2 1])
    xticks([0.1,0.25,1,2.5,5,10,20,100])
    yticks([1e-2 1e-1, 0.25, 0.5, 1])%cdf
    set(gca,'YScale','log')
    set(gca,'XScale','log')
%     title(strcat(subregions_well(i)," IBI Well"))
    hold off
    exportgraphics(gcf, strcat('./three stim share figs/all_well_5SD_500maxSD-ibi-',subregions_well(i),'.png'),'Resolution',1500)
    
    prob_fifty_ISI_1{i}{1}(1,:)=[1];
    prob_fifty_ISI_1{i}{1}(2,:)=[2];
    prob_fifty_ISI_1{i}{1}(3,:)=[3];
    
    prob=nostim_plot.points_plotted{i}{2}(abs(nostim_plot.points_plotted{i}{1}-0.5)==min(abs(nostim_plot.points_plotted{i}{1}-0.5)));
    prob_fifty_ISI_1{i}{2}(1)=prob(1);
    prob=theta_plot.points_plotted{i}{2}(abs(theta_plot.points_plotted{i}{1}-0.5)==min(abs(theta_plot.points_plotted{i}{1}-0.5)));
    prob_fifty_ISI_1{i}{2}(2)=prob(1);
    prob=hfs_plot.points_plotted{i}{2}(abs(hfs_plot.points_plotted{i}{1}-0.5)==min(abs(hfs_plot.points_plotted{i}{1}-0.5)));
    prob_fifty_ISI_1{i}{2}(3)=prob(1);
    
    prob_fifty_ISI_2{i}{1}(1,:)=[1];
    prob_fifty_ISI_2{i}{1}(2,:)=[2];
    prob_fifty_ISI_2{i}{1}(3,:)=[3];
    
    prob=nostim_plot.points_plotted{i}{2}(abs(nostim_plot.points_plotted{i}{1}-0.5)==min(abs(nostim_plot.points_plotted{i}{1}-0.5)));
    prob_fifty_ISI_2{i}{2}(1)=prob(1);
    prob=theta_plot.points_plotted{i}{2}(abs(theta_plot.points_plotted{i}{1}-0.5)==min(abs(theta_plot.points_plotted{i}{1}-0.5)));
    prob_fifty_ISI_2{i}{2}(2)=prob(1);
    prob=hfs_plot.points_plotted{i}{2}(abs(hfs_plot.points_plotted{i}{1}-0.5)==min(abs(hfs_plot.points_plotted{i}{1}-0.5)));
    prob_fifty_ISI_2{i}{2}(3)=prob(1);
end

%% Calculate N
nostim=load('D:\Brewer lab data\HFS\No Stim\02-Oct-2022_A\spike_burst_dyn_table_stim.mat'); % 5 SD
hfs5=load('D:\Brewer lab data\HFS\Theta Stim\27-Sep-2022_A\spike_burst_dyn_table_stim.mat'); % 5SD
hfs40=load('D:\Brewer lab data\HFS\HFS Stim\03-Oct-2022_A\spike_burst_dyn_table_stim.mat'); % 5 SD

nostim_x_poi = [4.78341,2.12803,1.88458,2.82549];
hfs5_x_poi=[2.30757,2.94227,1.1132,0.873059]; %5SD
hfs40_x_poi = [1.96247,2.12803,2.21598,3.32235];

n_array=[9,6,6];
for i=1:4
    low_bursts_nostim_ff=cell2mat(nostim.spike_burst_dyn_table_stim.IBI(nostim.spike_burst_dyn_table_stim.regi==i....
        & nostim.spike_burst_dyn_table_stim.if_ff==1));
    low_bursts_hfs5_ff=cell2mat(hfs5.spike_burst_dyn_table_stim.IBI(hfs5.spike_burst_dyn_table_stim.regi==i....
        & hfs5.spike_burst_dyn_table_stim.if_ff==1));
    low_bursts_hfs40_ff=cell2mat(hfs40.spike_burst_dyn_table_stim.IBI(hfs40.spike_burst_dyn_table_stim.regi==i....
        & hfs40.spike_burst_dyn_table_stim.if_ff==1));
    n_low_ff{i}=[length(low_bursts_nostim_ff(low_bursts_nostim_ff<=nostim_x_poi(i))),...
        length(low_bursts_hfs5_ff(low_bursts_hfs5_ff<=hfs5_x_poi(i))),length(low_bursts_hfs40_ff(low_bursts_hfs40_ff<=hfs40_x_poi(i)))];
    n_low_ff_norm{i}=n_low_ff{i}./n_array;

    low_bursts_nostim_fb=cell2mat(nostim.spike_burst_dyn_table_stim.IBI(nostim.spike_burst_dyn_table_stim.regi==i....
        & nostim.spike_burst_dyn_table_stim.if_ff==0));
    low_bursts_hfs5_fb=cell2mat(hfs5.spike_burst_dyn_table_stim.IBI(hfs5.spike_burst_dyn_table_stim.regi==i....
        & hfs5.spike_burst_dyn_table_stim.if_ff==0));
    low_bursts_hfs40_fb=cell2mat(hfs40.spike_burst_dyn_table_stim.IBI(hfs40.spike_burst_dyn_table_stim.regi==i....
        & hfs40.spike_burst_dyn_table_stim.if_ff==0));
    n_low_fb{i}=[length(low_bursts_nostim_fb(low_bursts_nostim_fb<=nostim_x_poi(i))),...
        length(low_bursts_hfs5_fb(low_bursts_hfs5_fb<=hfs5_x_poi(i))),length(low_bursts_hfs40_fb(low_bursts_hfs40_fb<=hfs40_x_poi(i)))];
    n_low_fb_norm{i}=n_low_fb{i}./n_array;

    hi_bursts_nostim_ff=cell2mat(nostim.spike_burst_dyn_table_stim.IBI(nostim.spike_burst_dyn_table_stim.regi==i....
        & nostim.spike_burst_dyn_table_stim.if_ff==1));
    hi_bursts_hfs5_ff=cell2mat(hfs5.spike_burst_dyn_table_stim.IBI(hfs5.spike_burst_dyn_table_stim.regi==i....
        & hfs5.spike_burst_dyn_table_stim.if_ff==1));
    hi_bursts_hfs40_ff=cell2mat(hfs40.spike_burst_dyn_table_stim.IBI(hfs40.spike_burst_dyn_table_stim.regi==i....
        & hfs40.spike_burst_dyn_table_stim.if_ff==1));
    n_hi_ff{i}=[length(hi_bursts_nostim_ff(hi_bursts_nostim_ff>nostim_x_poi(i))),...
        length(hi_bursts_hfs5_ff(hi_bursts_hfs5_ff>hfs5_x_poi(i))),length(hi_bursts_hfs40_ff(hi_bursts_hfs40_ff>hfs40_x_poi(i)))];
    n_hi_ff_norm{i}=n_hi_ff{i}./n_array;

    hi_bursts_nostim_fb=cell2mat(nostim.spike_burst_dyn_table_stim.IBI(nostim.spike_burst_dyn_table_stim.regi==i....
        & nostim.spike_burst_dyn_table_stim.if_ff==0));
    hi_bursts_hfs5_fb=cell2mat(hfs5.spike_burst_dyn_table_stim.IBI(hfs5.spike_burst_dyn_table_stim.regi==i....
        & hfs5.spike_burst_dyn_table_stim.if_ff==0));
    hi_bursts_hfs40_fb=cell2mat(hfs40.spike_burst_dyn_table_stim.IBI(hfs40.spike_burst_dyn_table_stim.regi==i....
        & hfs40.spike_burst_dyn_table_stim.if_ff==0));
    n_hi_fb{i}=[length(hi_bursts_nostim_fb(hi_bursts_nostim_fb>nostim_x_poi(i))),...
        length(hi_bursts_hfs5_fb(hi_bursts_hfs5_fb>hfs5_x_poi(i))),length(hi_bursts_hfs40_fb(hi_bursts_hfs40_fb>hfs40_x_poi(i)))];
    n_hi_fb_norm{i}=n_hi_fb{i}./n_array;
end

%% lower vs upper ff

%normalized ff and fb slopes
for i=1:4
    ff_lower_norm(:,i)=ff_lower(:,i)./abs(max(ff_lower(:,i)));
    ff_upper_norm(:,i)=ff_upper(:,i)./abs(max(ff_upper(:,i)));
end

lo_hi_delta=(abs(ff_lower_norm)-abs(ff_upper_norm))./(abs(ff_lower_norm)+abs(ff_upper_norm));

labels={'EC-DG','DG-CA3','CA3-CA1','CA1-EC'};
regs=categorical(labels);
regs=reordercats(regs,labels);

figure( 'Position', [100 100 700 600])
% b=barh(regs,ff_fb_perc);

b=plot_horzbar_CohenD(lo_hi_delta',zeros(size(lo_hi_delta')),zeros(size(lo_hi_delta')),labels);
b(1).BaseValue=0;
b(1).BaseLine.LineStyle = "--";
b(1).BaseLine.LineWidth=2;
xlim([-1,1])
xticks([-2:0.2:2])
xticklabels([-2:0.2:2])

xlabel("<- Slow Bursting            Fast Bursting ->")
set(gca,'LineWidth',2,'TickLength',[0.02 0.02]);
set(gca,'FontSize',18)
exportgraphics(gcf,'.\three stim share figs\ff_slow_fast_IBI_perc.png','Resolution',1500)
%% lower vs upper fb

%normalized ff and fb slopes
for i=1:4
    fb_lower_norm(:,i)=fb_lower(:,i)./abs(max(fb_lower(:,i)));
    fb_upper_norm(:,i)=fb_upper(:,i)./abs(max(fb_upper(:,i)));
end

lo_hi_delta=(abs(fb_lower_norm)-abs(fb_upper_norm))./(abs(fb_lower_norm)+abs(fb_upper_norm));

labels={'EC-DG','DG-CA3','CA3-CA1','CA1-EC'};
regs=categorical(labels);
regs=reordercats(regs,labels);

figure( 'Position', [100 100 700 600])
% b=barh(regs,ff_fb_perc);

b=plot_horzbar_CohenD(lo_hi_delta',zeros(size(lo_hi_delta')),zeros(size(lo_hi_delta')),labels);
b(1).BaseValue=0;
b(1).BaseLine.LineStyle = "--";
b(1).BaseLine.LineWidth=2;
xlim([-1,1])
xticks([-2:0.2:2])
xticklabels([-2:0.2:2])

xlabel("<- Slow Bursting            Fast Bursting ->")
set(gca,'LineWidth',2,'TickLength',[0.02 0.02]);
set(gca,'FontSize',18)
exportgraphics(gcf,'.\three stim share figs\fb_slow_fast_IBI_perc.png','Resolution',1500)

%% lower vs upper well

%normalized well slopes
for i=1:4
    well_lower_norm(:,i)=well_lower(:,i)./abs(max(well_lower(:,i)));
    well_upper_norm(:,i)=well_upper(:,i)./abs(max(well_upper(:,i)));
end

lo_hi_delta=(abs(well_lower_norm)-abs(well_upper_norm))./(abs(well_lower_norm)+abs(well_upper_norm));

labels={'EC','DG','CA3','CA1'};
regs=categorical(labels);
regs=reordercats(regs,labels);

figure( 'Position', [100 100 700 600])
% b=barh(regs,ff_well_perc);

b=plot_horzbar_CohenD(lo_hi_delta',zeros(size(lo_hi_delta')),zeros(size(lo_hi_delta')),labels);
b(1).BaseValue=0;
b(1).BaseLine.LineStyle = "--";
b(1).BaseLine.LineWidth=2;
xlim([-1,1])
xticks([-2:0.2:2])
xticklabels([-2:0.2:2])

xlabel("<- Slow Bursting            Fast Bursting ->")
set(gca,'LineWidth',2,'TickLength',[0.02 0.02]);
set(gca,'FontSize',18)
exportgraphics(gcf,'.\three stim share figs\well_slow_fast_IBI_perc.png','Resolution',1500)

%% Functions
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

function bff=plot_horzbar_CohenD...
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