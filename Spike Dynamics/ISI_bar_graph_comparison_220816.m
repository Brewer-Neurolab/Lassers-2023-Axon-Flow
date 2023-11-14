%% ISI comparison bar graphs and anova
clear 
clc
load('D:\Brewer lab data\HFS\No Stim\02-Oct-2022_A\spike_burst_dyn_table_stim.mat') %5 SD
nostim_table=spike_burst_dyn_table_stim;
load('D:\Brewer lab data\HFS\Theta Stim\27-Sep-2022_A\spike_burst_dyn_table_stim.mat') % 5SD
hfs5_table=spike_burst_dyn_table_stim;
load('D:\Brewer lab data\HFS\HFS Stim\03-Oct-2022_A\spike_burst_dyn_table_stim.mat') % 5 SD
hfs40_table=spike_burst_dyn_table_stim;

ff_colors=["#3869EE","#FF2303","#EE9B00"];
fb_colors=["#4fb7ff","#FF9B9B","#FFDA49"];

n_array=[9;6;6];
%% Feedforward

nostim=load('no_stim_ff_ancova_vars.mat');
theta=load('theta_stim_ff_ancova_vars.mat');
hfs=load('HFS_ff_ancova_vars.mat');

subregions_ff=["EC-DG","DG-CA3","CA3-CA1","CA1-EC"];
ttest_struct=[];
perc_change=[];
for i=1:length(subregions_ff)
    nostim_idx=contains(string(nostim.regLabel_post_ff),subregions_ff(i));
    theta_idx=contains(string(theta.regLabel_post_ff),subregions_ff(i));
    hfs_idx=contains(string(hfs.regLabel_post_ff),subregions_ff(i));
    sp_feature=[nostim.sp_feature_post_ff(nostim_idx);theta.sp_feature_post_ff(theta_idx);hfs.sp_feature_post_ff(hfs_idx)];
    prob_vec=[nostim.prob_vec_post_ff(nostim_idx);theta.prob_vec_post_ff(theta_idx);hfs.prob_vec_post_ff(hfs_idx)];
    regLabel=[nostim.regLabel_post_ff(nostim_idx);theta.regLabel_post_ff(theta_idx);hfs.regLabel_post_ff(hfs_idx)];
    n_ff{i}=[length(nostim.regLabel_post_ff(nostim_idx));length(theta.regLabel_post_ff(theta_idx));length(hfs.regLabel_post_ff(hfs_idx))];

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
    [~,pval]=ttest2(nostim.prob_vec_post_ff(nostim_idx),theta.prob_vec_post_ff(theta_idx));
    ttest_struct{i}.pval(1)=pval;
    [~,pval]=ttest2(nostim.prob_vec_post_ff(nostim_idx),hfs.prob_vec_post_ff(hfs_idx));
    ttest_struct{i}.pval(2)=pval;
    [~,pval]=ttest2(theta.prob_vec_post_ff(theta_idx),hfs.prob_vec_post_ff(hfs_idx));
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
ax.FontSize=30;
%legend('No Stimulation','Theta Stimulation','HFS Stimulation')
% title("Feedforward ISI")
ylabel 'Slope of ISI (s^{-1})'
ylim([-1,0])
hold off
ax.TickLength=[.02,.02];
ax.LineWidth=2;
% exportgraphics(gcf,'.\three stim share figs\ff_isi_bar.png','Resolution',1500)


for i=1:length(subregions_ff)
    ax.XAxis.Limits=categorical({char(subregions_ff(i)),char(subregions_ff(i))});
%     exportgraphics(gcf,strcat('.\three stim share figs\ff_isi_bar',subregions_ff(i),'.png'),'Resolution',1500)
end

% n isi
n_isi_ff_nostim=[];
n_isi_ff_hfs5=[];
n_isi_ff_hfs40=[];
for i=1:length(subregions_ff)
    n_isi_ff_nostim(i)=length(cell2mat(nostim_table.ISI(nostim_table.regi==i & nostim_table.if_ff==1)));
    n_isi_ff_hfs5(i)=length(cell2mat(hfs5_table.ISI(hfs5_table.regi==i & hfs5_table.if_ff==1)));
    n_isi_ff_hfs40(i)=length(cell2mat(hfs40_table.ISI(hfs40_table.regi==i & hfs40_table.if_ff==1)));
end
%% Feedback

nostim=load('no_stim_fb_ancova_vars.mat');
theta=load('theta_stim_fb_ancova_vars.mat');
hfs=load('HFS_fb_ancova_vars.mat');
ttest_struct=[];
subregions_fb=["EC-DG","DG-CA3","CA3-CA1","CA1-EC"];
figure( 'Position', [100 100 700 600])
for i=1:length(subregions_fb)
    nostim_idx=contains(string(nostim.regLabel_post_fb),subregions_fb(i));
    theta_idx=contains(string(theta.regLabel_post_fb),subregions_fb(i));
    hfs_idx=contains(string(hfs.regLabel_post_fb),subregions_fb(i));
    sp_feature=[nostim.sp_feature_post_fb(nostim_idx);theta.sp_feature_post_fb(theta_idx);hfs.sp_feature_post_fb(hfs_idx)];
    prob_vec=[nostim.prob_vec_post_fb(nostim_idx);theta.prob_vec_post_fb(theta_idx);hfs.prob_vec_post_fb(hfs_idx)];
    regLabel=[nostim.regLabel_post_fb(nostim_idx);theta.regLabel_post_fb(theta_idx);hfs.regLabel_post_fb(hfs_idx)];
    n_fb{i}=[length(nostim.regLabel_post_fb(nostim_idx));length(theta.regLabel_post_fb(theta_idx));length(hfs.regLabel_post_fb(hfs_idx))];

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
    [~,pval]=ttest2(nostim.prob_vec_post_fb(nostim_idx),theta.prob_vec_post_fb(theta_idx));
    ttest_struct{i}.pval(1)=pval;
    [~,pval]=ttest2(nostim.prob_vec_post_fb(nostim_idx),hfs.prob_vec_post_fb(hfs_idx));
    ttest_struct{i}.pval(2)=pval;
    [~,pval]=ttest2(theta.prob_vec_post_fb(theta_idx),hfs.prob_vec_post_fb(hfs_idx));
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
errorbar(x',ydata,errordata,'k','linestyle','none');
ax=gca;
ax.FontSize=30;
%legend('No Stimulation','Theta Stimulation','HFS Stimulation')
% title("Feedback ISI")
ylabel 'Slope of ISI (s^{-1})'
ylim([-1,0])
hold off
ax.TickLength=[.02,.02];
ax.LineWidth=2;
% exportgraphics(gcf,'.\three stim share figs\fb_isi_bar.png','Resolution',1500)

for i=1:length(subregions_fb)
    ax.XAxis.Limits=categorical({char(subregions_fb(i)),char(subregions_fb(i))});
%     exportgraphics(gcf,strcat('.\three stim share figs\fb_isi_bar',subregions_ff(i),'.png'),'Resolution',1500)
end

% n isi
n_isi_fb_nostim=[];
n_isi_fb_hfs5=[];
n_isi_fb_hfs40=[];
for i=1:length(subregions_ff)
    n_isi_fb_nostim(i)=length(cell2mat(nostim_table.ISI(nostim_table.regi==i & nostim_table.if_ff==0)));
    n_isi_fb_hfs5(i)=length(cell2mat(hfs5_table.ISI(hfs5_table.regi==i & hfs5_table.if_ff==0)));
    n_isi_fb_hfs40(i)=length(cell2mat(hfs40_table.ISI(hfs40_table.regi==i & hfs40_table.if_ff==0)));
end
%% Horzbar Slopes

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
xlabel("FB Slope of ISI (s^{-1})                                      FF Slope of ISI (s^{-1})")
xlim([-1.2,1.2])
xticks([-1.5:.1:1.5])
% ylim([0,5])
set(gca,'LineWidth',2,'TickLength',[0.02 0.02]);
set(gca,'FontSize',12)

% exportgraphics(gcf,strcat('.\three stim share figs\ff_fb_isi_horzbar.png'),'Resolution',1500)


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
%     exportgraphics(gcf,strcat('.\three stim share figs\ff_fb_isi_horzbar',labels{i},'.png'),'Resolution',1500)
end
%% FF FB ANCOVA

nostim_ff=load("no_stim_ff_ancova_vars.mat");
nostim_fb=load("no_stim_fb_ancova_vars.mat");
nostim_reg=unique(nostim_ff.regLabel_post_ff);
hfs5_ff=load("theta_stim_ff_ancova_vars.mat");
hfs5_fb=load("theta_stim_fb_ancova_vars.mat");
hfs40_ff=load("hfs_ff_ancova_vars.mat");
hfs40_fb=load("hfs_fb_ancova_vars.mat");

for i=1:4
    
    nostim_sp_feat=[nostim_ff.sp_feature_post_ff]
%     [~,~,~,stats]

end
%% FF/FB% 

ff_fb_perc=((dat2plot_ff-dat2plot_fb)./(dat2plot_ff+dat2plot_fb))*100;

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
xlim([-45,45])
xticks([-40:10:40])
xticklabels([-40:10:40])

xlabel("<- Excess FB            Excess FF ->")
set(gca,'LineWidth',2,'TickLength',[0.02 0.02]);
set(gca,'FontSize',18)
% exportgraphics(gcf,'.\three stim share figs\ff_fb_ISI_perc.png','Resolution',1500)
%% Wells

nostim=load('no_stim_well_ancova_vars_5SD_500maxSD.mat');
theta=load('theta_well_ancova_vars_5SD_500maxSD.mat');
hfs=load('HFS_well_ancova_vars_5SD_500maxSD.mat');

subregions_well=["EC","DG","CA3","CA1"];
ttest_struct=[];
perc_change=[];
for i=1:length(subregions_well)
    nostim_idx=contains(string(nostim.regLabel_post_well),subregions_well(i));
    theta_idx=contains(string(theta.regLabel_post_well),subregions_well(i));
    hfs_idx=contains(string(hfs.regLabel_post_well),subregions_well(i));
    sp_feature=[nostim.sp_feature_post_well(nostim_idx);theta.sp_feature_post_well(theta_idx);hfs.sp_feature_post_well(hfs_idx)];
    prob_vec=[nostim.prob_vec_post_well(nostim_idx);theta.prob_vec_post_well(theta_idx);hfs.prob_vec_post_well(hfs_idx)];
    regLabel=[nostim.regLabel_post_well(nostim_idx);theta.regLabel_post_well(theta_idx);hfs.regLabel_post_well(hfs_idx)];
    n_well{i}=[length(nostim.regLabel_post_well(nostim_idx));length(theta.regLabel_post_well(theta_idx));length(hfs.regLabel_post_well(hfs_idx))];

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
    [~,pval]=ttest2(nostim.prob_vec_post_well(nostim_idx),theta.prob_vec_post_well(theta_idx));
    ttest_struct{i}.pval(1)=pval;
    [~,pval]=ttest2(nostim.prob_vec_post_well(nostim_idx),hfs.prob_vec_post_well(hfs_idx));
    ttest_struct{i}.pval(2)=pval;
    [~,pval]=ttest2(theta.prob_vec_post_well(theta_idx),hfs.prob_vec_post_well(hfs_idx));
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
ax.FontSize=30;
%legend('No Stimulation','Theta Stimulation','HFS Stimulation')
% title("Well ISI")
ylabel 'Slope of ISI (s^{-1})'
ylim([-1,0])
hold off
ax.TickLength=[.02,.02];
ax.LineWidth=2;
exportgraphics(gcf,'.\three stim share figs\well_isi_bar_5SD_500maxSD.png','Resolution',1500)


for i=1:length(subregions_well)
    ax.XAxis.Limits=categorical({char(subregions_well(i)),char(subregions_well(i))});
    exportgraphics(gcf,strcat('.\three stim share figs\well_isi_bar_5SD_500maxSD',subregions_well(i),'.png'),'Resolution',1500)
end
%% Bar graphs by subregion ff

nostim_plot=load('ff_ISI_points_plotted_nostim.mat');
nostim_fit=load('ff_ISI_fit_plotted_nostim.mat');
theta_plot=load('ff_ISI_points_plotted_thetastim.mat');
theta_fit=load('ff_ISI_fit_plotted_thetastim.mat');
hfs_plot=load('ff_ISI_points_plotted_HFS.mat');
hfs_fit=load('ff_ISI_fit_plotted_HFS.mat');

subregions_ff=["EC-DG","DG-CA3","CA3-CA1","CA1-EC"];
subregions_fb=["DG-EC","CA3-DG","CA1-CA3","EC-CA1"];

prob_fifty_ISI=[];
for i=1:length(subregions_ff)
    figure( 'Position', [100 100 700 600])
    hold on
    plot(nostim_plot.points_plotted{i}{2},nostim_plot.points_plotted{i}{1},'Color',ff_colors(1),'LineWidth',4)
    plot(nostim_fit.plotted_fit{i}{1},nostim_fit.plotted_fit{i}{2},'--k','LineWidth',4)
    plot(theta_plot.points_plotted{i}{2},theta_plot.points_plotted{i}{1},'Color',ff_colors(2),'LineWidth',4)
    plot(theta_fit.plotted_fit{i}{1},theta_fit.plotted_fit{i}{2},'--k','LineWidth',4)
    plot(hfs_plot.points_plotted{i}{2},hfs_plot.points_plotted{i}{1},'Color',ff_colors(3),'LineWidth',4)
    plot(hfs_fit.plotted_fit{i}{1},hfs_fit.plotted_fit{i}{2},'--k','LineWidth',4)
    set(gca,'XScale','log')
    ax = gca;
    ax.FontSize = 40;
    ax.YScale = 'log';
    ax.XLim = [1e-2 5];
    ax.YLim = [1e-2,1];
    set(gca,'linewidth',3)
    ax.TickLength=[0.05,0.05];
    axis square
%     xticks(round(logspace(-2,0.6990,5),2))
    xticks([0.01,0.1,0.5,1,5])
    grid on
    exportgraphics(gcf, strcat('./three stim share figs/all_ff-isi-z1-',subregions_ff(i),'.png'),'Resolution',1500)
    ax.FontSize = 30;
    xlabel('Interspike Interval (s)')
    ylabel('Pr(X \geq x)')
    xlim([1e-2 2e-1]);
    ylim([1e-1 1])
    xticks([1e-2, 0.25e-1,0.5e-1,0.1,0.2,0.3,0.6e0,1,2])
    yticks([1e-1, 0.25, 0.5, 1])
    set(gca,'YScale','log')
    set(gca,'XScale','log')
%     title(strcat(subregions_ff(i)," ISI Feedforward"))
    hold off
    exportgraphics(gcf, strcat('./three stim share figs/all_ff-isi-',subregions_ff(i),'.png'),'Resolution',1500)
    
    prob_fifty_ISI{i}{1}(1,:)=[1,2];
    prob_fifty_ISI{i}{1}(2,:)=[1,3];
    prob_fifty_ISI{i}{1}(3,:)=[2,3];
    
    prob_fifty_ISI{i}{2}(1)=nostim_plot.points_plotted{i}{2}(abs(nostim_plot.points_plotted{i}{1}-0.5)==min(abs(nostim_plot.points_plotted{i}{1}-0.5)));
    prob_fifty_ISI{i}{2}(2)=theta_plot.points_plotted{i}{2}(abs(theta_plot.points_plotted{i}{1}-0.5)==min(abs(theta_plot.points_plotted{i}{1}-0.5)));
    prob_fifty_ISI{i}{2}(3)=hfs_plot.points_plotted{i}{2}(abs(hfs_plot.points_plotted{i}{1}-0.5)==min(abs(hfs_plot.points_plotted{i}{1}-0.5)));
end

%% fb Bar graphs by subregion

nostim_plot=load('fb_ISI_points_plotted_nostim.mat');
nostim_fit=load('fb_ISI_fit_plotted_nostim.mat');
theta_plot=load('fb_ISI_points_plotted_thetastim.mat');
theta_fit=load('fb_ISI_fit_plotted_thetastim.mat');
hfs_plot=load('fb_ISI_points_plotted_HFS.mat');
hfs_fit=load('fb_ISI_fit_plotted_HFS.mat');
prob_fifty_ISI=[];
for i=1:length(subregions_fb)
    figure( 'Position', [100 100 700 600])
    hold on
    plot(nostim_plot.points_plotted{i}{2},nostim_plot.points_plotted{i}{1},'Color',fb_colors(1),'LineWidth',4)
    plot(nostim_fit.plotted_fit{i}{1},nostim_fit.plotted_fit{i}{2},'--k','LineWidth',4)
    plot(theta_plot.points_plotted{i}{2},theta_plot.points_plotted{i}{1},'Color',fb_colors(2),'LineWidth',4)
    plot(theta_fit.plotted_fit{i}{1},theta_fit.plotted_fit{i}{2},'--k','LineWidth',4)
    plot(hfs_plot.points_plotted{i}{2},hfs_plot.points_plotted{i}{1},'Color',fb_colors(3),'LineWidth',4)
    plot(hfs_fit.plotted_fit{i}{1},hfs_fit.plotted_fit{i}{2},'--k','LineWidth',4)
    set(gca,'XScale','log')
    ax = gca;
    ax.FontSize = 40;
    ax.YScale = 'log';
    ax.XLim = [1e-2 5];
    ax.YLim = [1e-2,1];
    set(gca,'linewidth',3)
    ax.TickLength=[0.05,0.05];
    axis square
    xticks([0.01,0.1,0.5,1,5])
    grid on
    exportgraphics(gcf, strcat('./three stim share figs/all_fb-isi-z1-',subregions_fb(i),'.png'),'Resolution',1500)
    ax.FontSize = 30;
    xlabel('Interspike Interval (s)')
    ylabel('Pr(X \geq x)')
    xlim([1e-2 2e-1]);
    ylim([1e-1 1])
    xticks([1e-2, 0.25e-1,0.5e-1,0.1,0.2,0.3,0.6e0,1,2])
    yticks([1e-1, 0.25, 0.5, 1])
    set(gca,'YScale','log')
    set(gca,'XScale','log')
%     title(strcat(subregions_fb(i)," ISI Feedback"))
    hold off
    exportgraphics(gcf, strcat('./three stim share figs/all_fb-isi-',subregions_fb(i),'.png'),'Resolution',1500)
    
    prob_fifty_ISI{i}{1}(1,:)=[1,2];
    prob_fifty_ISI{i}{1}(2,:)=[1,3];
    prob_fifty_ISI{i}{1}(3,:)=[2,3];
    
    prob=nostim_plot.points_plotted{i}{2}(abs(nostim_plot.points_plotted{i}{1}-0.5)==min(abs(nostim_plot.points_plotted{i}{1}-0.5)));
    prob_fifty_ISI{i}{2}(1)=prob(1);
    prob=theta_plot.points_plotted{i}{2}(abs(theta_plot.points_plotted{i}{1}-0.5)==min(abs(theta_plot.points_plotted{i}{1}-0.5)));
    prob_fifty_ISI{i}{2}(2)=prob(1);
    prob=hfs_plot.points_plotted{i}{2}(abs(hfs_plot.points_plotted{i}{1}-0.5)==min(abs(hfs_plot.points_plotted{i}{1}-0.5)));
    prob_fifty_ISI{i}{2}(3)=prob(1);
end

%% Well Bar Graphs
nostim_plot=load('well_ISI_points_plotted_nostim_5SD_500maxSD.mat');
nostim_fit=load('well_ISI_fit_plotted_nostim_5SD_500maxSD.mat');
theta_plot=load('well_ISI_points_plotted_theta_5SD_500maxSD.mat');
theta_fit=load('well_ISI_fit_plotted_theta_5SD_500maxSD.mat');
hfs_plot=load('well_ISI_points_plotted_HFS_5SD_500maxSD.mat');
hfs_fit=load('well_ISI_fit_plotted_HFS_5SD_500maxSD.mat');
prob_fifty_ISI=[];
subregions_well=["EC","DG","CA3","CA1"];
for i=1:length(subregions_well)
    figure( 'Position', [100 100 700 600])
    hold on
    plot(nostim_plot.points_plotted{i}{2},nostim_plot.points_plotted{i}{1},'Color',[0 0.4470 0.7410],'LineWidth',4)
    plot(nostim_fit.plotted_fit{i}{1},nostim_fit.plotted_fit{i}{2},'--k','LineWidth',4)
    plot(theta_plot.points_plotted{i}{2},theta_plot.points_plotted{i}{1},'Color',[0.8500 0.3250 0.0980],'LineWidth',4)
    plot(theta_fit.plotted_fit{i}{1},theta_fit.plotted_fit{i}{2},'--k','LineWidth',4)
    plot(hfs_plot.points_plotted{i}{2},hfs_plot.points_plotted{i}{1},'Color',[0.9290 0.6940 0.1250],'LineWidth',4)
    plot(hfs_fit.plotted_fit{i}{1},hfs_fit.plotted_fit{i}{2},'--k','LineWidth',4)
    set(gca,'XScale','log')
    ax = gca;
    ax.FontSize = 40;
    ax.YScale = 'log';
    ax.XLim = [1e-2 5];
    ax.YLim = [1e-2,1];
    set(gca,'linewidth',3)
    ax.TickLength=[0.05,0.05];
    axis square
    xticks([0.01,0.1,0.5,1,5])
    grid on
    exportgraphics(gcf, strcat('./three stim share figs/all_well_5SD_500maxSD-isi-z1-',subregions_well(i),'.png'),'Resolution',1500)
    ax.FontSize = 30;
    xlabel('Interspike Interval (s)')
    ylabel('Pr(X \geq x)')
    xlim([1e-2 2e-1]);
    ylim([1e-1 1])
    xticks([1e-2, 0.25e-1,0.5e-1,0.1,0.2,0.3,0.6e0,1,2])
    yticks([1e-1, 0.25, 0.5, 1])
    set(gca,'YScale','log')
    set(gca,'XScale','log')
%     title(strcat(subregions_well(i)," ISI Well"))
    hold off
    exportgraphics(gcf, strcat('./three stim share figs/all_well_5SD_500maxSD-isi-',subregions_well(i),'.png'),'Resolution',1500)
    
    prob_fifty_ISI{i}{1}(1,:)=[1,2];
    prob_fifty_ISI{i}{1}(2,:)=[1,3];
    prob_fifty_ISI{i}{1}(3,:)=[2,3];
    
    prob=nostim_plot.points_plotted{i}{2}(abs(nostim_plot.points_plotted{i}{1}-0.5)==min(abs(nostim_plot.points_plotted{i}{1}-0.5)));
    prob_fifty_ISI{i}{2}(1)=prob(1);
    prob=theta_plot.points_plotted{i}{2}(abs(theta_plot.points_plotted{i}{1}-0.5)==min(abs(theta_plot.points_plotted{i}{1}-0.5)));
    prob_fifty_ISI{i}{2}(2)=prob(1);
    prob=hfs_plot.points_plotted{i}{2}(abs(hfs_plot.points_plotted{i}{1}-0.5)==min(abs(hfs_plot.points_plotted{i}{1}-0.5)));
    prob_fifty_ISI{i}{2}(3)=prob(1);
end
%% Functions
function plot_horzbar_from_anova...
    (dat2plot_ff,dat2plot_ff_err,dat2plot_fb,dat2plot_fb_err,labels,ff_colors,fb_colors)

% labels_cat=categorical(labels);
% % flip_labels=string(flip(labels));
% labels_cat=reordercats(labels_cat,labels);
% labels_cat=flip(labels_cat);

dat2plot_ff=flip(flip(dat2plot_ff,1),2);
dat2plot_ff_err=flip(flip(dat2plot_ff_err,1),2);
dat2plot_fb=flip(flip(dat2plot_fb,1),2);
dat2plot_fb_err=flip(flip(dat2plot_fb_err,1),2);

% labels=flip(categorical(labels));

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
