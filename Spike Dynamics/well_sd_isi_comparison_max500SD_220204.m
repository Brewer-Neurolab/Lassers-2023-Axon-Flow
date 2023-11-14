%% nostim
clear
clc
load('D:\Brewer lab data\HFS\No Stim\Wells_5SD\well_spike_dynamics_table_nostim.mat')
well_spike_burst_dyn_table = well_spike_dynamics_table;
clear spike_burst_dyn_table
well_regList = {'EC','DG','CA3','CA1'}; 
figOrder = [1 2 4 3 ];
col_array = ["#D95319","#77AC30","#4DBEEE","#7E2F8E"];
col_array_rgb = [0.8500 0.3250 0.0980; 0.4660 0.6740 0.1880; ...
    0.3010 0.7450 0.9330; 0.4940 0.1840 0.5560];
powerlawfitfun = @(b,x) 10.^(b(2)*log10(x) + b(1)); 

%% find n
ec_n=vertcat(well_spike_burst_dyn_table(well_spike_burst_dyn_table.regi==1,:).ISI{:});
dg_n=vertcat(well_spike_burst_dyn_table(well_spike_burst_dyn_table.regi==2,:).ISI{:});
ca3_n=vertcat(well_spike_burst_dyn_table(well_spike_burst_dyn_table.regi==3,:).ISI{:});
ca1_n=vertcat(well_spike_burst_dyn_table(well_spike_burst_dyn_table.regi==4,:).ISI{:});

%% ISI
fit_table = table([],[],[],[],'variablenames',{'regi','slope','intercept','rsquared'}); rowi=1;
figure( 'Position', [100 100 700 600])
binEdge = logspace(-2,1,100); binCenter = convert_edges_2_centers(binEdge);
xlimits = [ 1e-2, 2e-1];
sp_feature=[]; prob_vec=[]; regLabel=[]; pl=[];
hold on
disp('CCDF for Interspike Interval for FF axons')
for regi=1:4
    vec2plot = cell2mat(well_spike_burst_dyn_table(well_spike_burst_dyn_table.regi==regi,:).ISI);
    [ hist_object, hist_plot, Prob, xbin] = log_binned_histogram( vec2plot, binEdge, "cdf" );  
    hist_plot.Color = col_array(regi);
    pl(regi)=hist_plot;
    
    % fitting data between [xmin, xmax]
    [ coeff, stats, opt_lim, logBinC,logHistProb ] = find_powerlawfit_with_grid_search(binCenter, ...
        Prob, xlimits);
    sp_feature = [sp_feature; logBinC']; prob_vec=[prob_vec; logHistProb']; 
    regLabel = [regLabel; repmat(categorical(well_regList(regi)),length(logBinC),1)];
%     fitbinCenter = binCenter( binCenter > opt_lim(1) & binCenter< opt_lim(end));
    fitbinCenter = binCenter;
    regfit = powerlawfitfun(coeff,fitbinCenter);
    plot(fitbinCenter , regfit ,'--k')
    
    fit_table.regi(rowi) = regi;
    fit_table.intercept(rowi) = coeff(1);
    fit_table.slope(rowi) = coeff(2);
    fit_table.rsquared(rowi) = stats.Rsquared;
    rowi = rowi+1;
    points_plotted_5SD{regi}=[{Prob},{xbin}];
    plotted_fit_5SD{regi}=[{fitbinCenter},{regfit}];
end
ax = gca;
% ax.FontSize = 20;
ax.YScale = 'log';
ax.XLim = [-inf 2];
axis square
grid on
hold off
% legend(pl, well_regList)
% saveas(gcf, 'D:\Brewer lab data\HFS\figs_nostim\well-isi-z1.png')
xlabel('Interspike Interval (s)')
ylabel('Pr(X \geq x)')
xlim([1e-2 2e-1]); ylim([1e-1 1])
xticks([1e-2, 0.25e-1,0.5e-1, 1e-1, 2e-1])
yticks([1e-1, 0.25, 0.5, 1])
% ax.FontSize = 16;
% saveas(gcf, 'D:\Brewer lab data\HFS\figs_nostim\well-isi.png')
% save('well_ISI_points_plotted_nostim','points_plotted')
% save('well_ISI_fit_plotted_nostim','plotted_fit')
%% stats and ANCOVA
fit_table 
%ANOVA
[~,~,~,stats] = aoctool(sp_feature,prob_vec,regLabel,0.05,'','','','off');
[c,means_5SD]=multcompare(stats,0.05,'on','','s');

sp_feature_post_well=sp_feature;
prob_vec_post_well=prob_vec;
post_stim_string_vec=repmat("5SD",length(regLabel),1);
regLabel_string=string(regLabel);
regLabel_post_well=categorical(join([regLabel_string,post_stim_string_vec]));

%save('no_stim_well_ancova_vars','sp_feature_post_well','prob_vec_post_well','regLabel_post_well')
%save('well_means_nostim','means')
%% Error bar graph
figure( 'Position', [100 100 700 600])
xdata = 1:4;
ydata = means_5SD(:,1);
errdata = means_5SD(:,2);
b = bar(xdata,ydata);
b.FaceColor = 'flat';
for regi=1:4
    b.CData(regi,:) = col_array_rgb(regi,:);
end
hold on
errorbar(xdata, ydata, errdata, '.k')
xticks(1:4), xticklabels(well_regList) 
ylim([-1 0])
yticks(linspace(-1.5,0,7))
% yticks(log10(ytick_nums)); yticklabels(ytick_nums)
ylabel 'Slope of ISI (s^{-1})'
set(gca,'fontsize',17.5)
%saveas(gcf, 'D:\Brewer lab data\HFS\figs_nostim\well-isi-bar.png')
%% no stim 500 sd
load('D:\Brewer lab data\HFS\No Stim\Wells_5SD_500maxSD\well_spike_dynamics_table_nostim.mat')
well_spike_burst_dyn_table = well_spike_dynamics_table;
clear spike_burst_dyn_table
well_regList = {'EC','DG','CA3','CA1'}; 
figOrder = [1 2 4 3 ];
col_array = ["#D95319","#77AC30","#4DBEEE","#7E2F8E"];
col_array_rgb = [0.8500 0.3250 0.0980; 0.4660 0.6740 0.1880; ...
    0.3010 0.7450 0.9330; 0.4940 0.1840 0.5560];
powerlawfitfun = @(b,x) 10.^(b(2)*log10(x) + b(1)); 
%% find n
ec_n=vertcat(well_spike_burst_dyn_table(well_spike_burst_dyn_table.regi==1,:).ISI{:});
dg_n=vertcat(well_spike_burst_dyn_table(well_spike_burst_dyn_table.regi==2,:).ISI{:});
ca3_n=vertcat(well_spike_burst_dyn_table(well_spike_burst_dyn_table.regi==3,:).ISI{:});
ca1_n=vertcat(well_spike_burst_dyn_table(well_spike_burst_dyn_table.regi==4,:).ISI{:});
%% ISI
fit_table = table([],[],[],[],'variablenames',{'regi','slope','intercept','rsquared'}); rowi=1;
figure( 'Position', [100 100 700 600])
binEdge = logspace(-2,1,100); binCenter = convert_edges_2_centers(binEdge);
xlimits = [ 1e-2, 2e-1];
sp_feature=[]; prob_vec=[]; regLabel=[]; pl=[];
hold on
disp('CCDF for Interspike Interval for FF axons')
for regi=1:4
    vec2plot = cell2mat(well_spike_burst_dyn_table(well_spike_burst_dyn_table.regi==regi,:).ISI);
    [ hist_object, hist_plot, Prob, xbin] = log_binned_histogram( vec2plot, binEdge, "cdf" );  
    hist_plot.Color = col_array(regi);
    pl(regi)=hist_plot;
    
    % fitting data between [xmin, xmax]
    [ coeff, stats, opt_lim, logBinC,logHistProb ] = find_powerlawfit_with_grid_search(binCenter, ...
        Prob, xlimits);
    sp_feature = [sp_feature; logBinC']; prob_vec=[prob_vec; logHistProb']; 
    regLabel = [regLabel; repmat(categorical(well_regList(regi)),length(logBinC),1)];
%     fitbinCenter = binCenter( binCenter > opt_lim(1) & binCenter< opt_lim(end));
    fitbinCenter = binCenter;
    regfit = powerlawfitfun(coeff,fitbinCenter);
    plot(fitbinCenter , regfit ,'--k')
    
    fit_table.regi(rowi) = regi;
    fit_table.intercept(rowi) = coeff(1);
    fit_table.slope(rowi) = coeff(2);
    fit_table.rsquared(rowi) = stats.Rsquared;
    rowi = rowi+1;
    points_plotted_5SD_500max{regi}=[{Prob},{xbin}];
    plotted_fit_5SD_500max{regi}=[{fitbinCenter},{regfit}];
end
ax = gca;
% ax.FontSize = 20;
ax.YScale = 'log';
ax.XLim = [-inf 2];
axis square
grid on
hold off
% legend(pl, well_regList)
% saveas(gcf, 'D:\Brewer lab data\HFS\figs_nostim\well-isi-z1.png')
xlabel('Interspike Interval (s)')
ylabel('Pr(X \geq x)')
xlim([1e-2 2e-1]); ylim([1e-1 1])
xticks([1e-2, 0.25e-1,0.5e-1, 1e-1, 2e-1])
yticks([1e-1, 0.25, 0.5, 1])
% ax.FontSize = 16;
% saveas(gcf, 'D:\Brewer lab data\HFS\figs_nostim\well-isi.png')
% save('well_ISI_points_plotted_nostim','points_plotted')
% save('well_ISI_fit_plotted_nostim','plotted_fit')
%% stats and ANCOVA
fit_table 
%ANOVA
[~,~,~,stats] = aoctool(sp_feature,prob_vec,regLabel,0.05,'','','','off');
[c,means_5SD_500max]=multcompare(stats,0.05,'on','','s');

sp_feature_post_well=[sp_feature_post_well; sp_feature];
prob_vec_post_well=[prob_vec_post_well; prob_vec];
post_stim_string_vec=repmat("5SD-500SD",length(regLabel),1);
regLabel_string=string(regLabel);
regLabel_post_well=[regLabel_post_well; categorical(join([regLabel_string,post_stim_string_vec]))];

%save('no_stim_well_ancova_vars','sp_feature_post_well','prob_vec_post_well','regLabel_post_well')
%save('well_means_nostim','means')
%% Error bar graph
figure( 'Position', [100 100 700 600])
xdata = 1:4;
ydata = means_5SD_500max(:,1);
errdata = means_5SD_500max(:,2);
b = bar(xdata,ydata);
b.FaceColor = 'flat';
for regi=1:4
    b.CData(regi,:) = col_array_rgb(regi,:);
end
hold on
errorbar(xdata, ydata, errdata, '.k')
xticks(1:4), xticklabels(well_regList) 
ylim([-1 0])
yticks(linspace(-1.5,0,7))
% yticks(log10(ytick_nums)); yticklabels(ytick_nums)
ylabel 'Slope of ISI (s^{-1})'
set(gca,'fontsize',17.5)
%saveas(gcf, 'D:\Brewer lab data\HFS\figs_nostim\well-isi-bar.png')

%% Plotting together

subregions=["EC","DG","CA3","CA1"];
close all
matlab_colors=colororder;
for i=1:4
    to_plot={points_plotted_5SD{i},points_plotted_5SD_500max{i}};
    to_plot_fit={plotted_fit_5SD{i},plotted_fit_5SD_500max{i}};
    figure( 'Position', [100 100 700 600])
    hold on
    ax=gca;
    ax.XScale='log';
    for j=1:2
        plot(to_plot{j}{2},to_plot{j}{1},'LineWidth',2,'Color',matlab_colors(j,:))
        plot(to_plot_fit{j}{1},to_plot_fit{j}{2},'--k','LineWidth',2)
    end
    %ax = gca;
    % ax.FontSize = 20;
    ax.YScale = 'log';
    ax.XLim = [-inf 2];
    axis square
    grid on
    hold off
    % legend(pl, well_regList)
    % saveas(gcf, 'D:\Brewer lab data\HFS\figs_nostim\well-isi-z1.png')
    saveas(gcf,strcat('D:\Brewer lab data\HFS\well_SD_Comparison_pics_2\',subregions(i),'_plot_z1_nostim.png'))
    xlabel('Interspike Interval (s)')
    ylabel('Pr(X \geq x)')
    xlim([1e-2 2e-1]); ylim([1e-1 1])
    xticks([1e-2, 0.25e-1,0.5e-1, 1e-1, 2e-1])
    yticks([1e-1, 0.25, 0.5, 1])
    hold off
    set(gca,'fontsize',16)
    title(subregions(i),'FontSize',17)
    saveas(gcf,strcat('D:\Brewer lab data\HFS\well_SD_Comparison_pics_2\',subregions(i),'_plot_nostim.png'))
end
%% all stats and ancova
subregions=["EC","DG","CA3","CA1"];
close all
matlab_colors=colororder;
anove_struct=[];
perc_change=[];
for i=1:4
    figure( 'Position', [100 100 700 600])
    to_compare=(contains(string(regLabel_post_well),subregions(i)));
    [~,~,~,stats] = aoctool(sp_feature_post_well(to_compare),...
        prob_vec_post_well(to_compare),regLabel_post_well(to_compare),0.05,'','','','off');
    [c,means]=multcompare(stats,0.05,'on','','s');
    anova_struct{i}=c;
    xdata = 1:2;
    ydata = means(:,1);
    errdata = means(:,2);
    b = bar(xdata,ydata);
    b.FaceColor = 'flat';
    for regi=1:2
        b.CData(regi,:) = matlab_colors(regi,:);
    end
    hold on
    errorbar(xdata, ydata, errdata, '.k')
    xticks(1:4), xticklabels(unique(string(regLabel_post_well(to_compare)),'stable')) 
    ylim([-1.5 0])
    yticks(linspace(-1.5,0,7))
    % yticks(log10(ytick_nums)); yticklabels(ytick_nums)
    ylabel 'Slope of ISI (s^{-1})'
    title(subregions(i))
    set(gca,'fontsize',16)
    saveas(gcf,strcat('D:\Brewer lab data\HFS\well_SD_Comparison_pics_2\',subregions(i),'_bar_nostim.png'))
    combos=nchoosek([1:2],2);
    for j=1:length(combos(:,1))
        diff_vec(j,1)=((means(combos(j,2),1)-means(combos(j,1),1))/means(combos(j,1),1))*100;
    end
    perc_change{i}=[combos,diff_vec];
end

%% 5 HFS

load('D:\Brewer lab data\HFS\Theta Stim\Wells_5SD\well_spike_dynamics_table_theta.mat')
well_spike_burst_dyn_table = well_spike_dynamics_table;
clear spike_burst_dyn_table
well_regList = {'EC','DG','CA3','CA1'}; 
figOrder = [1 2 4 3 ];
col_array = ["#D95319","#77AC30","#4DBEEE","#7E2F8E"];
col_array_rgb = [0.8500 0.3250 0.0980; 0.4660 0.6740 0.1880; ...
    0.3010 0.7450 0.9330; 0.4940 0.1840 0.5560];
powerlawfitfun = @(b,x) 10.^(b(2)*log10(x) + b(1)); 

%% find n
ec_n=vertcat(well_spike_burst_dyn_table(well_spike_burst_dyn_table.regi==1,:).ISI{:});
dg_n=vertcat(well_spike_burst_dyn_table(well_spike_burst_dyn_table.regi==2,:).ISI{:});
ca3_n=vertcat(well_spike_burst_dyn_table(well_spike_burst_dyn_table.regi==3,:).ISI{:});
ca1_n=vertcat(well_spike_burst_dyn_table(well_spike_burst_dyn_table.regi==4,:).ISI{:});

%% ISI
fit_table = table([],[],[],[],'variablenames',{'regi','slope','intercept','rsquared'}); rowi=1;
figure( 'Position', [100 100 700 600])
binEdge = logspace(-2,1,100); binCenter = convert_edges_2_centers(binEdge);
xlimits = [ 1e-2, 2e-1];
sp_feature=[]; prob_vec=[]; regLabel=[]; pl=[];
hold on
disp('CCDF for Interspike Interval for FF axons')
for regi=1:4
    vec2plot = cell2mat(well_spike_burst_dyn_table(well_spike_burst_dyn_table.regi==regi,:).ISI);
    [ hist_object, hist_plot, Prob, xbin] = log_binned_histogram( vec2plot, binEdge, "cdf" );  
    hist_plot.Color = col_array(regi);
    pl(regi)=hist_plot;
    
    % fitting data between [xmin, xmax]
    [ coeff, stats, opt_lim, logBinC,logHistProb ] = find_powerlawfit_with_grid_search(binCenter, ...
        Prob, xlimits);
    sp_feature = [sp_feature; logBinC']; prob_vec=[prob_vec; logHistProb']; 
    regLabel = [regLabel; repmat(categorical(well_regList(regi)),length(logBinC),1)];
%     fitbinCenter = binCenter( binCenter > opt_lim(1) & binCenter< opt_lim(end));
    fitbinCenter = binCenter;
    regfit = powerlawfitfun(coeff,fitbinCenter);
    plot(fitbinCenter , regfit ,'--k')
    
    fit_table.regi(rowi) = regi;
    fit_table.intercept(rowi) = coeff(1);
    fit_table.slope(rowi) = coeff(2);
    fit_table.rsquared(rowi) = stats.Rsquared;
    rowi = rowi+1;
    points_plotted_5SD{regi}=[{Prob},{xbin}];
    plotted_fit_5SD{regi}=[{fitbinCenter},{regfit}];
end
ax = gca;
% ax.FontSize = 20;
ax.YScale = 'log';
ax.XLim = [-inf 2];
axis square
grid on
hold off
% legend(pl, well_regList)
% saveas(gcf, 'D:\Brewer lab data\HFS\figs_nostim\well-isi-z1.png')
xlabel('Interspike Interval (s)')
ylabel('Pr(X \geq x)')
xlim([1e-2 2e-1]); ylim([1e-1 1])
xticks([1e-2, 0.25e-1,0.5e-1, 1e-1, 2e-1])
yticks([1e-1, 0.25, 0.5, 1])
% ax.FontSize = 16;
% saveas(gcf, 'D:\Brewer lab data\HFS\figs_nostim\well-isi.png')
% save('well_ISI_points_plotted_nostim','points_plotted')
% save('well_ISI_fit_plotted_nostim','plotted_fit')
%% stats and ANCOVA
fit_table 
%ANOVA
[~,~,~,stats] = aoctool(sp_feature,prob_vec,regLabel,0.05,'','','','off');
[c,means_5SD]=multcompare(stats,0.05,'on','','s');

sp_feature_post_well=sp_feature;
prob_vec_post_well=prob_vec;
post_stim_string_vec=repmat("5SD",length(regLabel),1);
regLabel_string=string(regLabel);
regLabel_post_well=categorical(join([regLabel_string,post_stim_string_vec]));

%save('no_stim_well_ancova_vars','sp_feature_post_well','prob_vec_post_well','regLabel_post_well')
%save('well_means_nostim','means')
%% Error bar graph
figure( 'Position', [100 100 700 600])
xdata = 1:4;
ydata = means_5SD(:,1);
errdata = means_5SD(:,2);
b = bar(xdata,ydata);
b.FaceColor = 'flat';
for regi=1:4
    b.CData(regi,:) = col_array_rgb(regi,:);
end
hold on
errorbar(xdata, ydata, errdata, '.k')
xticks(1:4), xticklabels(well_regList) 
ylim([-1 0])
yticks(linspace(-1.5,0,7))
% yticks(log10(ytick_nums)); yticklabels(ytick_nums)
ylabel 'Slope of ISI (s^{-1})'
set(gca,'fontsize',17.5)
%saveas(gcf, 'D:\Brewer lab data\HFS\figs_nostim\well-isi-bar.png')
%% 5 HFS 500 SD
load('D:\Brewer lab data\HFS\Theta Stim\Wells_5SD_500maxSD\well_spike_dynamics_table_theta.mat')
well_spike_burst_dyn_table = well_spike_dynamics_table;
clear spike_burst_dyn_table
well_regList = {'EC','DG','CA3','CA1'}; 
figOrder = [1 2 4 3 ];
col_array = ["#D95319","#77AC30","#4DBEEE","#7E2F8E"];
col_array_rgb = [0.8500 0.3250 0.0980; 0.4660 0.6740 0.1880; ...
    0.3010 0.7450 0.9330; 0.4940 0.1840 0.5560];
powerlawfitfun = @(b,x) 10.^(b(2)*log10(x) + b(1)); 
%% find n
ec_n=vertcat(well_spike_burst_dyn_table(well_spike_burst_dyn_table.regi==1,:).ISI{:});
dg_n=vertcat(well_spike_burst_dyn_table(well_spike_burst_dyn_table.regi==2,:).ISI{:});
ca3_n=vertcat(well_spike_burst_dyn_table(well_spike_burst_dyn_table.regi==3,:).ISI{:});
ca1_n=vertcat(well_spike_burst_dyn_table(well_spike_burst_dyn_table.regi==4,:).ISI{:});
%% ISI
fit_table = table([],[],[],[],'variablenames',{'regi','slope','intercept','rsquared'}); rowi=1;
figure( 'Position', [100 100 700 600])
binEdge = logspace(-2,1,100); binCenter = convert_edges_2_centers(binEdge);
xlimits = [ 1e-2, 2e-1];
sp_feature=[]; prob_vec=[]; regLabel=[]; pl=[];
hold on
disp('CCDF for Interspike Interval for FF axons')
for regi=1:4
    vec2plot = cell2mat(well_spike_burst_dyn_table(well_spike_burst_dyn_table.regi==regi,:).ISI);
    [ hist_object, hist_plot, Prob, xbin] = log_binned_histogram( vec2plot, binEdge, "cdf" );  
    hist_plot.Color = col_array(regi);
    pl(regi)=hist_plot;
    
    % fitting data between [xmin, xmax]
    [ coeff, stats, opt_lim, logBinC,logHistProb ] = find_powerlawfit_with_grid_search(binCenter, ...
        Prob, xlimits);
    sp_feature = [sp_feature; logBinC']; prob_vec=[prob_vec; logHistProb']; 
    regLabel = [regLabel; repmat(categorical(well_regList(regi)),length(logBinC),1)];
%     fitbinCenter = binCenter( binCenter > opt_lim(1) & binCenter< opt_lim(end));
    fitbinCenter = binCenter;
    regfit = powerlawfitfun(coeff,fitbinCenter);
    plot(fitbinCenter , regfit ,'--k')
    
    fit_table.regi(rowi) = regi;
    fit_table.intercept(rowi) = coeff(1);
    fit_table.slope(rowi) = coeff(2);
    fit_table.rsquared(rowi) = stats.Rsquared;
    rowi = rowi+1;
    points_plotted_5SD_500max{regi}=[{Prob},{xbin}];
    plotted_fit_5SD_500max{regi}=[{fitbinCenter},{regfit}];
end
ax = gca;
% ax.FontSize = 20;
ax.YScale = 'log';
ax.XLim = [-inf 2];
axis square
grid on
hold off
% legend(pl, well_regList)
% saveas(gcf, 'D:\Brewer lab data\HFS\figs_nostim\well-isi-z1.png')
xlabel('Interspike Interval (s)')
ylabel('Pr(X \geq x)')
xlim([1e-2 2e-1]); ylim([1e-1 1])
xticks([1e-2, 0.25e-1,0.5e-1, 1e-1, 2e-1])
yticks([1e-1, 0.25, 0.5, 1])
% ax.FontSize = 16;
% saveas(gcf, 'D:\Brewer lab data\HFS\figs_nostim\well-isi.png')
% save('well_ISI_points_plotted_nostim','points_plotted')
% save('well_ISI_fit_plotted_nostim','plotted_fit')
%% stats and ANCOVA
fit_table 
%ANOVA
[~,~,~,stats] = aoctool(sp_feature,prob_vec,regLabel,0.05,'','','','off');
[c,means_5SD_500max]=multcompare(stats,0.05,'on','','s');

sp_feature_post_well=[sp_feature_post_well; sp_feature];
prob_vec_post_well=[prob_vec_post_well; prob_vec];
post_stim_string_vec=repmat("5SD-500SD",length(regLabel),1);
regLabel_string=string(regLabel);
regLabel_post_well=[regLabel_post_well; categorical(join([regLabel_string,post_stim_string_vec]))];

%save('no_stim_well_ancova_vars','sp_feature_post_well','prob_vec_post_well','regLabel_post_well')
%save('well_means_nostim','means')
%% Error bar graph
figure( 'Position', [100 100 700 600])
xdata = 1:4;
ydata = means_5SD_500max(:,1);
errdata = means_5SD_500max(:,2);
b = bar(xdata,ydata);
b.FaceColor = 'flat';
for regi=1:4
    b.CData(regi,:) = col_array_rgb(regi,:);
end
hold on
errorbar(xdata, ydata, errdata, '.k')
xticks(1:4), xticklabels(well_regList) 
ylim([-1 0])
yticks(linspace(-1.5,0,7))
% yticks(log10(ytick_nums)); yticklabels(ytick_nums)
ylabel 'Slope of ISI (s^{-1})'
set(gca,'fontsize',17.5)
%saveas(gcf, 'D:\Brewer lab data\HFS\figs_nostim\well-isi-bar.png')

%% Plotting together

subregions=["EC","DG","CA3","CA1"];
close all
matlab_colors=colororder;
for i=1:4
    to_plot={points_plotted_5SD{i},points_plotted_5SD_500max{i}};
    to_plot_fit={plotted_fit_5SD{i},plotted_fit_5SD_500max{i}};
    figure( 'Position', [100 100 700 600])
    hold on
    ax=gca;
    ax.XScale='log';
    for j=1:2
        plot(to_plot{j}{2},to_plot{j}{1},'LineWidth',2,'Color',matlab_colors(j,:))
        plot(to_plot_fit{j}{1},to_plot_fit{j}{2},'--k','LineWidth',2)
    end
    %ax = gca;
    % ax.FontSize = 20;
    ax.YScale = 'log';
    ax.XLim = [-inf 2];
    axis square
    grid on
    hold off
    % legend(pl, well_regList)
    % saveas(gcf, 'D:\Brewer lab data\HFS\figs_nostim\well-isi-z1.png')
    saveas(gcf,strcat('D:\Brewer lab data\HFS\well_SD_Comparison_pics_2\',subregions(i),'_plot_z1_5hfs.png'))
    xlabel('Interspike Interval (s)')
    ylabel('Pr(X \geq x)')
    xlim([1e-2 2e-1]); ylim([1e-1 1])
    xticks([1e-2, 0.25e-1,0.5e-1, 1e-1, 2e-1])
    yticks([1e-1, 0.25, 0.5, 1])
    hold off
    set(gca,'fontsize',16)
    title(subregions(i),'FontSize',17)
    saveas(gcf,strcat('D:\Brewer lab data\HFS\well_SD_Comparison_pics_2\',subregions(i),'_plot_5hfs.png'))
end
%% all stats and ancova
subregions=["EC","DG","CA3","CA1"];
close all
matlab_colors=colororder;
anove_struct=[];
perc_change=[];
for i=1:4
    figure( 'Position', [100 100 700 600])
    to_compare=(contains(string(regLabel_post_well),subregions(i)));
    [~,~,~,stats] = aoctool(sp_feature_post_well(to_compare),...
        prob_vec_post_well(to_compare),regLabel_post_well(to_compare),0.05,'','','','off');
    [c,means]=multcompare(stats,0.05,'on','','s');
    anova_struct{i}=c;
    xdata = 1:2;
    ydata = means(:,1);
    errdata = means(:,2);
    b = bar(xdata,ydata);
    b.FaceColor = 'flat';
    for regi=1:2
        b.CData(regi,:) = matlab_colors(regi,:);
    end
    hold on
    errorbar(xdata, ydata, errdata, '.k')
    xticks(1:4), xticklabels(unique(string(regLabel_post_well(to_compare)),'stable')) 
    ylim([-1.5 0])
    yticks(linspace(-1.5,0,7))
    % yticks(log10(ytick_nums)); yticklabels(ytick_nums)
    ylabel 'Slope of ISI (s^{-1})'
    title(subregions(i))
    set(gca,'fontsize',16)
    saveas(gcf,strcat('D:\Brewer lab data\HFS\well_SD_Comparison_pics_2\',subregions(i),'_bar_5hfs.png'))
    combos=nchoosek([1:2],2);
    for j=1:length(combos(:,1))
        diff_vec(j,1)=((means(combos(j,2),1)-means(combos(j,1),1))/means(combos(j,1),1))*100;
    end
    perc_change{i}=[combos,diff_vec];
end

%% 40 HFS

load('D:\Brewer lab data\HFS\HFS Stim\Wells_5SD\well_spike_dynamics_table_HFS.mat')
well_spike_burst_dyn_table = well_spike_dynamics_table;
clear spike_burst_dyn_table
well_regList = {'EC','DG','CA3','CA1'}; 
figOrder = [1 2 4 3 ];
col_array = ["#D95319","#77AC30","#4DBEEE","#7E2F8E"];
col_array_rgb = [0.8500 0.3250 0.0980; 0.4660 0.6740 0.1880; ...
    0.3010 0.7450 0.9330; 0.4940 0.1840 0.5560];
powerlawfitfun = @(b,x) 10.^(b(2)*log10(x) + b(1)); 

%% find n
ec_n=vertcat(well_spike_burst_dyn_table(well_spike_burst_dyn_table.regi==1,:).ISI{:});
dg_n=vertcat(well_spike_burst_dyn_table(well_spike_burst_dyn_table.regi==2,:).ISI{:});
ca3_n=vertcat(well_spike_burst_dyn_table(well_spike_burst_dyn_table.regi==3,:).ISI{:});
ca1_n=vertcat(well_spike_burst_dyn_table(well_spike_burst_dyn_table.regi==4,:).ISI{:});

%% ISI
fit_table = table([],[],[],[],'variablenames',{'regi','slope','intercept','rsquared'}); rowi=1;
figure( 'Position', [100 100 700 600])
binEdge = logspace(-2,1,100); binCenter = convert_edges_2_centers(binEdge);
xlimits = [ 1e-2, 2e-1];
sp_feature=[]; prob_vec=[]; regLabel=[]; pl=[];
hold on
disp('CCDF for Interspike Interval for FF axons')
for regi=1:4
    vec2plot = cell2mat(well_spike_burst_dyn_table(well_spike_burst_dyn_table.regi==regi,:).ISI);
    [ hist_object, hist_plot, Prob, xbin] = log_binned_histogram( vec2plot, binEdge, "cdf" );  
    hist_plot.Color = col_array(regi);
    pl(regi)=hist_plot;
    
    % fitting data between [xmin, xmax]
    [ coeff, stats, opt_lim, logBinC,logHistProb ] = find_powerlawfit_with_grid_search(binCenter, ...
        Prob, xlimits);
    sp_feature = [sp_feature; logBinC']; prob_vec=[prob_vec; logHistProb']; 
    regLabel = [regLabel; repmat(categorical(well_regList(regi)),length(logBinC),1)];
%     fitbinCenter = binCenter( binCenter > opt_lim(1) & binCenter< opt_lim(end));
    fitbinCenter = binCenter;
    regfit = powerlawfitfun(coeff,fitbinCenter);
    plot(fitbinCenter , regfit ,'--k')
    
    fit_table.regi(rowi) = regi;
    fit_table.intercept(rowi) = coeff(1);
    fit_table.slope(rowi) = coeff(2);
    fit_table.rsquared(rowi) = stats.Rsquared;
    rowi = rowi+1;
    points_plotted_5SD{regi}=[{Prob},{xbin}];
    plotted_fit_5SD{regi}=[{fitbinCenter},{regfit}];
end
ax = gca;
% ax.FontSize = 20;
ax.YScale = 'log';
ax.XLim = [-inf 2];
axis square
grid on
hold off
% legend(pl, well_regList)
% saveas(gcf, 'D:\Brewer lab data\HFS\figs_nostim\well-isi-z1.png')
xlabel('Interspike Interval (s)')
ylabel('Pr(X \geq x)')
xlim([1e-2 2e-1]); ylim([1e-1 1])
xticks([1e-2, 0.25e-1,0.5e-1, 1e-1, 2e-1])
yticks([1e-1, 0.25, 0.5, 1])
% ax.FontSize = 16;
% saveas(gcf, 'D:\Brewer lab data\HFS\figs_nostim\well-isi.png')
% save('well_ISI_points_plotted_nostim','points_plotted')
% save('well_ISI_fit_plotted_nostim','plotted_fit')
%% stats and ANCOVA
fit_table 
%ANOVA
[~,~,~,stats] = aoctool(sp_feature,prob_vec,regLabel,0.05,'','','','off');
[c,means_5SD]=multcompare(stats,0.05,'on','','s');

sp_feature_post_well=sp_feature;
prob_vec_post_well=prob_vec;
post_stim_string_vec=repmat("5SD",length(regLabel),1);
regLabel_string=string(regLabel);
regLabel_post_well=categorical(join([regLabel_string,post_stim_string_vec]));

%save('no_stim_well_ancova_vars','sp_feature_post_well','prob_vec_post_well','regLabel_post_well')
%save('well_means_nostim','means')
%% Error bar graph
figure( 'Position', [100 100 700 600])
xdata = 1:4;
ydata = means_5SD(:,1);
errdata = means_5SD(:,2);
b = bar(xdata,ydata);
b.FaceColor = 'flat';
for regi=1:4
    b.CData(regi,:) = col_array_rgb(regi,:);
end
hold on
errorbar(xdata, ydata, errdata, '.k')
xticks(1:4), xticklabels(well_regList) 
ylim([-1 0])
yticks(linspace(-1.5,0,7))
% yticks(log10(ytick_nums)); yticklabels(ytick_nums)
ylabel 'Slope of ISI (s^{-1})'
set(gca,'fontsize',17.5)
%saveas(gcf, 'D:\Brewer lab data\HFS\figs_nostim\well-isi-bar.png')
%% 40 HFS 500 SD
load('D:\Brewer lab data\HFS\HFS Stim\Wells_5SD_500maxSD\well_spike_dynamics_table_hfs.mat')
well_spike_burst_dyn_table = well_spike_dynamics_table;
clear spike_burst_dyn_table
well_regList = {'EC','DG','CA3','CA1'}; 
figOrder = [1 2 4 3 ];
col_array = ["#D95319","#77AC30","#4DBEEE","#7E2F8E"];
col_array_rgb = [0.8500 0.3250 0.0980; 0.4660 0.6740 0.1880; ...
    0.3010 0.7450 0.9330; 0.4940 0.1840 0.5560];
powerlawfitfun = @(b,x) 10.^(b(2)*log10(x) + b(1)); 
%% find n
ec_n=vertcat(well_spike_burst_dyn_table(well_spike_burst_dyn_table.regi==1,:).ISI{:});
dg_n=vertcat(well_spike_burst_dyn_table(well_spike_burst_dyn_table.regi==2,:).ISI{:});
ca3_n=vertcat(well_spike_burst_dyn_table(well_spike_burst_dyn_table.regi==3,:).ISI{:});
ca1_n=vertcat(well_spike_burst_dyn_table(well_spike_burst_dyn_table.regi==4,:).ISI{:});
%% ISI
fit_table = table([],[],[],[],'variablenames',{'regi','slope','intercept','rsquared'}); rowi=1;
figure( 'Position', [100 100 700 600])
binEdge = logspace(-2,1,100); binCenter = convert_edges_2_centers(binEdge);
xlimits = [ 1e-2, 2e-1];
sp_feature=[]; prob_vec=[]; regLabel=[]; pl=[];
hold on
disp('CCDF for Interspike Interval for FF axons')
for regi=1:4
    vec2plot = cell2mat(well_spike_burst_dyn_table(well_spike_burst_dyn_table.regi==regi,:).ISI);
    [ hist_object, hist_plot, Prob, xbin] = log_binned_histogram( vec2plot, binEdge, "cdf" );  
    hist_plot.Color = col_array(regi);
    pl(regi)=hist_plot;
    
    % fitting data between [xmin, xmax]
    [ coeff, stats, opt_lim, logBinC,logHistProb ] = find_powerlawfit_with_grid_search(binCenter, ...
        Prob, xlimits);
    sp_feature = [sp_feature; logBinC']; prob_vec=[prob_vec; logHistProb']; 
    regLabel = [regLabel; repmat(categorical(well_regList(regi)),length(logBinC),1)];
%     fitbinCenter = binCenter( binCenter > opt_lim(1) & binCenter< opt_lim(end));
    fitbinCenter = binCenter;
    regfit = powerlawfitfun(coeff,fitbinCenter);
    plot(fitbinCenter , regfit ,'--k')
    
    fit_table.regi(rowi) = regi;
    fit_table.intercept(rowi) = coeff(1);
    fit_table.slope(rowi) = coeff(2);
    fit_table.rsquared(rowi) = stats.Rsquared;
    rowi = rowi+1;
    points_plotted_5SD_500max{regi}=[{Prob},{xbin}];
    plotted_fit_5SD_500max{regi}=[{fitbinCenter},{regfit}];
end
ax = gca;
% ax.FontSize = 20;
ax.YScale = 'log';
ax.XLim = [-inf 2];
axis square
grid on
hold off
% legend(pl, well_regList)
% saveas(gcf, 'D:\Brewer lab data\HFS\figs_nostim\well-isi-z1.png')
xlabel('Interspike Interval (s)')
ylabel('Pr(X \geq x)')
xlim([1e-2 2e-1]); ylim([1e-1 1])
xticks([1e-2, 0.25e-1,0.5e-1, 1e-1, 2e-1])
yticks([1e-1, 0.25, 0.5, 1])
% ax.FontSize = 16;
% saveas(gcf, 'D:\Brewer lab data\HFS\figs_nostim\well-isi.png')
% save('well_ISI_points_plotted_nostim','points_plotted')
% save('well_ISI_fit_plotted_nostim','plotted_fit')
%% stats and ANCOVA
fit_table 
%ANOVA
[~,~,~,stats] = aoctool(sp_feature,prob_vec,regLabel,0.05,'','','','off');
[c,means_5SD_500max]=multcompare(stats,0.05,'on','','s');

sp_feature_post_well=[sp_feature_post_well; sp_feature];
prob_vec_post_well=[prob_vec_post_well; prob_vec];
post_stim_string_vec=repmat("5SD-500SD",length(regLabel),1);
regLabel_string=string(regLabel);
regLabel_post_well=[regLabel_post_well; categorical(join([regLabel_string,post_stim_string_vec]))];

%save('no_stim_well_ancova_vars','sp_feature_post_well','prob_vec_post_well','regLabel_post_well')
%save('well_means_nostim','means')
%% Error bar graph
figure( 'Position', [100 100 700 600])
xdata = 1:4;
ydata = means_5SD_500max(:,1);
errdata = means_5SD_500max(:,2);
b = bar(xdata,ydata);
b.FaceColor = 'flat';
for regi=1:4
    b.CData(regi,:) = col_array_rgb(regi,:);
end
hold on
errorbar(xdata, ydata, errdata, '.k')
xticks(1:4), xticklabels(well_regList) 
ylim([-1 0])
yticks(linspace(-1.5,0,7))
% yticks(log10(ytick_nums)); yticklabels(ytick_nums)
ylabel 'Slope of ISI (s^{-1})'
set(gca,'fontsize',17.5)
%saveas(gcf, 'D:\Brewer lab data\HFS\figs_nostim\well-isi-bar.png')

%% Plotting together

subregions=["EC","DG","CA3","CA1"];
close all
matlab_colors=colororder;
for i=1:4
    to_plot={points_plotted_5SD{i},points_plotted_5SD_500max{i}};
    to_plot_fit={plotted_fit_5SD{i},plotted_fit_5SD_500max{i}};
    figure( 'Position', [100 100 700 600])
    hold on
    ax=gca;
    ax.XScale='log';
    for j=1:2
        plot(to_plot{j}{2},to_plot{j}{1},'LineWidth',2,'Color',matlab_colors(j,:))
        plot(to_plot_fit{j}{1},to_plot_fit{j}{2},'--k','LineWidth',2)
    end
    %ax = gca;
    % ax.FontSize = 20;
    ax.YScale = 'log';
    ax.XLim = [-inf 2];
    axis square
    grid on
    hold off
    % legend(pl, well_regList)
    % saveas(gcf, 'D:\Brewer lab data\HFS\figs_nostim\well-isi-z1.png')
    saveas(gcf,strcat('D:\Brewer lab data\HFS\well_SD_Comparison_pics_2\',subregions(i),'_plot_z1_40hfs.png'))
    xlabel('Interspike Interval (s)')
    ylabel('Pr(X \geq x)')
    xlim([1e-2 2e-1]); ylim([1e-1 1])
    xticks([1e-2, 0.25e-1,0.5e-1, 1e-1, 2e-1])
    yticks([1e-1, 0.25, 0.5, 1])
    hold off
    set(gca,'fontsize',16)
    title(subregions(i),'FontSize',17)
    saveas(gcf,strcat('D:\Brewer lab data\HFS\well_SD_Comparison_pics_2\',subregions(i),'_plot_40hfs.png'))
end
%% all stats and ancova
subregions=["EC","DG","CA3","CA1"];
close all
matlab_colors=colororder;
anove_struct=[];
perc_change=[];
for i=1:4
    figure( 'Position', [100 100 700 600])
    to_compare=(contains(string(regLabel_post_well),subregions(i)));
    [~,~,~,stats] = aoctool(sp_feature_post_well(to_compare),...
        prob_vec_post_well(to_compare),regLabel_post_well(to_compare),0.05,'','','','off');
    [c,means]=multcompare(stats,0.05,'on','','s');
    anova_struct{i}=c;
    xdata = 1:2;
    ydata = means(:,1);
    errdata = means(:,2);
    b = bar(xdata,ydata);
    b.FaceColor = 'flat';
    for regi=1:2
        b.CData(regi,:) = matlab_colors(regi,:);
    end
    hold on
    errorbar(xdata, ydata, errdata, '.k')
    xticks(1:4), xticklabels(unique(string(regLabel_post_well(to_compare)),'stable')) 
    ylim([-1.5 0])
    yticks(linspace(-1.5,0,7))
    % yticks(log10(ytick_nums)); yticklabels(ytick_nums)
    ylabel 'Slope of ISI (s^{-1})'
    title(subregions(i))
    set(gca,'fontsize',16)
    saveas(gcf,strcat('D:\Brewer lab data\HFS\well_SD_Comparison_pics_2\',subregions(i),'_bar_40hfs.png'))
    combos=nchoosek([1:2],2);
    for j=1:length(combos(:,1))
        diff_vec(j,1)=((means(combos(j,2),1)-means(combos(j,1),1))/means(combos(j,1),1))*100;
    end
    perc_change{i}=[combos,diff_vec];
end
