clear
clc
%load('D:\Brewer lab data\HFS\No Stim\Wells_5SD\well_spike_dynamics_table_nostim.mat')
load('D:\Brewer lab data\HFS\No Stim\Wells_5SD_500maxSD\well_spike_dynamics_table_nostim.mat')
well_spike_burst_dyn_table = well_spike_dynamics_table;
clear spike_burst_dyn_table
well_regList = {'EC','DG','CA3','CA1'}; 
figOrder = [1 2 4 3 ];
col_array = ["#D95319","#77AC30","#4DBEEE","#7E2F8E"];
col_array_rgb = [0.8500 0.3250 0.0980; 0.4660 0.6740 0.1880; ...
    0.3010 0.7450 0.9330; 0.4940 0.1840 0.5560];
powerlawfitfun = @(b,x) 10.^(b(2)*log10(x) + b(1)); 
cd("D:\Brewer lab data\HFS\")
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
    points_plotted{regi}=[{Prob},{xbin}];
    plotted_fit{regi}=[{fitbinCenter},{regfit}];
end
ax = gca;
ax.FontSize = 20;
ax.YScale = 'log';
ax.XLim = [-inf 2];
axis square
grid on
hold off
% legend(pl, well_regList)
%saveas(gcf, 'D:\Brewer lab data\HFS\figs_nostim\well_5SD-isi-z1.png')
% saveas(gcf, 'D:\Brewer lab data\HFS\figs_nostim\well_5SD_500maxSD-isi-z1.png')
xlabel('Interspike Interval (s)')
ylabel('Pr(X \geq x)')
xlim([1e-2 2e-1]); ylim([1e-1 1])
xticks([1e-2, 0.25e-1,0.5e-1, 1e-1, 2e-1])
yticks([1e-1, 0.25, 0.5, 1])
% ax.FontSize = 16;
% saveas(gcf, 'D:\Brewer lab data\HFS\figs_nostim\well_5SD-isi.png')
% save('well_ISI_points_plotted_nostim_5SD','points_plotted')
% save('well_ISI_fit_plotted_nostim_5SD','plotted_fit')

saveas(gcf, 'D:\Brewer lab data\HFS\figs_nostim\well_5SD_500maxSD-isi.png')
save('well_ISI_points_plotted_nostim_5SD_500maxSD','points_plotted')
save('well_ISI_fit_plotted_nostim_5SD_500maxSD','plotted_fit')
%% stats and ANCOVA
fit_table 
%ANOVA
[~,~,~,stats] = aoctool(sp_feature,prob_vec,regLabel,0.05,'','','','off');
[c,means]=multcompare(stats,0.05,'on','','s');

sp_feature_post_well=sp_feature;
prob_vec_post_well=prob_vec;
post_stim_string_vec=repmat("No Stim",length(regLabel),1);
regLabel_string=string(regLabel);
regLabel_post_well=categorical(join([regLabel_string,post_stim_string_vec]));

% save('no_stim_well_ancova_vars_5SD','sp_feature_post_well','prob_vec_post_well','regLabel_post_well')
% save('well_means_nostim_5SD','means')

save('no_stim_well_ancova_vars_5SD_500maxSD','sp_feature_post_well','prob_vec_post_well','regLabel_post_well')
save('well_means_nostim_5SD_500maxSD','means')
%% Error bar graph
figure( 'Position', [100 100 700 600])
xdata = 1:4;
ydata = means(:,1);
errdata = means(:,2);
b = bar(xdata,ydata);
b.FaceColor = 'flat';
for regi=1:4
    b.CData(regi,:) = col_array_rgb(regi,:);
end
hold on
errorbar(xdata, ydata, errdata, '.k')
xticks(1:4), xticklabels(well_regList) 
ylim([-1 0])
yticks([-1.0 -0.9 -0.8, -0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0])
% yticks(log10(ytick_nums)); yticklabels(ytick_nums)
ylabel 'Slope of ISI (s^{-1})'
set(gca,'fontsize',17.5)
set(gca,'linewidth',2)
% saveas(gcf, 'D:\Brewer lab data\HFS\figs_nostim\well_5SD-isi-bar.png')

saveas(gcf, 'D:\Brewer lab data\HFS\figs_nostim\well_5SD_500maxSD-isi-bar.png')

%% IBI
fit_table = table([],{},{},{},{},'variablenames',{'regi','slope','intercept','rsquared','limits'}); rowi=1;
figure( 'Position', [100 100 700 600])
low_limit = 1.1e0;
xlimits = [1e-1, 6; 1e-1, 5; 1e-1 1.1; 1e-1 1.5]; x_limit_final = [3e1,2e1,1e1,1.2e1];
x_poi = [4.91389,4.91389,0.86881,1.7793];
sp_feature_1=[]; prob_vec_1=[]; regLabel_1=[]; pl=[];
sp_feature_2=[]; prob_vec_2=[]; regLabel_2=[];
hold on
binEdge = logspace(-1,2.5,120); binCenter = convert_edges_2_centers(binEdge);
for regi=1:4
    vec2plot = cell2mat(well_spike_burst_dyn_table( well_spike_burst_dyn_table.regi==regi,:).IBI);
    [ hist_object, hist_plot, Prob, xbin] = log_binned_histogram( vec2plot, binEdge, "cdf" );  
    hist_plot.Color = col_array(regi);
    pl(regi)=hist_plot;
    
    % fitting data between [xmin_lower, xmax_lower] for lower fit
    [ low_coeff, low_stats, opt_lim, logBinC,logHistProb ] = find_powerlawfit_with_grid_search(binCenter, ...
        Prob, xlimits(regi,:));
    sp_feature_1 = [sp_feature_1; logBinC']; prob_vec_1=[prob_vec_1; logHistProb']; 
    regLabel_1 = [regLabel_1; repmat(categorical(well_regList(regi)),length(logBinC),1)];
    fitbinCenter = binCenter(binCenter > opt_lim(1) & binCenter < x_poi(regi));
    regfit = powerlawfitfun(low_coeff,fitbinCenter);
    plot(fitbinCenter , regfit ,'--k','LineWidth',3)
    
    points_plotted{regi}=[{xbin},{Prob}];
    
    plotted_fit_1{regi}=[{fitbinCenter},{regfit}];
    if regi>0
        % fitting data between [x_max_lower, abs_xmax]
        [ high_coeff, high_stats, opt_lim_high, logBinC,logHistProb ] = find_powerlawfit_with_grid_search(binCenter, ...
            Prob, [max(opt_lim), x_limit_final(regi)],[0.5 0.5]  );
        sp_feature_2 = [sp_feature_2; logBinC']; prob_vec_2=[prob_vec_2; logHistProb']; 
        regLabel_2 = [regLabel_2; repmat(categorical(well_regList(regi)),length(logBinC),1)];
        fitbinCenter = binCenter(binCenter > x_poi(regi) & binCenter < opt_lim_high(2));
        regfit = powerlawfitfun(high_coeff,fitbinCenter);
        plot(fitbinCenter , regfit ,':','LineWidth',3,'Color',[0.5 0.5 0.5])
    else
        high_coeff = [0,0];
        high_stats.Rsquared = 0;
    end
    plotted_fit_2{regi}=[{fitbinCenter},{regfit}];
%     
    fit_table.regi(rowi) = regi;
    fit_table.intercept{rowi} = [low_coeff(1), high_coeff(1)];
    fit_table.slope{rowi} = [low_coeff(2), high_coeff(2)];
    fit_table.rsquared{rowi} = [low_stats.Rsquared, high_stats.Rsquared];
    fit_table.limits{rowi} = [opt_lim(2),opt_lim_high(2)];
    rowi = rowi+1;
end
% full scale
ax = gca;
ax.FontSize = 20;
ax.YScale = 'log';
ax.XLim = [-inf 2.5e1]; ax.YLim = [1e-2, 1];
xticks([0.1,0.25,1,2.5,5,10,20])
yticks([1e-2, 1e-1, 0.25, 0.5, 1])%cdf
axis square
hold off
% legend(pl, well_regList)
grid on 
xlabel('Interburst Interval (s)')
ylabel('Pr(X \geq x)')
% saveas(gcf, 'D:\Brewer lab data\HFS\figs_nostim\well_5SD-ibi.png')
saveas(gcf, 'D:\Brewer lab data\HFS\figs_nostim\well_5SD_500maxSD-ibi.png')
% % lower fits
% ax.XLim = [-inf 1e0]; ax.YLim = [0.5, 1];
% xticks([0.1, 0.25 0.5 1])
% yticks(0.5:0.1:1)
% saveas(gcf, './figs/well-ibi-z1.png')
% % upper fits
% ax.XLim = [1e0 1.1e1]; ax.YLim = [1e-2, 1e-1];
% xticks([0.1, 0.25 0.5 1])
% yticks([1e-2, 0.25e-1, 0.5e-1, 1e-1]);
% xticks([1e0, 0.25e1, 0.5e1, 1e1])
% saveas(gcf, './figs/well-ibi-z2.png')
fit_table.x_poi = x_poi'

% save('well_IBI_points_plotted_nostim_5SD','points_plotted')
% save('well_IBI_fit_plotted_nostim_1_5SD','plotted_fit_1')
% save('well_IBI_fit_plotted_nostim_2_5SD','plotted_fit_2')

save('well_IBI_points_plotted_nostim_5SD_500maxSD','points_plotted')
save('well_IBI_fit_plotted_nostim_1_5SD_500maxSD','plotted_fit_1')
save('well_IBI_fit_plotted_nostim_2_5SD_500maxSD','plotted_fit_2')
%% low ANOVA

[~,~,~,stats] = aoctool(sp_feature_1,prob_vec_1,regLabel_1,0.05,'','','','off');
[c,means_lower]=multcompare(stats,0.05,'on','','s');

% save("IBI_nostim1_well_means_5SD","means_lower")
save("IBI_nostim1_well_means_5SD_500maxSD","means_lower")

sp_feature_post_well_1=sp_feature_1;
prob_vec_post_well_1=prob_vec_1;
pre_stim_string_vec=repmat("No Stim",length(regLabel_1),1);
regLabel_string=string(regLabel_1);
regLabel_post_well_1=categorical(join([regLabel_string,pre_stim_string_vec]));
save('IBI_no_stim_well_ancova_vars_1_5SD','sp_feature_post_well_1','prob_vec_post_well_1','regLabel_post_well_1')
save('IBI_no_stim_well_ancova_vars_1_5SD_500maxSD','sp_feature_post_well_1','prob_vec_post_well_1','regLabel_post_well_1')

%% high ANOVA

[~,~,~,stats] = aoctool(sp_feature_2,prob_vec_2,regLabel_2,0.05,'','','','off');
[c,means_higher]=multcompare(stats,0.05,'on','','s');

% save("IBI_nostim2_well_means_5SD","means_higher")
save("IBI_nostim2_well_means_5SD_500maxSD","means_higher")


sp_feature_post_well_2=sp_feature_2;
prob_vec_post_well_2=prob_vec_2;
pre_stim_string_vec=repmat("No Stim",length(regLabel_2),1);
regLabel_string=string(regLabel_2);
regLabel_post_well_2=categorical(join([regLabel_string,pre_stim_string_vec]));
% save('IBI_no_stim_well_ancova_vars_2_5SD','sp_feature_post_well_2','prob_vec_post_well_2','regLabel_post_well_2')
save('IBI_no_stim_well_ancova_vars_2_5SD_500maxSD','sp_feature_post_well_2','prob_vec_post_well_2','regLabel_post_well_2')

%% lower Error bar graph
figure( 'Position', [100 100 700 600])
xdata = 1:4;
ydata = means_lower(:,1);
errdata = means_lower(:,2);
b = bar(xdata,ydata);
b.FaceColor = 'flat';
for regi=1:4
    b.CData(regi,:) = col_array_rgb(regi,:);
end
hold on
errorbar(xdata, ydata, errdata, '.k')
xticks(1:4), xticklabels(well_regList) 
ylim([-1 0])
ytick_nums = [0.1, 0.25, 0.5, 1];
% yticks(log10(ytick_nums)); yticklabels(ytick_nums)
ylabel 'Slope of IBI (s^{-1})'
set(gca,'fontsize',17.5)
% saveas(gcf, 'D:\Brewer lab data\HFS\figs_nostim\well_5SD-ibi-bar.png')
saveas(gcf, 'D:\Brewer lab data\HFS\figs_nostim\well_5SD_500maxSD-ibi-bar.png')

%% higher Error bar graph
figure( 'Position', [100 100 700 600])
xdata = 1:4;
ydata = means_higher(:,1);
errdata = means_higher(:,2);
width=0.35;
b = bar(xdata,ydata);
b.FaceColor = 'flat';
for regi=1:4
    b.CData(regi,:) = col_array_rgb(regi,:);
end
hold on
errorbar(xdata, ydata, errdata, '.k')
xticks(1:4), xticklabels(well_regList) 
ytick_nums = [0.005, 0.01,0.025, 0.05, 0.1, 0.25, 0.5, 1];
% yticks(log10(ytick_nums)); yticklabels(ytick_nums)
ylim([-2.5 0])
ylabel 'Slope of IBI (s^{-1})'
set(gca,'fontsize',17.5)
% saveas(gcf, 'D:\Brewer lab data\HFS\figs_nostim\well_5SD-ibi-bar-2.png')
saveas(gcf, 'D:\Brewer lab data\HFS\figs_nostim\well_5SD_500maxSD-ibi-bar-2.png')

%% SPNB
xstart = [4,4,4,4]; xstop = [1e2,1e2,2e1,1e2]
fit_table = table([],[],[],[],'variablenames',{'regi','slope','intercept','rsquared'}); rowi=1;
figure( 'Position', [100 100 700 600])
sp_feature=[]; prob_vec=[]; regLabel=[]; pl=[];
hold on
binEdge = logspace(log10(4),3,50); binCenter = convert_edges_2_centers(binEdge);
for regi=1:4
    vec2plot = cell2mat(well_spike_burst_dyn_table(well_spike_burst_dyn_table.regi==regi,:).SpikeperBurst);
    [ hist_object, hist_plot, Prob, xbin] = log_binned_histogram( vec2plot, binEdge, "cdf" );  
    hist_plot.Color = col_array(regi);
    pl(regi)=hist_plot;
    
    % fitting data between [xmin, xmax]
    xlimits = [xstart(regi), xstop(regi)];
    [ coeff, stats, opt_lim, logBinC,logHistProb ] = find_powerlawfit_with_grid_search(binCenter, ...
        Prob, xlimits);
    sp_feature = [sp_feature; logBinC']; prob_vec=[prob_vec; logHistProb']; 
    regLabel = [regLabel; repmat(categorical(well_regList(regi)),length(logBinC),1)];
%     fitbinCenter = binCenter( binCenter > opt_lim(1) & binCenter< opt_lim(end));
    fitbinCenter = binCenter;
    regfit = powerlawfitfun(coeff,fitbinCenter);
    plot(fitbinCenter , regfit ,'--k')
    
    points_plotted{regi}=[{Prob},{xbin}];
    plotted_fit{regi}=[{fitbinCenter},{regfit}];
    
    fit_table.regi(rowi) = regi;
    fit_table.intercept(rowi) = coeff(1);
    fit_table.slope(rowi) = coeff(2);
    fit_table.rsquared(rowi) = stats.Rsquared;
    rowi = rowi+1;
end
ax = gca;
% ax.FontSize = 20;
ax.YScale = 'log';
% ax.LineWidth = 2*ax.LineWidth;
% ax.TickLength= 2*ax.TickLength;
ax.YLim = [1e-1 1]; ax.XLim = [-inf 50];
xticks([6 12 25 50])
yticks([0.1 0.25 0.5 1])
axis square
hold off
% legend(pl, well_regList)
grid on
xlabel('Spikes per Burst')
ylabel('Pr(X \geq x)')
% saveas(gcf, 'D:\Brewer lab data\HFS\figs_nostim\well_5SD-spnb-z1.png')
saveas(gcf, 'D:\Brewer lab data\HFS\figs_nostim\well_5SD_500maxSD-spnb-z1.png')

xlim([-inf 5e2]); ylim([1e-2 1])
xticks([6 12 25 50 100 150 250 500])
yticks([0.01 0.1 0.25 0.5 1])
% ax.FontSize = 16;
% saveas(gcf, 'D:\Brewer lab data\HFS\figs_nostim\well_5SD-spnb.png')
saveas(gcf, 'D:\Brewer lab data\HFS\figs_nostim\well_5SD_500maxSD-spnb.png')
hold off

% save('well_spnb_points_plotted_nostim_5SD','points_plotted')
% save('well_spnb_fit_plotted_nostim_5SD','plotted_fit')
save('well_spnb_points_plotted_nostim_5SD_500maxSD','points_plotted')
save('well_spnb_fit_plotted_nostim_5SD_500maxSD','plotted_fit')
%% ANOVA

[~,~,~,stats] = aoctool(sp_feature,prob_vec,regLabel,0.05,'','','','off');
[c,means]=multcompare(stats,0.05,'off','','s');
fit_table 

spnb_sp_feature_post_well=sp_feature;
spnb_prob_vec_post_well=prob_vec;
spnb_post_stim_string_vec=repmat("No Stim",length(regLabel),1);
spnb_regLabel_string=string(regLabel);
spnb_regLabel_post_well=categorical(join([spnb_regLabel_string,spnb_post_stim_string_vec]));
% save('no_stim_well_ancova_vars_spnb_5SD','spnb_sp_feature_post_well','spnb_prob_vec_post_well','spnb_regLabel_post_well')
% save('well_means_nostim_spnb_5SD','means')
save('no_stim_well_ancova_vars_spnb_5SD_500maxSD','spnb_sp_feature_post_well','spnb_prob_vec_post_well','spnb_regLabel_post_well')
save('well_means_nostim_spnb_5SD_500maxSD','means')
%% Error bar graph
figure( 'Position', [100 100 700 600])
xdata = 1:4;
ydata = means(:,1);
errdata = means(:,2);
width=0.35;
b = bar(xdata,ydata);
b.FaceColor = 'flat';
for regi=1:4
    b.CData(regi,:) = col_array_rgb(regi,:);
end
hold on
errorbar(xdata, ydata, errdata, '.k')
xticks(1:4), xticklabels(well_regList) 
ylim([-3.5 0])
yticks(linspace(-3.5,0,15))
ylabel 'Slope of spnb'
set(gca,'fontsize',17.5)
% saveas(gcf, 'D:\Brewer lab data\HFS\figs_nostim\well_5SD-spnb-bar.png')
saveas(gcf, 'D:\Brewer lab data\HFS\figs_nostim\well_5SD_500maxSD-spnb-bar.png')


%% Burst duration
moment_table = table([],[],[],[],[],[],[], 'VariableNames',{'regi','mean','median','mode','sd','emp_mode','emp_sd'});
figure( 'Position', [100 100 700 600])
%figure( 'Position', [100 100 2800 600])
binEdge = logspace(0,4,200); binCenter = convert_edges_2_centers(binEdge); 
for regi = 1:4
    coli=1; histcount=[];
    subplot(2,2,figOrder(regi))
    for fi=1:10
        vec2plot = cell2mat(well_spike_burst_dyn_table(... 
            well_spike_burst_dyn_table.regi==regi & well_spike_burst_dyn_table.fi==fi,:).BurstDuration);
        [hist_object,hist_plot,histcount(coli,:)] = log_binned_histogram(vec2plot, binEdge, "pdf",0);
        coli=coli+1;
    end
    histcount(~any(~isnan(histcount), 2),:)=[];
    meanC = mean(histcount,1); stdC = stdErr(histcount,1);
   
    % calculating mean and sd 
    f = fit(log(binCenter)', meanC', 'gauss1');
    
%     if regi==3
%         f = fit(log(binCenter)', meanC', 'gauss1','Exclude', find(ismember(meanC',maxk(meanC',2))));
%     end
    
    mu = f.b1; sigma = f.c1/sqrt(2);
    moment_table.regi(regi)=regi;
    
    %fitted
    moment_table.mean(regi) =exp(mu+(sigma^2)/2);
    moment_table.mode(regi) = exp(mu-(sigma^2)/2);
    moment_table.sd(regi) = sqrt((exp(sigma^2)-1)*(moment_table.mean(regi) ^2));
    moment_table.median(regi)=exp(mu);
    n=length(cell2mat(well_spike_burst_dyn_table(... 
            well_spike_burst_dyn_table.regi==regi,:).BurstDuration));
    moment_table.n(regi)=n;
    moment_table.se(regi) = sqrt((exp(sigma^2)-1)*(moment_table.mean(regi) ^2))/sqrt(n);
           
    % plotting data
    errorbar(binCenter, meanC, stdC,'k');
    title(well_regList(regi))
    xlabel 'Burst Duration (ms)', ylabel 'Probability'
    set(gca, 'xscale','log', 'fontsize',18)
    xlim([1e1 1e3]), ylim([0 0.04])
    format_2x2_plot(regi)
        
    % Overlaying fit
    hold on
    fitvals = f(log(binCenter)');
    pf = plot(binCenter, fitvals);
    pf.LineStyle = '--';
    pf.LineWidth = 3;
    pf.Color = 'r';
    moment_table.emp_mode(regi) = binCenter(fitvals == max(fitvals));
    moment_table.emp_sd(regi) = sqrt(var(binCenter,fitvals));
%     moment_table.n(regi)=length(fitvals);
    moment_table.lnmu(regi)=mu;
    moment_table.lnsd(regi)=sigma;
    moment_table.n_fit(regi)=length(meanC);
end
moment_table
% saveas(gcf, 'D:\Brewer lab data\HFS\figs_nostim\well_5SD-bd.png')
saveas(gcf, 'D:\Brewer lab data\HFS\figs_nostim\well_5SD_500maxSD-bd.png')
% save("nostim_BD_well_moment_table","moment_table")

%% Testing fit
gof_table = table([],[],'variablenames',{'regi','rsquare'}); rowi=1;
figure( 'Position', [100 100 700 600])
for regi = 1:4
    coli=1; histcount=[];
    subplot(2,2,figOrder(regi))
    for fi=1:10
        vec2plot = cell2mat( well_spike_burst_dyn_table(... 
            well_spike_burst_dyn_table.regi==regi & well_spike_burst_dyn_table.fi==fi,:).BurstDuration);
        [hist_object,hist_plot,histcount(coli,:)] = log_binned_histogram(vec2plot, binEdge, "pdf",0);
        coli=coli+1;
    end
    histcount(~any(~isnan(histcount), 2),:)=[];
    meanC = mean(histcount,1); stdC = stdErr(histcount,1);
    
    % calculating mean and sd 
    [f,gof] = fit(log(binCenter)', meanC', 'gauss1');
    
%     if regi==3
%         [f,gof] = fit(log(binCenter)', meanC', 'gauss1','Exclude', find(ismember(meanC',maxk(meanC',2))));
%     end
    
    mu = f.b1; sigma = f.c1/sqrt(2);       
    plot(f, log(binCenter), meanC)
    title(well_regList(regi))
    xlabel 'ln(Burst Duration (ms) )'; ylabel 'Probability'
    legend off
    format_2x2_plot(regi)
%     axis square
    ylim([0 0.05])
    
    % populating GOF table
    gof_table.regi(rowi)=regi;
    gof_table.rsquare(rowi)=gof.rsquare;
    rowi=rowi+1;
end
% saveas(gcf,'D:\Brewer lab data\HFS\figs_nostim\well_5SD-bd-test-fit.png')
saveas(gcf,'D:\Brewer lab data\HFS\figs_nostim\well_5SD_500maxSD-bd-test-fit.png')
gof_table
% save('nostim_well_BD_5SD','moment_table','gof_table')
save('nostim_well_BD_5SD_500maxSD','moment_table','gof_table')
%% Burst duration long format
moment_table = table([],[],[],[],[],[],[], 'VariableNames',{'regi','mean','median','mode','sd','emp_mode','emp_sd'});
%figure( 'Position', [100 100 700 600])
figure( 'Position', [100 100 2800 600])
binEdge = logspace(0,4,200); binCenter = convert_edges_2_centers(binEdge); 
tiledlayout(1,4, 'Padding', 'none', 'TileSpacing', 'compact');
for regi = 1:4
    coli=1; histcount=[];
    nexttile
    %subplot(1,4,(regi))
    for fi=1:10
        vec2plot = cell2mat(well_spike_burst_dyn_table(... 
            well_spike_burst_dyn_table.regi==regi & well_spike_burst_dyn_table.fi==fi,:).BurstDuration);
        [hist_object,hist_plot,histcount(coli,:)] = log_binned_histogram(vec2plot, binEdge, "pdf",0);
        coli=coli+1;
    end
    histcount(~any(~isnan(histcount), 2),:)=[];
    meanC = mean(histcount,1); stdC = stdErr(histcount,1);
   
    % calculating mean and sd 
    f = fit(log(binCenter)', meanC', 'gauss1');
    
%     if regi==3
%         f = fit(log(binCenter)', meanC', 'gauss1','Exclude', find(ismember(meanC',maxk(meanC',2))));
%     end
    
    mu = f.b1; sigma = f.c1/sqrt(2);
    moment_table.regi(regi)=regi;
    
    %fitted
    moment_table.mean(regi) =exp(mu+(sigma^2)/2);
    moment_table.mode(regi) = exp(mu-(sigma^2)/2);
    moment_table.sd(regi) = sqrt((exp(sigma^2)-1)*(moment_table.mean(regi) ^2));
    moment_table.median(regi)=exp(mu);
           
    % plotting data
    errorbar(binCenter, meanC, stdC,'k');
    title(well_regList(regi))
    xlabel 'Burst Duration (ms)', ylabel 'Probability'
    set(gca, 'xscale','log', 'fontsize',18)
    xlim([1e1 1e3]), ylim([0 0.04])
    %format_2x2_plot(regi)
        
    % Overlaying fit
    hold on
    fitvals = f(log(binCenter)');
    pf = plot(binCenter, fitvals);
    pf.LineStyle = '--';
    pf.LineWidth = 3;
    pf.Color = 'r';
    moment_table.emp_mode(regi) = binCenter(fitvals == max(fitvals));
    moment_table.emp_sd(regi) = sqrt(var(binCenter,fitvals));
    moment_table.n(regi)=length(fitvals);
    set(gca,'FontSize',24)
end
moment_table
% saveas(gcf, 'D:\Brewer lab data\HFS\figs_nostim\well_5SD-bd.png')
saveas(gcf, 'D:\Brewer lab data\HFS\figs_nostim\well_5SD_500maxSD-bd_long.png')
%% IBSR
moment_table = table([],[],[],[],[],[],[], 'VariableNames',{'regi','mean','median','mode','sd','emp_mode','emp_sd'});
figure( 'Position', [100 100 700 600])
binEdge = logspace(0,3,150); binCenter = convert_edges_2_centers(binEdge); 
for regi = 1:4
    coli=1; histcount=[];
    subplot(2,2,figOrder(regi))
    for fi=1:10
        vec2plot = cell2mat( well_spike_burst_dyn_table(... 
            well_spike_burst_dyn_table.regi==regi & well_spike_burst_dyn_table.fi==fi,:).IntraBurstSpikeRate);
        [hist_object,hist_plot,histcount(coli,:)] = log_binned_histogram(vec2plot, binEdge, "pdf",0);
        coli=coli+1;
    end
    histcount(~any(~isnan(histcount), 2),:)=[];
    meanC = mean(histcount,1); stdC = stdErr(histcount,1);
    
    % calculating mean and sd 
    f = fit(log(binCenter)', meanC', 'gauss1');
    
%     if regi==3
%         f = fit(log(binCenter)', meanC', 'gauss1','Exclude', find(ismember(meanC',maxk(meanC',2))));
%     end
    
    mu = f.b1; sigma = f.c1/sqrt(2);
    moment_table.regi(regi)=regi;
    
    %fitted
    moment_table.mean(regi) =exp(mu+(sigma^2)/2);
    moment_table.mode(regi) = exp(mu-((sigma^2)));
    moment_table.sd(regi) = sqrt((exp(sigma^2)-1)*(moment_table.mean(regi) ^2));
    moment_table.median(regi)=exp(mu);
    n=length(cell2mat(well_spike_burst_dyn_table(... 
            well_spike_burst_dyn_table.regi==regi,:).IntraBurstSpikeRate));
    moment_table.n(regi)=n;
    moment_table.se(regi) = sqrt((exp(sigma^2)-1)*(moment_table.mean(regi) ^2))/sqrt(n);
       
    errorbar(binCenter, meanC, stdC,'k');
    title(well_regList(regi))
    xlabel 'Intraburst spike rate (Hz)', ylabel 'Probability'
    set(gca, 'xscale','log', 'fontsize',18,'XMinorTick','off')
    ylim([0 0.06]), xlim([20 620])
    xticks([25 50 100 200 400])
    format_2x2_plot(regi)
        
    % Overlaying fit
    hold on
    fitvals = f(log(binCenter)');
    pf = plot(binCenter, fitvals);
    pf.LineStyle = '--';
    pf.LineWidth = 3;
    pf.Color = 'r';
    moment_table.emp_mode(regi) = binCenter(fitvals == max(fitvals));
    moment_table.emp_sd(regi) = sqrt(var(binCenter,fitvals));
%     moment_table.n(regi)=length(fitvals);
    moment_table.lnmu(regi)=mu;
    moment_table.lnsd(regi)=sigma;
    moment_table.n_fit(regi)=length(meanC);
end
% saveas(gcf,'D:\Brewer lab data\HFS\figs_nostim\well_5SD-ibsr.png')
saveas(gcf,'D:\Brewer lab data\HFS\figs_nostim\well_5SD_500maxSD-ibsr.png')
moment_table
% save("nostim_IBSR_well_moment_table","moment_table")
%% Testing fit
gof_table = table([],[],'variablenames',{'regi','rsquare'}); rowi=1;
figure( 'Position', [100 100 700 600])
binEdge = logspace(0,3,150); binCenter = convert_edges_2_centers(binEdge); 
for regi = 1:4
    coli=1; histcount=[];
    subplot(2,2,figOrder(regi))
    for fi=1:10
        vec2plot = cell2mat( well_spike_burst_dyn_table(... 
            well_spike_burst_dyn_table.regi==regi & well_spike_burst_dyn_table.fi==fi,:).IntraBurstSpikeRate);
        [hist_object,hist_plot,histcount(coli,:)] = log_binned_histogram(vec2plot, binEdge, "pdf",0);
        coli=coli+1;
    end
    histcount(~any(~isnan(histcount), 2),:)=[];
    meanC = mean(histcount,1); stdC = stdErr(histcount,1);
    
    % calculating mean and sd 
    [f,gof] = fit(log(binCenter)', meanC', 'gauss1');
    
%     if regi==3
%         [f,gof] = fit(log(binCenter)', meanC', 'gauss1','Exclude', find(ismember(meanC',maxk(meanC',2))));
%     end
    
    mu = f.b1; sigma = f.c1/sqrt(2);       
    plot(f, log(binCenter), meanC)
    title(well_regList(regi))
    xlabel 'ln(Burst Duration (ms) )'; ylabel 'Probability'
    legend off
    format_2x2_plot(regi)
%     axis square
    ylim([0 0.06])
    
    % populating GOF table
    gof_table.regi(rowi)=regi;
    gof_table.rsquare(rowi)=gof.rsquare;
    rowi=rowi+1;
end
%saveas(gcf,'D:\Brewer lab data\HFS\figs_nostim\well_5SD-ibsr-test-fit.png')
saveas(gcf,'D:\Brewer lab data\HFS\figs_nostim\well_5SD_500maxSD-ibsr-test-fit.png')
gof_table
% save('nostim_well_IBSR_5SD','moment_table','gof_table')
save('nostim_well_IBSR_5SD_500maxSD','moment_table','gof_table')
%% IBSR long format
moment_table = table([],[],[],[],[],[],[], 'VariableNames',{'regi','mean','median','mode','sd','emp_mode','emp_sd'});
figure( 'Position', [100 100 2800 600])
binEdge = logspace(0,3,150); binCenter = convert_edges_2_centers(binEdge);
tiledlayout(1,4, 'Padding', 'none', 'TileSpacing', 'compact');
for regi = 1:4
    coli=1; histcount=[];
    nexttile
    %subplot(1,4,(regi))
    for fi=1:10
        vec2plot = cell2mat( well_spike_burst_dyn_table(... 
            well_spike_burst_dyn_table.regi==regi & well_spike_burst_dyn_table.fi==fi,:).IntraBurstSpikeRate);
        [hist_object,hist_plot,histcount(coli,:)] = log_binned_histogram(vec2plot, binEdge, "pdf",0);
        coli=coli+1;
    end
    histcount(~any(~isnan(histcount), 2),:)=[];
    meanC = mean(histcount,1); stdC = stdErr(histcount,1);
    
    % calculating mean and sd 
    f = fit(log(binCenter)', meanC', 'gauss1');
    
%     if regi==3
%         f = fit(log(binCenter)', meanC', 'gauss1','Exclude', find(ismember(meanC',maxk(meanC',2))));
%     end
    
    mu = f.b1; sigma = f.c1/sqrt(2);
    moment_table.regi(regi)=regi;
    
    %fitted
    moment_table.mean(regi) =exp(mu+(sigma^2)/2);
    moment_table.mode(regi) = exp(mu-((sigma^2)));
    moment_table.sd(regi) = sqrt((exp(sigma^2)-1)*(moment_table.mean(regi) ^2));
    moment_table.median(regi)=exp(mu);
       
    errorbar(binCenter, meanC, stdC,'k');
    title(well_regList(regi))
    xlabel 'Intraburst spike rate (Hz)', ylabel 'Probability'
    set(gca, 'xscale','log', 'fontsize',18,'XMinorTick','off')
    ylim([0 0.06]), xlim([20 620])
    xticks([25 50 100 200 400])
    %format_2x2_plot(regi)
        
    % Overlaying fit
    hold on
    fitvals = f(log(binCenter)');
    pf = plot(binCenter, fitvals);
    pf.LineStyle = '--';
    pf.LineWidth = 3;
    pf.Color = 'r';
    moment_table.emp_mode(regi) = binCenter(fitvals == max(fitvals));
    moment_table.emp_sd(regi) = sqrt(var(binCenter,fitvals));
    moment_table.n(regi)=length(fitvals);
    set(gca,'FontSize',24)
end
% saveas(gcf,'D:\Brewer lab data\HFS\figs_nostim\well_5SD-ibsr.png')
saveas(gcf,'D:\Brewer lab data\HFS\figs_nostim\well_5SD_500maxSD-ibsr_long.png')
moment_table
