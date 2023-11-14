clear
clc
%load('D:\Brewer lab data\HFS\HFS Stim\24-Nov-2021_A\spike_burst_dyn_table_stim.mat') % 11 SD
load('D:\Brewer lab data\HFS\HFS Stim\03-Oct-2022_A\spike_burst_dyn_table_stim.mat') % 5 SD
cd 'D:\Brewer lab data\HFS'
tunnel_spike_burst_dyn_table = spike_burst_dyn_table_stim;
clear spike_burst_dyn_table
ff_regList = ["EC-DG","DG-CA3","CA3-CA1","CA1-EC"]; figOrder = [1 2 4 3];
fb_regList = ["DG-EC","CA3-DG","CA1-CA3","EC-CA1"]; 
powerlawfitfun = @(b,x) 10.^(b(2)*log10(x) + b(1)); 
col_array = ["#D95319","#77AC30","#4DBEEE","#7E2F8E"];
col_array_rgb = [0.8500 0.3250 0.0980; 0.4660 0.6740 0.1880; ...
    0.3010 0.7450 0.9330; 0.4940 0.1840 0.5560];
%% ISI histograms  
% FF histograms

%FF
fit_table = table([],[],[],[],'variablenames',{'regi','slope','intercept','rsquared'}); rowi=1;
figure( 'Position', [100 100 700 600])
binEdge = logspace(-2,1,100); binCenter = convert_edges_2_centers(binEdge);
% xlimits = [ 2e-2, 0.8e0]; %original
% xlimits = [ 1e-2, 2e-1]; % 11 SD
xlimits = [ 1e-2, 1e0]; %new
sp_feature=[]; prob_vec=[]; regLabel=[]; pl=[];
hold on
disp('CCDF for Interspike Interval for FF axons')
points_plotted=[];
for regi=1:4
    vec2plot = cell2mat(tunnel_spike_burst_dyn_table(tunnel_spike_burst_dyn_table.if_ff==1 &...
        tunnel_spike_burst_dyn_table.regi==regi,:).ISI);
    [ hist_object, hist_plot, prob, xbin] = log_binned_histogram( vec2plot, binEdge, "cdf" );  
    hist_plot.Color = col_array(regi);
    pl(regi)=hist_plot;
    
    % fitting data between [xmin, xmax]
    [ coeff, stats, opt_lim, logBinC,logHistProb ] = find_powerlawfit_with_grid_search(binCenter, ...
        prob, xlimits);
    if isempty(opt_lim)
        sp_feature = [sp_feature; NaN]; prob_vec=[prob_vec; NaN];
        regLabel = [regLabel; repmat(categorical(ff_regList(regi)),1,1)];
        continue
    end
    sp_feature = [sp_feature; logBinC']; prob_vec=[prob_vec; logHistProb']; 
    regLabel = [regLabel; repmat(categorical(ff_regList(regi)),length(logBinC),1)];
%     sp_feature = [sp_feature; log10(binCenter)']; prob_vec=[prob_vec; log10(prob)'];
%     regLabel = [regLabel; repmat(categorical(ff_regList(regi)),length(log10(binCenter)),1)];
    
     fitbinCenter = binCenter( binCenter > opt_lim(1) & binCenter< opt_lim(end));
%    fitbinCenter = binCenter;
    regfit = powerlawfitfun(coeff,fitbinCenter);
    plot(fitbinCenter , regfit ,'--k','LineWidth',2)
    fit_table.regi(rowi) = regi;
    fit_table.intercept(rowi) = coeff(1);
    fit_table.slope(rowi) = coeff(2);
    fit_table.rsquared(rowi) = stats.Rsquared;
    rowi = rowi+1;
    points_plotted{regi}=[{prob},{xbin}];
    plotted_fit{regi}=[{fitbinCenter},{regfit}];
end
ax = gca;
ax.FontSize = 20;
ax.YScale = 'log';
ax.XLim = [-inf 2];
ax.YLim=[1e-2 1];
set(gca,'linewidth',2)
axis square
hold off
%legend(pl, ff_regList)
grid on
saveas(gcf, './figs_HFS/ff-isi-z1.png')
xlabel('Interspike Interval (s)')
ylabel('Pr(X \geq x)')
xlim([1e-2 2e-1]); %original
%xlim([1e-2 2])
ylim([1e-1 1])
xticks([1e-2, 0.25e-1,0.5e-1,0.1,0.2,0.3,0.6e0,1,2])
yticks([1e-1, 0.25, 0.5, 1])
set(gca,'YScale','log')
set(gca,'XScale','log')
title("ISI Feedforward")
ax.FontSize = 16;
saveas(gcf, './figs_HFS/ff-isi.png')
save('ff_ISI_points_plotted_HFS','points_plotted')
save('ff_ISI_fit_plotted_HFS','plotted_fit')
%% FF stats and ANCOVA

fit_table 
%ANOVA
[~,~,~,stats] = aoctool(sp_feature,prob_vec,regLabel,0.05,'','','','off');
[ff_c,ff_means]=multcompare(stats,0.05,'on','','s');


%temporary block of code
%ff_means=[ff_means(1,:);0,0;ff_means(2:end,:)];

sp_feature_post_ff=sp_feature;
prob_vec_post_ff=prob_vec;
post_stim_string_vec=repmat("HFS",length(regLabel),1);
regLabel_string=string(regLabel);
regLabel_post_ff=categorical(join([regLabel_string,post_stim_string_vec]));
title('Feedforward ISI Slope ANCOVA')
xlabel('Slope')
saveas(gcf,'./figs_HFS/ff-isi-ancova.png')
save('HFS_ff_ancova_vars','sp_feature_post_ff','prob_vec_post_ff','regLabel_post_ff')
save('ff_means_HFS','ff_means')
%% FB Histograms

% FB
fit_table = table([],[],[],[],'variablenames',{'regi','slope','intercept','rsquared'}); rowi=1;
figure( 'Position', [100 100 700 600])
binEdge = logspace(-2,1,100); binCenter = convert_edges_2_centers(binEdge);
% Play around with these values to see what happens, talk with Ruiyi
%xlimits = [ 2e-2, 12e-1];
% xlimits = [ 2e-2, 0.8e0]; %original
% xlimits = [ 1e-2, 2e-1]; % 11 SD
xlimits = [ 1e-2, 1e0]; %new
sp_feature=[]; prob_vec=[]; regLabel=[]; pl=[];
hold on
disp('CCDF for Interspike Interval for FB axons')
for regi=1:4
    vec2plot = cell2mat(tunnel_spike_burst_dyn_table(tunnel_spike_burst_dyn_table.if_ff==0 &...
        tunnel_spike_burst_dyn_table.regi==regi,:).ISI);
    [ hist_object, hist_plot, prob, xbin] = log_binned_histogram( vec2plot, binEdge, "cdf" );  
    hist_plot.Color = col_array(regi);
    pl(regi)=hist_plot;
    
    % fitting data between [xmin, xmax]
    [ coeff, stats, opt_lim, logBinC,logHistProb ] = find_powerlawfit_with_grid_search(binCenter, ...
        prob, xlimits);
    if isempty(opt_lim)
        sp_feature = [sp_feature; NaN]; prob_vec=[prob_vec; NaN];
        regLabel = [regLabel; repmat(categorical(ff_regList(regi)),1,1)];
        continue
    end
    sp_feature = [sp_feature; logBinC']; prob_vec=[prob_vec; logHistProb']; 
    regLabel = [regLabel; repmat(categorical(ff_regList(regi)),length(logBinC),1)];
%     sp_feature = [sp_feature; log10(binCenter)']; prob_vec=[prob_vec; log10(prob)'];
%     regLabel = [regLabel; repmat(categorical(ff_regList(regi)),length(log10(binCenter)),1)];    
    
    fitbinCenter = binCenter( binCenter > opt_lim(1) & binCenter< opt_lim(end));
%     fitbinCenter = binCenter;
    regfit = powerlawfitfun(coeff,fitbinCenter);
    plot(fitbinCenter , regfit ,'--k','LineWidth',2)
    fit_table.regi(rowi) = regi;
    fit_table.intercept(rowi) = coeff(1);
    fit_table.slope(rowi) = coeff(2);
    fit_table.rsquared(rowi) = stats.Rsquared;
    rowi = rowi+1;
    points_plotted{regi}=[{prob},{xbin}];
    plotted_fit{regi}=[{fitbinCenter},{regfit}];
end
ax = gca;
ax.FontSize = 20;
ax.YScale = 'log';
ax.XLim = [-inf 5];
set(gca,'linewidth',2)
axis square
hold off
%legend(pl, fb_regList)
grid on
saveas(gcf, './figs_HFS/fb-isi-z1.png')
xlabel('Interspike Interval (s)')
ylabel('Pr(X \geq x)')
xlim([1e-2 2e-1]); %original
%xlim([1e-2,2])
ylim([1e-1 1])
xticks([1e-2, 0.25e-1,0.5e-1,0.1,0.2,0.3,0.6e0,1,2])
yticks([1e-1, 0.25, 0.5, 1])
set(gca,'YScale','log')
set(gca,'XScale','log')
title("ISI Feedbackward")
ax.FontSize = 16;
saveas(gcf, './figs_HFS/fb-isi.png')
save('fb_ISI_points_plotted_HFS','points_plotted')
save('fb_ISI_fit_plotted_HFS','plotted_fit')
%% FB stats and ANCOVA

fit_table 
[~,~,~,stats] = aoctool(sp_feature,prob_vec,regLabel,0.05,'','','','off');
[fb_c,fb_means]=multcompare(stats,0.05,'on','','s');

%Saving for ancova in a different script
sp_feature_post_fb=sp_feature;
prob_vec_post_fb=prob_vec;
post_stim_string_vec=repmat("HFS",length(regLabel),1);
regLabel_string=string(regLabel);
regLabel_post_fb=categorical(join([regLabel_string,post_stim_string_vec]));
title('Feedback ISI Slope ANCOVA')
xlabel('Slope')
saveas(gcf,'./figs_poststim/fb-isi-ancova.png')
save('HFS_fb_ancova_vars','sp_feature_post_fb','prob_vec_post_fb','regLabel_post_fb')
save('fb_means_HFS','fb_means')
%% FF error bar graph

figure( 'Position', [100 100 700 600])
xdata = 1:4;
ydata = [ff_means(:,1); fb_means(:,1)]';
errdata = [ff_means(:,2); fb_means(:,2)]';

b = bar(xdata,ydata(1:4));
b.FaceColor = 'flat';
for regi=1:4
    b.CData(regi,:) = col_array_rgb(regi,:);
end

hold on
errorbar(1:4, ydata(1:4), errdata(1:4), '.k')
xticks(1:4), xticklabels(ff_regList) 
ylabel 'Slope of ISI (s^{-1})'
set(gca,'fontsize',17.5)
set(gca,'linewidth',2)
hold off
ylim([-1 0])
yticks([-1.0 -0.9 -0.8, -0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0])
title("Feedforward ISI Error")
saveas(gcf,'./figs_HFS/ff-isi-bar.png')
%% FB error bar graph

figure( 'Position', [0 0 700 600])
b = bar(xdata,ydata(5:end));
b.FaceColor = 'flat';
for regi=1:4
    b.CData(regi,:) = col_array_rgb(regi,:);
end

hold on
errorbar(1:4, ydata(5:end), errdata(5:end), '.k')
xticks(1:4), xticklabels(fb_regList) 
ylabel 'Slope of ISI (s^{-1})'
set(gca,'fontsize',17.5)
set(gca,'linewidth',2)
hold off
%yticks(linspace(-1.5,0,7))
ylim([-1 0])
yticks([-1 -0.9 -0.8, -0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0])
% title("Feedback ISI Error")
saveas(gcf,'./figs_HFS/fb-isi-bar.png')
%% IBI histograms
% FF histograms

disp 'CCDF for Interburst Interval for FF axons'
fit_table = table([],{},{},{},{},'variablenames',{'regi','slope','intercept','rsquared','limits'}); rowi=1;
figure( 'Position', [100 100 700 600])
ylim([1e-2,1])%cdf
%ylim([1,1000])%count
%xlim([1e-1,5e1])
xlim([1e-1,2.5e1])
% 11 SD
% xlimits = [1e-1, 5; 1e-1, 2; 1e-1 2.5; 1e-1 0.75]; x_limit_final = [2.5e1,15,20,7];
% x_poi = [4.9811,2.12803,2.60566,1.47804]
% 5 SD
xlimits = [1e-1, 2; 1e-1, 2; 1e-1 2.5; 1e-1 3]; x_limit_final = [2.5e1,15,20,20];
x_poi = [1.96247,2.12803,2.21598,3.32235]
% Test if x_poi can be replaced with second column of xlimits, may be
% redundant
sp_feature_1=[]; prob_vec_1=[]; regLabel_1=[]; pl=[];
sp_feature_2=[]; prob_vec_2=[]; regLabel_2=[];
hold on
binEdge = logspace(-1,2.5,200); binCenter = convert_edges_2_centers(binEdge);
for regi=1:4
    vec2plot = cell2mat(tunnel_spike_burst_dyn_table( tunnel_spike_burst_dyn_table.if_ff == 1 & tunnel_spike_burst_dyn_table.regi==regi,:).IBI);
    [ hist_object, hist_plot, Prob, xbin] = log_binned_histogram( vec2plot, binEdge, "cdf" );  
    hist_plot.Color = col_array(regi);
    pl(regi)=hist_plot;
    
    points_plotted{regi}=[{xbin},{Prob}];
    
    % fitting data between [xmin_lower, xmax_lower] for lower fit
    [ low_coeff, low_stats, opt_lim, logBinC,logHistProb ] = find_powerlawfit_with_grid_search(binCenter, ...
        Prob, xlimits(regi,:));
    %modification for low population vec2plot by Sam Lassers 9/1/21
    if ~isempty(opt_lim)
        sp_feature_1 = [sp_feature_1; logBinC']; prob_vec_1=[prob_vec_1; logHistProb']; 
        regLabel_1 = [regLabel_1; repmat(categorical(ff_regList(regi)),length(logBinC),1)];
    %     fitbinCenter = binCenter(binCenter > opt_lim(1) & binCenter < opt_lim(2));
    %     fitbinCenter = binCenter;
        fitbinCenter = binCenter(binCenter > opt_lim(1) & binCenter < x_poi(regi));
        regfit = powerlawfitfun(low_coeff,fitbinCenter);
        plot(fitbinCenter , regfit ,'--k','LineWidth',3)
        
        plotted_fit_1{regi}=[{fitbinCenter},{regfit}];
    
        if regi>0
            % fitting data between [x_max_lower, abs_xmax]
            [ high_coeff, high_stats, opt_lim_high, logBinC,logHistProb ] = find_powerlawfit_with_grid_search(binCenter, ...
                Prob, [max(opt_lim), x_limit_final(regi)],[0.5 0.5]  );
            sp_feature_2 = [sp_feature_2; logBinC']; prob_vec_2=[prob_vec_2; logHistProb']; 
            regLabel_2 = [regLabel_2; repmat(categorical(ff_regList(regi)),length(logBinC),1)];
    %         fitbinCenter = binCenter(binCenter > opt_lim_high(1) & binCenter < opt_lim_high(2));
    %         fitbinCenter = binCenter;
            fitbinCenter = binCenter(binCenter > x_poi(regi) & binCenter < opt_lim_high(2));
            regfit = powerlawfitfun(high_coeff,fitbinCenter);
            plot(fitbinCenter , regfit ,':','LineWidth',3,'Color',[0.5 0.5 0.5])
        else
            high_coeff = [0,0];
            high_stats.Rsquared = 0;
        end
        plotted_fit_2{regi}=[{fitbinCenter},{regfit}];
    end
%   
    if ~isempty(opt_lim)
        fit_table.regi(rowi) = regi;
        fit_table.intercept{rowi} = [low_coeff(1), high_coeff(1)];
        fit_table.slope{rowi} = [low_coeff(2), high_coeff(2)];
        fit_table.rsquared{rowi} = [low_stats.Rsquared, high_stats.Rsquared];
        fit_table.limits{rowi} = [opt_lim(2),opt_lim_high(2)];
    else
        fit_table.regi(rowi) = regi;
        fit_table.intercept{rowi} = [];
        fit_table.slope{rowi} = [];
        fit_table.rsquared{rowi} = [];
        fit_table.limits{rowi} = [];
    end
    rowi = rowi+1;
end
% full scale
ax = gca;
ax.FontSize = 20;
ax.YScale = 'log';
set(gca,'linewidth',2)
% ax.LineWidth = 2*ax.LineWidth;
% ax.TickLength = 4*ax.TickLength;
%ax.XLim = [-inf 1e2]; ax.YLim = [1e-2, 1];
xticks([0.1,0.25,1,2.5,5,10,20,100])
yticks([1e-2 1e-1, 0.25, 0.5, 1])%cdf
%yticks([1,5,10,50,100, 200, 500, 1000]) %count
axis square
grid on
hold off
%legend(pl, ff_regList)
xlabel('Interburst Interval (s)')
%ylabel('Pr(X \geq x)')%cdf
ylabel('Spike Count Per Bin')%count
title('Feedforward IBI')
saveas(gcf, './figs_HFS/ff-ibi.png')
% % lower fits
% ax.XLim = [-inf 1e0]; ax.YLim = [0.5, 1];
% xticks([0.1, 0.25 0.5 1])
% yticks([0.1, 0.25 0.5 1])
% yticks(0.5:0.1:1)
% saveas(gcf, './figs_stim/ff-ibi-z1.png')
% % upper fits
% ax.XLim = [1e0 1.1e1]; ax.YLim = [1e-2, 1e-1];
% xticks([0.1, 0.25 0.5 1])
% yticks([1e-2, 0.25e-1, 0.5e-1, 1e-1]);
% xticks([1e0, 0.25e1, 0.5e1, 1e1])
% saveas(gcf, './figs_stim/ff-ibi-z2.png')

%delete empty rows
del_fit_row=[];
for i=1:4
    if isempty(fit_table(i,:).slope{1})
        del_fit_row=[del_fit_row,i];
    end
end

fit_table(del_fit_row,:)=[];
x_poi(del_fit_row)=[];

fit_table.x_poi = x_poi';
% lower fits FF stats and ANCOVA

save('ff_IBI_points_plotted_HFS','points_plotted')
save('ff_IBI_fit_plotted_HFS_1','plotted_fit_1')
save('ff_IBI_fit_plotted_HFS_2','plotted_fit_2')

% save('ff_IBI_points_plotted_poststim_count','points_plotted')
% save('ff_IBI_fit_plotted_poststim_1_count','plotted_fit_1')
% save('ff_IBI_fit_plotted_poststim_2_count','plotted_fit_2')
%% ANOVA
[~,~,~,stats] = aoctool(sp_feature_1,prob_vec_1,regLabel_1,0.05,'','','','off');
[ff_c,ff_means_lower]=multcompare(stats,0.05,'on','','s');
% higher fits FF stats and ANCOVA
save("IBI_hfs1_ffmeans","ff_means_lower")

sp_feature_post_ff_1=sp_feature_1;
prob_vec_post_ff_1=prob_vec_1;
pre_stim_string_vec=repmat("HFS",length(regLabel_1),1);
regLabel_string=string(regLabel_1);
regLabel_post_ff_1=categorical(join([regLabel_string,pre_stim_string_vec]));
save('IBI_HFS_ff_ancova_vars_1','sp_feature_post_ff_1','prob_vec_post_ff_1','regLabel_post_ff_1')
%% ANOVA
[~,~,~,stats] = aoctool(sp_feature_2,prob_vec_2,regLabel_2,0.05,'','','','off');
[ff_c,ff_means_higher]=multcompare(stats,0.05,'on','','s');

save("IBI_HFS2_ffmeans","ff_means_higher")

sp_feature_post_ff_2=sp_feature_2;
prob_vec_post_ff_2=prob_vec_2;
pre_stim_string_vec=repmat("HFS",length(regLabel_2),1);
regLabel_string=string(regLabel_2);
regLabel_post_ff_2=categorical(join([regLabel_string,pre_stim_string_vec]));
save('IBI_HFS_ff_ancova_vars_2','sp_feature_post_ff_2','prob_vec_post_ff_2','regLabel_post_ff_2')
%% FB histograms

disp 'CCDF for Interburst Interval for FB axons'
fit_table = table([],{},{},{},{},'variablenames',{'regi','slope','intercept','rsquared','limits'}); rowi=1;
figure( 'Position', [100 100 700 600])
ylim([1e-2,1])%cdf
%ylim([1,1000])
xlim([1e-1,25])
low_limit = 1.1e0;
% 11 SD
% xlimits = [1e-1, 8; 1e-1, 3; 1e-1 2.5; 1e-1 4]; x_limit_final = [30,20,10,4e1];
% x_poi=[9.14417,2.30757,2.30757,3.06387];
%5SD
xlimits = [1e-1, 3; 1e-1, 3; 1e-1 2.5; 1e-1 4]; x_limit_final = [30,30,30,4e1];
x_poi=[1.96247,2.21598,2.30757,2.21598];
sp_feature_1=[]; prob_vec_1=[]; regLabel_1=[]; pl=[];
sp_feature_2=[]; prob_vec_2=[]; regLabel_2=[];
hold on
binEdge = logspace(-1,2.5,200); binCenter = convert_edges_2_centers(binEdge);
for regi=1:4
    vec2plot = cell2mat(tunnel_spike_burst_dyn_table( tunnel_spike_burst_dyn_table.if_ff == 0 & tunnel_spike_burst_dyn_table.regi==regi,:).IBI);
    [ hist_object, hist_plot, Prob, xbin] = log_binned_histogram( vec2plot, binEdge, "cdf" );  
    hist_plot.Color = col_array(regi);
    pl(regi)=hist_plot;
    
    points_plotted{regi}=[{xbin},{Prob}];
    
    % fitting data between [xmin_lower, xmax_lower] for lower fit
    [ low_coeff, low_stats, opt_lim, logBinC,logHistProb ] = find_powerlawfit_with_grid_search(binCenter, ...
        Prob, xlimits(regi,:));
    if ~isempty(opt_lim)
        sp_feature_1 = [sp_feature_1; logBinC']; prob_vec_1=[prob_vec_1; logHistProb']; 
        regLabel_1 = [regLabel_1; repmat(categorical(ff_regList(regi)),length(logBinC),1)];
    %     fitbinCenter = binCenter(binCenter > opt_lim(1) & binCenter < opt_lim(2));
        fitbinCenter = binCenter(binCenter > opt_lim(1) & binCenter < x_poi(regi));
    %     fitbinCenter = binCenter;
        regfit = powerlawfitfun(low_coeff,fitbinCenter);
        plot(fitbinCenter , regfit ,'--k','LineWidth',3)
        plotted_fit_1{regi}=[{fitbinCenter},{regfit}];
        
        if regi>0
            % fitting data between [x_max_lower, abs_xmax]
            [ high_coeff, high_stats, opt_lim_high, logBinC,logHistProb ] = find_powerlawfit_with_grid_search(binCenter, ...
                Prob, [max(opt_lim), x_limit_final(regi)],[0.5 0.5]  );
            sp_feature_2 = [sp_feature_2; logBinC']; prob_vec_2=[prob_vec_2; logHistProb']; 
            regLabel_2 = [regLabel_2; repmat(categorical(ff_regList(regi)),length(logBinC),1)];
    %         fitbinCenter = binCenter(binCenter > opt_lim_high(1) & binCenter < opt_lim_high(2));
            fitbinCenter = binCenter(binCenter > x_poi(regi) & binCenter < opt_lim_high(2));
    %         fitbinCenter = binCenter;
            regfit = powerlawfitfun(high_coeff,fitbinCenter);
            plot(fitbinCenter , regfit ,':','LineWidth',3,'Color',[0.5 0.5 0.5])
        else
            high_coeff = [0,0];
            high_stats.Rsquared = 0;
        end
        plotted_fit_2{regi}=[{fitbinCenter},{regfit}];
    end
%     
    if ~isempty(opt_lim)
        fit_table.regi(rowi) = regi;
        fit_table.intercept{rowi} = [low_coeff(1), high_coeff(1)];
        fit_table.slope{rowi} = [low_coeff(2), high_coeff(2)];
        fit_table.rsquared{rowi} = [low_stats.Rsquared, high_stats.Rsquared];
        fit_table.limits{rowi} = [opt_lim(2),opt_lim_high(2)];
    else
        fit_table.regi(rowi) = regi;
        fit_table.intercept{rowi} = [];
        fit_table.slope{rowi} = [];
        fit_table.rsquared{rowi} = [];
        fit_table.limits{rowi} = [];
    end
    
    rowi = rowi+1;
end
% full scale
ax = gca;
ax.FontSize = 20;
ax.YScale = 'log';
set(gca,'linewidth',2)
% ax.LineWidth = 2*ax.LineWidth;
% ax.TickLength = 4*ax.TickLength;
%ax.XLim = [-inf 1e2]; ax.YLim = [1e-2, 1];
xticks([0.1,0.25,1,2.5,5,10,20,100])
yticks([1e-2, 1e-1, 0.25, 0.5, 1])%cdf
%yticks([1,5,10,50,100, 200, 500, 1000]) %count
axis square
hold off
%legend(pl, fb_regList)
grid on
xlabel('Interburst Interval (s)')
ylabel('Pr(X \geq x)')%cdf
%ylabel('Spike Count Per Bin')%count
title('Feedback IBI')
% grid on
saveas(gcf, './figs_HFS/fb-ibi.png')
% % lower fits
% ax.XLim = [-inf 1e0]; ax.YLim = [0.5, 1];
% xticks([0.1, 0.25 0.5 1])
% yticks(0.5:0.1:1)
% saveas(gcf, './figs_stim/fb-ibi-z1.png')
% % upper fits
% ax.XLim = [1e0 1.1e1]; ax.YLim = [1e-2, 1e-1];
% xticks([0.1, 0.25 0.5 1])
% yticks([1e-2, 0.25e-1, 0.5e-1, 1e-1]);
% xticks([1e0, 0.25e1, 0.5e1, 1e1])
% saveas(gcf, './figs_stim/fb-ibi-z2.png')

%delete empty rows
del_fit_row=[];
for i=1:4
    if isempty(fit_table(i,:).slope{1})
        del_fit_row=[del_fit_row,i];
    end
end

fit_table(del_fit_row,:)=[];
x_poi(del_fit_row)=[];

fit_table.x_poi = x_poi'

save('fb_IBI_points_plotted_HFS','points_plotted')
save('fb_IBI_fit_plotted_HFS_1','plotted_fit_1')
save('fb_IBI_fit_plotted_HFS_2','plotted_fit_2')

% save('fb_IBI_points_plotted_poststim_count','points_plotted')
% save('fb_IBI_fit_plotted_poststim_1_count','plotted_fit_1')
% save('fb_IBI_fit_plotted_poststim_2_count','plotted_fit_2')
%% lower fits FB stats and ANCOVA

fit_table 
[~,~,~,stats] = aoctool(sp_feature_1,prob_vec_1,regLabel_1,0.05,'','','','off');
[fb_c,fb_means_lower]=multcompare(stats,0.05,'on','','s');
% higher fits FB stats and ANCOVA
save("IBI_HFS1_fbmeans","fb_means_lower")

sp_feature_post_fb_1=sp_feature_1;
prob_vec_post_fb_1=prob_vec_1;
pre_stim_string_vec=repmat("HFS",length(regLabel_1),1);
regLabel_string=string(regLabel_1);
regLabel_post_fb_1=categorical(join([regLabel_string,pre_stim_string_vec]));
save('IBI_HFS_fb_ancova_vars_1','sp_feature_post_fb_1','prob_vec_post_fb_1','regLabel_post_fb_1')
%% ANOVA
[~,~,~,stats] = aoctool(sp_feature_2,prob_vec_2,regLabel_2,0.05,'','','','off');
[fb_c,fb_means_higher]=multcompare(stats,0.05,'on','','s');
save("IBI_HFS2_fbmeans","fb_means_higher")

sp_feature_post_fb_2=sp_feature_2;
prob_vec_post_fb_2=prob_vec_2;
pre_stim_string_vec=repmat("HFS",length(regLabel_2),1);
regLabel_string=string(regLabel_2);
regLabel_post_fb_2=categorical(join([regLabel_string,pre_stim_string_vec]));
save('IBI_HFS_fb_ancova_vars_2','sp_feature_post_fb_2','prob_vec_post_fb_2','regLabel_post_fb_2')
%% lower fits FF error bar graph

figure( 'Position', [100 100 700 600])
xdata = 1:4;
ydata = [ff_means_lower(:,1); fb_means_lower(:,1)]';
errdata = [ff_means_lower(:,2); fb_means_lower(:,2)]';
b = bar(xdata,ydata(1:4));
b.FaceColor = 'flat';
for regi=1:4
    b.CData(regi,:) = col_array_rgb(regi,:);
end
hold on
errorbar(1:4, ydata(1:4), errdata(1:4), '.k')
xticks(1:4), xticklabels(ff_regList) 
ylabel 'Slope of IBI (s^{-1})'
set(gca,'fontsize',17.5)
hold off
ylim([-1 0])
ytick_nums = [0.1, 0.25, 0.5, 1];
title("Poststim Feedforward Lower Fits")
% yticks(log10(ytick_nums)); yticklabels(ytick_nums)
saveas(gcf,'./figs_hfs/ff-ibi-bar.png')
% lower FB error bar graph

figure( 'Position', [100 100 700 600])
b = bar(xdata,ydata(5:end));
b.FaceColor = 'flat';
for regi=1:4
    b.CData(regi,:) = col_array_rgb(regi,:);
end
hold on
errorbar(1:4, ydata(5:end), errdata(5:end), '.k')
xticks(1:4), xticklabels(fb_regList) 
ylabel 'Slope of IBI (s^{-1})'
set(gca,'fontsize',17.5)
hold off
ylim([-1 0])
ytick_nums = [0.1, 0.25, 0.5, 1];
title("Poststim Feedback Lower Fits")
% yticks(log10(ytick_nums)); yticklabels(ytick_nums)
saveas(gcf,'./figs_hfs/fb-ibi-bar.png')

%% higher fits FF error bar graph

figure( 'Position', [100 100 700 600])
xdata = 1:4;
ydata = [ff_means_higher(:,1); fb_means_higher(:,1)]';
errdata = [ff_means_higher(:,2); fb_means_higher(:,2)]';
b = bar(xdata,ydata(1:4));
b.FaceColor = 'flat';
for regi=1:4
    b.CData(regi,:) = col_array_rgb(regi,:);
end
hold on
errorbar(1:4, ydata(1:4), errdata(1:4), '.k')
xticks(1:4), xticklabels(ff_regList) 
ylabel 'Slope of IBI (s^{-1})'
set(gca,'fontsize',17.5)
hold off
ylim([-2.5 0])
ytick_nums = [0.005, 0.01,0.025, 0.05, 0.1, 0.25, 0.5, 1];
% yticks(log10(ytick_nums)); yticklabels(ytick_nums)
title("Feedforward Higher Fits")
saveas(gcf,'./figs_hfs/ff-ibi-bar-2.png')
% higher FB error bar graph

figure( 'Position', [100 100 700 600])
b = bar(xdata,ydata(5:end));
b.FaceColor = 'flat';
for regi=1:4
    b.CData(regi,:) = col_array_rgb(regi,:);
end
hold on
errorbar(1:4, ydata(5:end), errdata(5:end), '.k')
xticks(1:4), xticklabels(fb_regList) 
ylabel 'Slope of IBI (s^{-1})'
set(gca,'fontsize',17.5)
hold off
ylim([-2.5 0])
ytick_nums = [0.005, 0.01,0.025, 0.05, 0.1, 0.25, 0.5, 1];
title("Poststim Feedback Higher Fits")
% yticks(log10(ytick_nums)); yticklabels(ytick_nums)
saveas(gcf,'./figs_hfs/fb-ibi-bar-2.png')

%% SPNB histograms
% FF Histogram

disp 'CCDF for spikes per burst for FF axons'
% xstart = [3,4,4,4]; xstop = [40,25,2e1,40]; % 11 SD
xstart = [4,4,4,4]; xstop = [40,50,50,50];
fit_table = table([],[],[],[],'variablenames',{'regi','slope','intercept','rsquared'}); rowi=1;
figure( 'Position', [100 100 700 600])
sp_feature=[]; prob_vec=[]; regLabel=[]; pl=[];
hold on
binEdge = logspace(log10(4),3,50); binCenter = convert_edges_2_centers(binEdge);
for regi=1:4
    vec2plot = cell2mat(tunnel_spike_burst_dyn_table(tunnel_spike_burst_dyn_table.if_ff==1 & tunnel_spike_burst_dyn_table.regi==regi,:).SpikeperBurst);
    [ hist_object, hist_plot, prob, xbin] = log_binned_histogram( vec2plot, binEdge, "cdf" );  
    hist_plot.Color = col_array(regi);
    pl(regi)=hist_plot;
    
    points_plotted{regi}=[{prob},{xbin}];
    
    % fitting data between [xmin, xmax]
    xlimits = [xstart(regi), xstop(regi)];
    [ coeff, stats, opt_lim, logBinC,logHistProb ] = find_powerlawfit_with_grid_search(binCenter, ...
        prob, xlimits);
    sp_feature = [sp_feature; logBinC']; prob_vec=[prob_vec; logHistProb']; 
    regLabel = [regLabel; repmat(categorical(ff_regList(regi)),length(logBinC),1)];
     fitbinCenter = binCenter( binCenter > opt_lim(1) & binCenter< opt_lim(end));
    %fitbinCenter = binCenter;
    regfit = powerlawfitfun(coeff,fitbinCenter);
    plot(fitbinCenter , regfit ,'--k')
    
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
ax.YLim = [1e-2 1]; ax.XLim = [-inf 500];
xticks([6 12 25 50 100 150 250 500])
yticks([0.01 0.1 0.25 0.5 1])
axis square
hold off
% legend(pl, ff_regList)
% grid on
xlabel('Spikes per Burst')
ylabel('Pr(X \geq x)')
%axis([4 55 0.1 1]), yticks([0.1,0.25,0.5,1])
saveas(gcf, './figs_hfs/ff-spnb-z1.png')

xlim([-inf 50]); ylim([1e-1 1])
%ax.YTick = [1e-2, 1e-1, 1];
%ax.XTick = [1e1, 1e2, 1e3];
% ax.FontSize = 16;
saveas(gcf, './figs_hfs/ff-spnb.png')
hold off
save('ff_spnb_points_plotted_hfs','points_plotted')
save('ff_spnb_fit_plotted_hfs','plotted_fit')
%% ANOVA

fit_table 
%ANOVA
[~,~,~,stats] = aoctool(sp_feature,prob_vec,regLabel,0.05,'','','','off');
[ff_c,ff_means]=multcompare(stats,0.05,'off','','s');
ff_means

spnb_sp_feature_post_ff=sp_feature;
spnb_prob_vec_post_ff=prob_vec;
spnb_post_stim_string_vec=repmat("HFS",length(regLabel),1);
spnb_regLabel_string=string(regLabel);
spnb_regLabel_post_ff=categorical(join([spnb_regLabel_string,spnb_post_stim_string_vec]));
%title('Feedforward ISI Slope ANCOVA')
%xlabel('Slope')
saveas(gcf,'./figs_HFS/ff-isi-ancova_spnb.png')
save('HFS_ff_ancova_vars_spnb','spnb_sp_feature_post_ff','spnb_prob_vec_post_ff','spnb_regLabel_post_ff')
save('ff_means_HFS_spnb','ff_means')
%% FB Histogram

disp 'CCDF for spikes per burst for FF axons'
% xstart = [4,4,4,4]; xstop = [18,50,2e1,1e2]; %11SD
xstart = [4,4,4,4]; xstop = [50,50,50,50]; %5SD
fit_table = table([],[],[],[],'variablenames',{'regi','slope','intercept','rsquared'}); rowi=1;
sp_feature=[]; prob_vec=[]; regLabel=[]; pl=[];
figure( 'Position', [100 100 700 600])
hold on
binEdge = logspace(log10(4),3,50); binCenter = convert_edges_2_centers(binEdge);
for regi=1:4
    vec2plot = cell2mat(tunnel_spike_burst_dyn_table(tunnel_spike_burst_dyn_table.if_ff==0 & tunnel_spike_burst_dyn_table.regi==regi,:).SpikeperBurst);
    [ hist_object, hist_plot, prob, xbin] = log_binned_histogram( vec2plot, binEdge, "cdf" );  
    hist_plot.Color = col_array(regi);
    pl(regi)=hist_plot;
    
    points_plotted{regi}=[{prob},{xbin}];
    
    % fitting data between [xmin, xmax]
    xlimits = [xstart(regi), xstop(regi)];
    [ coeff, stats, opt_lim, logBinC,logHistProb ] = find_powerlawfit_with_grid_search(binCenter, ...
        prob, xlimits);
    sp_feature = [sp_feature; logBinC']; prob_vec=[prob_vec; logHistProb']; 
    regLabel = [regLabel; repmat(categorical(ff_regList(regi)),length(logBinC),1)];
     fitbinCenter = binCenter( binCenter > opt_lim(1) & binCenter< opt_lim(end));
    %fitbinCenter = binCenter;
    regfit = powerlawfitfun(coeff,fitbinCenter);
    plot(fitbinCenter , regfit ,'--k')
    
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
ax.YLim = [1e-2 1]; ax.XLim = [-inf 500];
xticks([6 12 25 50 100 150 250 500])
yticks([0.01 0.1 0.25 0.5 1])
axis square
hold off
% legend(pl, fb_regList)
% grid on
xlabel('Spikes per Burst')
ylabel('Pr(X \geq x)')
%axis([4 55 0.1 1]), yticks([0.1,0.25,0.5,1])
saveas(gcf, './figs_HFS/fb-spnb-z1.png')

xlim([-inf 50]); ylim([1e-1 1])
%ax.YTick = [1e-2, 1e-1, 1];
%ax.XTick = [1e1, 1e2, 1e3];
% ax.FontSize = 16;
saveas(gcf, './figs_HFS/fb-spnb.png')
hold off
save('fb_spnb_points_plotted_hfs','points_plotted')
save('fb_spnb_fit_plotted_hfs','plotted_fit')
%% ANOVA

fit_table 
[~,~,~,stats] = aoctool(sp_feature,prob_vec,regLabel,0.05,'','','','off');
[fb_c,fb_means]=multcompare(stats,0.05,'off','','s');
fb_means

spnb_sp_feature_post_fb=sp_feature;
spnb_prob_vec_post_fb=prob_vec;
spnb_post_stim_string_vec=repmat("HFS",length(regLabel),1);
spnb_regLabel_string=string(regLabel);
spnb_regLabel_post_fb=categorical(join([spnb_regLabel_string,spnb_post_stim_string_vec]));
%title('Feedforward ISI Slope ANCOVA')
%xlabel('Slope')
saveas(gcf,'./figs_HFS/fb-isi-ancova_spnb.png')
save('HFS_fb_ancova_vars_spnb','spnb_sp_feature_post_fb','spnb_prob_vec_post_fb','spnb_regLabel_post_fb')
save('fb_means_HFS_spnb','fb_means')
%% FF error bar graph

figure( 'Position', [100 100 700 600])
xdata = [(1:4)-0.2, (1:4)+0.2];
ydata = [ff_means(:,1); fb_means(:,1)]';
errdata = [ff_means(:,2); fb_means(:,2)]';
b = bar(1:4,ydata(1:4));
b.FaceColor = 'flat';
for regi=1:4
    b.CData(regi,:) = col_array_rgb(regi,:);
end
hold on
errorbar(1:4, ydata(1:4), errdata(1:4), '.k')
xticks(1:4), xticklabels(ff_regList) 
ylabel 'Slope of spnb'
yticks(linspace(-5.5,0,12))
set(gca,'fontsize',17.5)
hold off
ylim([-5.5 0])
saveas(gcf,'./figs_HFS/ff-spnb-bar.png')
%% FB error bar graph

figure( 'Position', [100 100 700 600])
b = bar(1:4,ydata(5:end));
b.FaceColor = 'flat';
for regi=1:4
    b.CData(regi,:) = col_array_rgb(regi,:);
end
hold on
errorbar(1:4, ydata(5:end), errdata(5:end), '.k')
xticks(1:4), xticklabels(fb_regList) 
ylabel 'Slope of spnb'
yticks(linspace(-5.5,0,12))
set(gca,'fontsize',17.5)
hold off
ylim([-5.5 0])
%ytick_nums = [0.035, 0.05, 0.1, 0.25, 0.5, 1];
% yticks(log10(ytick_nums)); yticklabels(ytick_nums)
saveas(gcf,'./figs_HFS/fb-spnb-bar.png')

%% Burst duration histogram
% FF Histograms

disp 'Probabiblity distribution for Burst duration for FF axons'
moment_table = table([],[],[],[],[],[],[], 'VariableNames',{'regi','mean','median','mode','sd','se','emp_mode'});
% binEdge = logspace(0,4,100); binCenter = convert_edges_2_centers(binEdge); %11SD
binEdge = logspace(0,4,200); binCenter = convert_edges_2_centers(binEdge); %5SD
figure( 'Position', [100 100 700 600])
for regi = 1:4
    coli=1; histcount=[];
    subplot(2,2,figOrder(regi))
    for fi=1:10
        vec2plot = cell2mat(tunnel_spike_burst_dyn_table(tunnel_spike_burst_dyn_table.if_ff==1 & ... 
            tunnel_spike_burst_dyn_table.regi==regi & tunnel_spike_burst_dyn_table.fi==fi,:).BurstDuration);
        [hist_object,hist_plot,histcount(coli,:)] = log_binned_histogram(vec2plot, binEdge, "pdf",0);
        coli=coli+1;
    end
    histcount(~any(~isnan(histcount), 2),:)=[];
    meanC = mean(histcount,1); stdC = stdErr(histcount,1);
    
    % calculating mean and sd 
    f = fit(log(binCenter)', meanC', 'gauss1');
    mu = f.b1; sigma = f.c1/sqrt(2);
    moment_table.regi(regi)=regi;
    
    %fitted
    moment_table.mean(regi) =exp(mu+(sigma^2)/2);
    moment_table.mode(regi) = exp(mu-(sigma^2)/2);
    moment_table.sd(regi) = sqrt((exp(sigma^2)-1)*(moment_table.mean(regi) ^2));
    moment_table.median(regi)=exp(mu);
    n=length(cell2mat(tunnel_spike_burst_dyn_table(tunnel_spike_burst_dyn_table.if_ff==1 & ... 
            tunnel_spike_burst_dyn_table.regi==regi,:).BurstDuration));
    moment_table.n(regi)=n;
    moment_table.se(regi) = sqrt((exp(sigma^2)-1)*(moment_table.mean(regi) ^2))/sqrt(n);
       
    errorbar(binCenter, meanC, stdC,'k');
    title(ff_regList(regi))
    xlabel 'Burst Duration (ms)', ylabel 'Probability'
    set(gca, 'xscale','log', 'fontsize',18)
%     xlim([1e1 1e3]), ylim([0 0.2]) %11SD
    xlim([1e1 1e3]), ylim([0 0.05])
%     xticks([25 100 400])
    format_2x2_plot(regi)
        
    % Overlaying fit
    hold on
    fitvals = f(log(binCenter)');
    pf = plot(binCenter, fitvals);
    pf.LineStyle = '--';
    pf.LineWidth = 3;
    pf.Color = 'r';
    moment_table.emp_mode(regi) = binCenter(fitvals == max(fitvals));
%     moment_table.n(regi)=length(fitvals);
    moment_table.lnmu(regi)=mu;
    moment_table.lnsd(regi)=sigma;
    moment_table.n_fit(regi)=length(meanC);
end
saveas(gcf, './figs_HFS/ff-bd.png') 
moment_table
% ln FF testing fit

gof_table = table([],[],'variablenames',{'regi','rsquare'}); rowi=1;
figure( 'Position', [100 100 700 600])
for regi = 1:4
    coli=1; histcount=[];
    subplot(2,2,figOrder(regi))
    for fi=1:10
        vec2plot = cell2mat( tunnel_spike_burst_dyn_table(tunnel_spike_burst_dyn_table.if_ff==1 & ... 
            tunnel_spike_burst_dyn_table.regi==regi & tunnel_spike_burst_dyn_table.fi==fi,:).BurstDuration);
        [hist_object,hist_plot,histcount(coli,:)] = log_binned_histogram(vec2plot, binEdge, "pdf",0);
        coli=coli+1;
    end
    histcount(~any(~isnan(histcount), 2),:)=[];
    meanC = mean(histcount,1); stdC = stdErr(histcount,1);
    
    % calculating mean and sd 
    [f,gof] = fit(log(binCenter)', meanC', 'gauss1');
    mu = f.b1; sigma = f.c1/sqrt(2);       
    plot(f, log(binCenter), meanC)
    title(ff_regList(regi))
    xlabel 'ln(Burst Duration (ms) )'; ylabel 'Probability'
    legend off
    format_2x2_plot(regi)
    ylim([0 0.05])
    
    % populating GOF table
    gof_table.regi(rowi)=regi;
    gof_table.rsquare(rowi)=gof.rsquare;
    rowi=rowi+1;
end
saveas(gcf,'./figs_HFS/ff-bd-test-fit.png')
gof_table
save('HFS_ff_bd','moment_table','gof_table')
%% Burst duration histogram long
% FF Histograms

disp 'Probabiblity distribution for Burst duration for FF axons'
moment_table = table([],[],[],[],[], 'VariableNames',{'regi','mean','sd','mode','emp_mode'});
% binEdge = logspace(0,4,100); binCenter = convert_edges_2_centers(binEdge); %11SD
binEdge = logspace(0,4,200); binCenter = convert_edges_2_centers(binEdge); %5SD 
figure( 'Position', [100 100 2800 600])
tiledlayout(1,4, 'Padding', 'none', 'TileSpacing', 'compact');
for regi = 1:4
    nexttile
    coli=1; histcount=[];
%     subplot(2,2,figOrder(regi))
    for fi=1:10
        vec2plot = cell2mat(tunnel_spike_burst_dyn_table(tunnel_spike_burst_dyn_table.if_ff==1 & ... 
            tunnel_spike_burst_dyn_table.regi==regi & tunnel_spike_burst_dyn_table.fi==fi,:).BurstDuration);
        [hist_object,hist_plot,histcount(coli,:)] = log_binned_histogram(vec2plot, binEdge, "pdf",0);
        coli=coli+1;
    end
    histcount(~any(~isnan(histcount), 2),:)=[];
    meanC = mean(histcount,1); stdC = stdErr(histcount,1);
    
    % calculating mean and sd 
    f = fit(log(binCenter)', meanC', 'gauss1');
    mu = f.b1; sigma = f.c1/sqrt(2);
    moment_table.regi(regi)=regi;
    
    %fitted
    moment_table.mean(regi) =exp(mu+(sigma^2)/2);
    moment_table.mode(regi) = exp(mu-(sigma^2)/2);
    moment_table.sd(regi) = sqrt((exp(sigma^2)-1)*(moment_table.mean(regi) ^2));
       
    errorbar(binCenter, meanC, stdC,'k');
    title(ff_regList(regi))
    xlabel 'Burst Duration (ms)', ylabel 'Probability'
    set(gca, 'xscale','log', 'fontsize',18)
    %     xlim([1e1 1e3]), ylim([0 0.2]) %11SD
    xlim([1e1 1e3]), ylim([0 0.05])
    ax=gca;
    ax.FontSize=35;
%     xticks([25 100 400])
%     format_2x2_plot(regi)
        
    % Overlaying fit
    hold on
    fitvals = f(log(binCenter)');
    pf = plot(binCenter, fitvals);
    pf.LineStyle = '--';
    pf.LineWidth = 3;
    pf.Color = 'r';
    moment_table.emp_mode(regi) = binCenter(fitvals == max(fitvals));
    moment_table.n(regi)=length(fitvals);
end
saveas(gcf, './figs_HFS/ff-bd-long.png') 
moment_table
%% FB histograms

disp 'Probabiblity distribution for Burst duration for FB axons'
moment_table = table([],[],[],[],[],[],[], 'VariableNames',{'regi','mean','median','mode','sd','se','emp_mode'});
figure( 'Position', [100 100 700 600])
for regi = 1:4
    coli=1; histcount=[];
    subplot(2,2,figOrder(regi))
    for fi=1:10
        vec2plot = cell2mat(tunnel_spike_burst_dyn_table(tunnel_spike_burst_dyn_table.if_ff==0 & ... 
            tunnel_spike_burst_dyn_table.regi==regi & tunnel_spike_burst_dyn_table.fi==fi,:).BurstDuration);
        [hist_object,hist_plot,histcount(coli,:)] = log_binned_histogram(vec2plot, binEdge, "pdf",0);
        coli=coli+1;
    end
    histcount(~any(~isnan(histcount), 2),:)=[];
    meanC = mean(histcount,1); stdC = stdErr(histcount,1);
    
    % calculating mean and sd 
%     if regi~=3 %11SD
%         f = fit(log(binCenter)', meanC', 'gauss1');
%         mu = f.b1; sigma = f.c1/sqrt(2);
%         moment_table.regi(regi)=regi;
%     else
%         f = fit(log(binCenter)', meanC', 'gauss1', 'Exclude', find((ismember(meanC,maxk(meanC,3)))));
%         mu = f.b1; sigma = f.c1/sqrt(2);
%         moment_table.regi(regi)=regi;
%     end
    % 5SD
    if regi~=3
        f = fit(log(binCenter)', meanC', 'gauss1');
        mu = f.b1; sigma = f.c1/sqrt(2);
        moment_table.regi(regi)=regi;
    else
        f = fit(log(binCenter)', meanC', 'gauss1', 'Exclude', find((ismember(meanC,maxk(meanC,2)))));
        mu = f.b1; sigma = f.c1/sqrt(2);
        moment_table.regi(regi)=regi;
    end
    
    %fitted
    moment_table.mean(regi) =exp(mu+(sigma^2)/2);
    moment_table.mode(regi) = exp(mu-(sigma^2)/2);
    moment_table.sd(regi) = sqrt((exp(sigma^2)-1)*(moment_table.mean(regi) ^2));
    moment_table.median(regi)=exp(mu);
    n=length(cell2mat(tunnel_spike_burst_dyn_table(tunnel_spike_burst_dyn_table.if_ff==0 & ... 
            tunnel_spike_burst_dyn_table.regi==regi,:).BurstDuration));
    moment_table.n(regi)=n;
    moment_table.se(regi) = sqrt((exp(sigma^2)-1)*(moment_table.mean(regi) ^2))/sqrt(n);
    
    % Plotting data
    errorbar(binCenter, meanC, stdC,'k');
    title(fb_regList(regi))
    xlabel 'Burst Duration (ms)', ylabel 'Probability'
    set(gca, 'xscale','log', 'fontsize',18)
    %     xlim([1e1 1e3]), ylim([0 0.2]) %11SD
    xlim([1e1 1e3]), ylim([0 0.05])
    format_2x2_plot(regi)
        
    % Overlaying fit
    hold on
    fitvals = f(log(binCenter)');
    pf = plot(binCenter, fitvals);
    pf.LineStyle = '--';
    pf.LineWidth = 3;
    pf.Color = 'r';
    moment_table.emp_mode(regi) = binCenter(fitvals == max(fitvals));
%     moment_table.n(regi)=length(fitvals);
    moment_table.lnmu(regi)=mu;
    moment_table.lnsd(regi)=sigma;
    moment_table.n_fit(regi)=length(meanC);
end
moment_table
saveas(gcf, './figs_HFS/fb-bd.png')
%% FB testing fit

gof_table = table([],[],'variablenames',{'regi','rsquare'}); rowi=1;
figure( 'Position', [100 100 700 600])
%moment_table = table([],[],[],[], 'VariableNames',{'regi','mean','sd','mode'});
for regi = 1:4
    coli=1; histcount=[];
    subplot(2,2,figOrder(regi))
    for fi=1:10
        vec2plot = cell2mat( tunnel_spike_burst_dyn_table(tunnel_spike_burst_dyn_table.if_ff==0 & ... 
            tunnel_spike_burst_dyn_table.regi==regi & tunnel_spike_burst_dyn_table.fi==fi,:).BurstDuration);
        [hist_object,hist_plot,histcount(coli,:)] = log_binned_histogram(vec2plot, binEdge, "pdf",0);
        coli=coli+1;
    end
    histcount(~any(~isnan(histcount), 2),:)=[];
    meanC = mean(histcount,1); stdC = stdErr(histcount,1);
    
    % calculating mean and sd 
%     if regi~=3 %11SD
%         f = fit(log(binCenter)', meanC', 'gauss1');
%         mu = f.b1; sigma = f.c1/sqrt(2);
%         moment_table.regi(regi)=regi;
%     else
%         f = fit(log(binCenter)', meanC', 'gauss1', 'Exclude', find((ismember(meanC,maxk(meanC,3)))));
%         mu = f.b1; sigma = f.c1/sqrt(2);
%         moment_table.regi(regi)=regi;
%     end
    % 5SD
    if regi~=3 %11SD
        [f,gof] = fit(log(binCenter)', meanC', 'gauss1');
        mu = f.b1; sigma = f.c1/sqrt(2);
        moment_table.regi(regi)=regi;
    else
        [f,gof] = fit(log(binCenter)', meanC', 'gauss1', 'Exclude', find((ismember(meanC,maxk(meanC,3)))));
        mu = f.b1; sigma = f.c1/sqrt(2);
        moment_table.regi(regi)=regi;
    end
    
    %fitted
%     moment_table.mean(regi) =exp(mu+(sigma^2)/2);
%     moment_table.mode(regi) = exp(mu-(sigma^2)/2);
%     moment_table.sd(regi) = sqrt((exp(sigma^2)-1)*(moment_table.mean(regi) ^2));
    
    % calculating mean and sd 
%     f = fit(log(binCenter)', meanC', 'gauss1');
    mu = f.b1; sigma = f.c1/sqrt(2);       
    plot(f, log(binCenter), meanC)
    title(fb_regList(regi))
    xlabel 'ln(Burst Duration (ms) )'; ylabel 'Probability'
    legend off
    format_2x2_plot(regi)
    ylim([0 0.05])

    % populating GOF table
    gof_table.regi(rowi)=regi;
    gof_table.rsquare(rowi)=gof.rsquare;
    rowi=rowi+1;
end
saveas(gcf,'./figs_HFS/fb-bd-test-fit.png')
gof_table
save('HFS_fb_bd','moment_table','gof_table')
%% FB histograms long

disp 'Probabiblity distribution for Burst duration for FB axons'
moment_table = table([],[],[],[],[], 'VariableNames',{'regi','mean','sd','mode','emp_mode'});
figure( 'Position', [100 100 2800 600])
tiledlayout(1,4, 'Padding', 'none', 'TileSpacing', 'compact');
for regi = 1:4
    nexttile
    coli=1; histcount=[];
%     subplot(2,2,figOrder(regi))
    for fi=1:10
        vec2plot = cell2mat(tunnel_spike_burst_dyn_table(tunnel_spike_burst_dyn_table.if_ff==0 & ... 
            tunnel_spike_burst_dyn_table.regi==regi & tunnel_spike_burst_dyn_table.fi==fi,:).BurstDuration);
        [hist_object,hist_plot,histcount(coli,:)] = log_binned_histogram(vec2plot, binEdge, "pdf",0);
        coli=coli+1;
    end
    histcount(~any(~isnan(histcount), 2),:)=[];
    meanC = mean(histcount,1); stdC = stdErr(histcount,1);
    
    % calculating mean and sd 
%     if regi~=3 %11SD
%         f = fit(log(binCenter)', meanC', 'gauss1');
%         mu = f.b1; sigma = f.c1/sqrt(2);
%         moment_table.regi(regi)=regi;
%     else
%         f = fit(log(binCenter)', meanC', 'gauss1', 'Exclude', find((ismember(meanC,maxk(meanC,3)))));
%         mu = f.b1; sigma = f.c1/sqrt(2);
%         moment_table.regi(regi)=regi;
%     end
    % 5SD
    if regi~=3 
        f = fit(log(binCenter)', meanC', 'gauss1');
        mu = f.b1; sigma = f.c1/sqrt(2);
        moment_table.regi(regi)=regi;
    else
        f = fit(log(binCenter)', meanC', 'gauss1', 'Exclude', find((ismember(meanC,maxk(meanC,2)))));
        mu = f.b1; sigma = f.c1/sqrt(2);
        moment_table.regi(regi)=regi;
    end
    
    %fitted
    moment_table.mean(regi) =exp(mu+(sigma^2)/2);
    moment_table.mode(regi) = exp(mu-(sigma^2)/2);
    moment_table.sd(regi) = sqrt((exp(sigma^2)-1)*(moment_table.mean(regi) ^2));
    
    % Plotting data
    errorbar(binCenter, meanC, stdC,'k');
    title(fb_regList(regi))
    xlabel 'Burst Duration (ms)', ylabel 'Probability'
    set(gca, 'xscale','log', 'fontsize',18)
    %     xlim([1e1 1e3]), ylim([0 0.2]) %11SD
    xlim([1e1 1e3]), ylim([0 0.05])
    ax=gca;
    ax.FontSize=35;
%     format_2x2_plot(regi)
        
    % Overlaying fit
    hold on
    fitvals = f(log(binCenter)');
    pf = plot(binCenter, fitvals);
    pf.LineStyle = '--';
    pf.LineWidth = 3;
    pf.Color = 'r';
    moment_table.emp_mode(regi) = binCenter(fitvals == max(fitvals));
    moment_table.n(regi)=length(fitvals);
end
moment_table
saveas(gcf, './figs_HFS/fb-bd-long.png')
%% Intraburst spike rate
% FF Histograms

moment_table = table([],[],[],[],[],[],[], 'VariableNames',{'regi','mean','median','mode','sd','se','emp_mode'});
disp 'Probabiblity distribution for Intraburst spike rate for FF axons'
% binEdge = logspace(0,3,75); binCenter = convert_edges_2_centers(binEdge); %11SD
binEdge = logspace(0,3,150); binCenter = convert_edges_2_centers(binEdge); %5SD
figure( 'Position', [100 100 700 600])
for regi = 1:4
    coli=1; histcount=[];
    subplot(2,2,figOrder(regi))
    for fi=1:10
        vec2plot = cell2mat(tunnel_spike_burst_dyn_table(tunnel_spike_burst_dyn_table.if_ff==1 & ... 
            tunnel_spike_burst_dyn_table.regi==regi & tunnel_spike_burst_dyn_table.fi==fi,:).IntraBurstSpikeRate);
        [hist_object,hist_plot,histcount(coli,:)] = log_binned_histogram(vec2plot, binEdge, "pdf",0);
        coli=coli+1;
    end
    histcount(~any(~isnan(histcount), 2),:)=[];
    meanC = mean(histcount,1); stdC = stdErr(histcount,1);
    
    % calculating mean and sd 
%     f = fit(log(binCenter)', meanC', 'gauss1'); % 11SD\
    if regi==2 %5SD
        [~,max_idx]=sort(meanC','descend');
        f = fit(log(binCenter)', meanC', 'gauss1','Exclude',max_idx(1:3));
    else
        f = fit(log(binCenter)', meanC', 'gauss1');
    end
    mu = f.b1; sigma = f.c1/sqrt(2);
    moment_table.regi(regi)=regi;


    %fitted
    moment_table.mean(regi) =exp(mu+(sigma^2)/2);
    moment_table.mode(regi) = exp(mu-(sigma^2)/2);
    moment_table.sd(regi) = sqrt((exp(sigma^2)-1)*(moment_table.mean(regi) ^2));
    n=length(cell2mat(tunnel_spike_burst_dyn_table(tunnel_spike_burst_dyn_table.if_ff==1 & ... 
            tunnel_spike_burst_dyn_table.regi==regi,:).IntraBurstSpikeRate));
    moment_table.n(regi)=n;
    moment_table.se(regi) = sqrt((exp(sigma^2)-1)*(moment_table.mean(regi) ^2))/sqrt(n);   
    moment_table.median(regi)=exp(mu);
    % Plotting data
    errorbar(binCenter, meanC, stdC,'k');
    title(ff_regList(regi))
    xlabel 'Intraburst spike rate (Hz)', ylabel 'Probability'
    set(gca, 'xscale','log', 'fontsize',18,'XMinorTick','off')
    ylim([0 0.10]), xlim([20 620])
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
    moment_table.lnmu(regi)=mu;
    moment_table.lnsd(regi)=sigma;
    moment_table.n_fit(regi)=length(meanC);
end
saveas(gcf, './figs_HFS/ff-ibsr.png')

moment_table
% FF testing fit

gof_table = table([],[],'variablenames',{'regi','rsquare'}); rowi=1;
figure( 'Position', [100 100 700 600])
for regi = 1:4
     coli=1; histcount=[];
    subplot(2,2,figOrder(regi))
    for fi=1:10
        vec2plot = cell2mat(tunnel_spike_burst_dyn_table(tunnel_spike_burst_dyn_table.if_ff==1 & ... 
            tunnel_spike_burst_dyn_table.regi==regi & tunnel_spike_burst_dyn_table.fi==fi,:).IntraBurstSpikeRate);
        [hist_object,hist_plot,histcount(coli,:)] = log_binned_histogram(vec2plot, binEdge, "pdf",0);
        coli=coli+1;
    end
    histcount(~any(~isnan(histcount), 2),:)=[];
    meanC = mean(histcount,1); stdC = stdErr(histcount,1);
    
    % calculating mean and sd 
%     [f,gof] = fit(log(binCenter)', meanC', 'gauss1'); % 11SD
    if regi==2 %5SD
        [~,max_idx]=sort(meanC','descend');
        [f,gof] = fit(log(binCenter)', meanC', 'gauss1','Exclude',max_idx(1:3));
    else
        [f,gof] = fit(log(binCenter)', meanC', 'gauss1');
    end
    mu = f.b1; sigma = f.c1/sqrt(2);       
    plot(f, log(binCenter), meanC)
    title(ff_regList(regi))
    xlabel 'ln(Intraburst spike rate (Hz))', ylabel 'Probability'
    legend off
    ylim([0 0.10])
    format_2x2_plot(regi)

    % populating GOF table
    gof_table.regi(rowi)=regi;
    gof_table.rsquare(rowi)=gof.rsquare;
    rowi=rowi+1;
end
saveas(gcf,'./figs_HFS/ff-ibsr-test-fit.png')
gof_table
save('HFS_ff_IBSR','moment_table','gof_table')
%% Intraburst spike rate Long
% FF Histograms

saveMean=[];
moment_table = table([],[],[],[],[], 'VariableNames',{'regi','mean','sd','mode','emp_mode'});
disp 'Probabiblity distribution for Intraburst spike rate for FF axons'
% binEdge = logspace(0,3,75); binCenter = convert_edges_2_centers(binEdge); %11SD
binEdge = logspace(0,3,150); binCenter = convert_edges_2_centers(binEdge); %5SD 
figure( 'Position', [100 100 2800 600])
tiledlayout(1,4, 'Padding', 'none', 'TileSpacing', 'compact');
for regi = 1:4
    nexttile
    coli=1; histcount=[];
%     subplot(2,2,figOrder(regi))
    for fi=1:10
        vec2plot = cell2mat(tunnel_spike_burst_dyn_table(tunnel_spike_burst_dyn_table.if_ff==1 & ... 
            tunnel_spike_burst_dyn_table.regi==regi & tunnel_spike_burst_dyn_table.fi==fi,:).IntraBurstSpikeRate);
        [hist_object,hist_plot,histcount(coli,:)] = log_binned_histogram(vec2plot, binEdge, "pdf",0);
        coli=coli+1;
    end
    histcount(~any(~isnan(histcount), 2),:)=[];
    meanC = mean(histcount,1); stdC = stdErr(histcount,1);
    
    % calculating mean and sd 
%     f = fit(log(binCenter)', meanC', 'gauss1'); %11SD
    if regi==2 %5SD
        [~,max_idx]=sort(meanC','descend');
        f = fit(log(binCenter)', meanC', 'gauss1','Exclude',max_idx(1:3));
    else
        f = fit(log(binCenter)', meanC', 'gauss1');
    end
    mu = f.b1; sigma = f.c1/sqrt(2);
    moment_table.regi(regi)=regi;
    %saveMean{regi}=number_samples_at_index_211116(log(binCenter),histcount);
    
    %fitted
    moment_table.mean(regi) =exp(mu+(sigma^2)/2);
    moment_table.mode(regi) = exp(mu-(sigma^2)/2);
    moment_table.sd(regi) = sqrt((exp(sigma^2)-1)*(moment_table.mean(regi) ^2));
    n=length(cell2mat(tunnel_spike_burst_dyn_table(tunnel_spike_burst_dyn_table.if_ff==1 & ... 
            tunnel_spike_burst_dyn_table.regi==regi,:).IntraBurstSpikeRate));
    moment_table.n(regi)=n;
    moment_table.se(regi) = sqrt((exp(sigma^2)-1)*(moment_table.mean(regi) ^2))/sqrt(n);
       
    % Plotting data
    errorbar(binCenter, meanC, stdC,'k');
    title(ff_regList(regi))
    xlabel 'Intraburst spike rate (Hz)', ylabel 'Probability'
    set(gca, 'xscale','log', 'fontsize',18,'XMinorTick','off')
    ylim([0 0.10]), xlim([20 620])
    xticks([25 50 100 200 400])
    ax=gca;
    ax.FontSize=30;
%     format_2x2_plot(regi)
        
    %ylim([0,1])
    %continue
    % Overlaying fit
    hold on
    fitvals = f(log(binCenter)');
    pf = plot(binCenter, fitvals);
    pf.LineStyle = '--';
    pf.LineWidth = 3;
    pf.Color = 'r';
    moment_table.emp_mode(regi) = binCenter(fitvals == max(fitvals));
    moment_table.n(regi)=length(fitvals);
end
saveas(gcf, './figs_HFS/ff-ibsr-long.png')

moment_table
%% FB Histograms

disp 'Probabiblity distribution for Intraburst spike rate for FB axons'
moment_table = table([],[],[],[],[],[],[], 'VariableNames',{'regi','mean','median','mode','sd','se','emp_mode'});
figure( 'Position', [100 100 700 600])
for regi = 1:4
    coli=1; histcount=[];
    subplot(2,2,figOrder(regi))
    for fi=1:10
        vec2plot = cell2mat(tunnel_spike_burst_dyn_table(tunnel_spike_burst_dyn_table.if_ff==0 & ... 
            tunnel_spike_burst_dyn_table.regi==regi & tunnel_spike_burst_dyn_table.fi==fi,:).IntraBurstSpikeRate);
        [hist_object,hist_plot,histcount(coli,:)] = log_binned_histogram(vec2plot, binEdge, "pdf",0);
        coli=coli+1;
    end
    histcount(~any(~isnan(histcount), 2),:)=[];
    meanC = mean(histcount,1); stdC = stdErr(histcount,1);
    
    % calculating mean and sd 
%     f = fit(log(binCenter)', meanC', 'gauss1'); %11SD
    if regi==3 %5SD
        [~,max_idx]=sort(meanC','descend');
        f = fit(log(binCenter)', meanC', 'gauss1','Exclude',max_idx(1:2));
    else
        f = fit(log(binCenter)', meanC', 'gauss1');
    end
    mu = f.b1; sigma = f.c1/sqrt(2);
    moment_table.regi(regi)=regi;
    
    %fitted
    moment_table.mean(regi) =exp(mu+(sigma^2)/2);
    moment_table.mode(regi) = exp(mu-(sigma^2)/2);
    moment_table.sd(regi) = sqrt((exp(sigma^2)-1)*(moment_table.mean(regi) ^2));
    n=length(cell2mat(tunnel_spike_burst_dyn_table(tunnel_spike_burst_dyn_table.if_ff==0 & ... 
            tunnel_spike_burst_dyn_table.regi==regi,:).IntraBurstSpikeRate));
    moment_table.n(regi)=n;
    moment_table.se(regi) = sqrt((exp(sigma^2)-1)*(moment_table.mean(regi) ^2))/sqrt(n);
    moment_table.median(regi)=exp(mu);
       
    % Plotting raw
    errorbar(binCenter, meanC, stdC,'k');
    title(fb_regList(regi))
    xlabel 'Intraburst spike rate (Hz)', ylabel 'Probability'
    set(gca, 'xscale','log', 'fontsize',18,'XMinorTick','off')
    ylim([0 0.10]), xlim([20 620])
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
%     axis square
    moment_table.lnmu(regi)=mu;
    moment_table.lnsd(regi)=sigma;
    moment_table.n_fit(regi)=length(meanC);
end
saveas(gcf, './figs_HFS/fb-ibsr.png')

moment_table
%% FB testing fit

gof_table = table([],[],'variablenames',{'regi','rsquare'}); rowi=1;
figure( 'Position', [100 100 700 600])
for regi = 1:4
    coli=1;histcount=[];
    subplot(2,2,figOrder(regi))
    for fi=1:10
        vec2plot = cell2mat( tunnel_spike_burst_dyn_table(tunnel_spike_burst_dyn_table.if_ff==0 &... 
            tunnel_spike_burst_dyn_table.regi==regi & tunnel_spike_burst_dyn_table.fi==fi,:).IntraBurstSpikeRate);
        [hist_object,hist_plot,histcount(coli,:)] = log_binned_histogram(vec2plot, binEdge, "pdf",0);
        coli=coli+1;
    end
    histcount(~any(~isnan(histcount), 2),:)=[];
    meanC = mean(histcount,1); stdC = stdErr(histcount,1);
    
    % calculating mean and sd 
%     [f,gof] = fit(log(binCenter)', meanC', 'gauss1'); % 11SD
    if regi==3 %5SD
        [~,max_idx]=sort(meanC','descend');
        [f,gof] = fit(log(binCenter)', meanC', 'gauss1','Exclude',max_idx(1:2));
    else
        [f,gof] = fit(log(binCenter)', meanC', 'gauss1');
    end
    mu = f.b1; sigma = f.c1/sqrt(2);       
    plot(f, log(binCenter), meanC)
    title(fb_regList(regi))
    ylim([0 0.10])
    xlabel 'ln(Intraburst spike rate (Hz))', ylabel 'Probability'
    legend off
    format_2x2_plot(regi)
    
    % populating GOF table
    gof_table.regi(rowi)=regi;
    gof_table.rsquare(rowi)=gof.rsquare;
    rowi=rowi+1;
end
saveas(gcf,'./figs_HFS/fb-ibsr-test-fit.png')
gof_table
save('HFS_fb_IBSR','moment_table','gof_table')
%% FB Histograms long

disp 'Probabiblity distribution for Intraburst spike rate for FB axons'
moment_table = table([],[],[],[],[], 'VariableNames',{'regi','mean','sd','mode','emp_mode'});
figure( 'Position', [100 100 2800 600])
tiledlayout(1,4, 'Padding', 'none', 'TileSpacing', 'compact');
for regi = 1:4
    nexttile
    coli=1; histcount=[];
%     subplot(2,2,figOrder(regi))
    for fi=1:10
        vec2plot = cell2mat(tunnel_spike_burst_dyn_table(tunnel_spike_burst_dyn_table.if_ff==0 & ... 
            tunnel_spike_burst_dyn_table.regi==regi & tunnel_spike_burst_dyn_table.fi==fi,:).IntraBurstSpikeRate);
        [hist_object,hist_plot,histcount(coli,:)] = log_binned_histogram(vec2plot, binEdge, "pdf",0);
        coli=coli+1;
    end
    histcount(~any(~isnan(histcount), 2),:)=[];
    meanC = mean(histcount,1); stdC = stdErr(histcount,1);
    
    % calculating mean and sd 
%     f = fit(log(binCenter)', meanC', 'gauss1'); %11SD
    if regi==3 %5SD
        [~,max_idx]=sort(meanC','descend');
        f = fit(log(binCenter)', meanC', 'gauss1','Exclude',max_idx(1:2));
    else
        f = fit(log(binCenter)', meanC', 'gauss1');
    end
    mu = f.b1; sigma = f.c1/sqrt(2);
    moment_table.regi(regi)=regi;
    
    %fitted
    moment_table.mean(regi) =exp(mu+(sigma^2)/2);
    moment_table.mode(regi) = exp(mu-(sigma^2)/2);
    moment_table.sd(regi) = sqrt((exp(sigma^2)-1)*(moment_table.mean(regi) ^2));
    n=length(cell2mat(tunnel_spike_burst_dyn_table(tunnel_spike_burst_dyn_table.if_ff==0 & ... 
            tunnel_spike_burst_dyn_table.regi==regi,:).IntraBurstSpikeRate));
    moment_table.n(regi)=n;
    moment_table.se(regi) = sqrt((exp(sigma^2)-1)*(moment_table.mean(regi) ^2))/sqrt(n);
    
    % Plotting raw
    errorbar(binCenter, meanC, stdC,'k');
    title(fb_regList(regi))
    xlabel 'Intraburst spike rate (Hz)', ylabel 'Probability'
    set(gca, 'xscale','log', 'fontsize',18,'XMinorTick','off')
    ylim([0 0.10]), xlim([20 620])
    xticks([25 50 100 200 400])
    ax=gca;
    ax.FontSize=30;
%     format_2x2_plot(regi)
        
    % Overlaying fit
    hold on
    fitvals = f(log(binCenter)');
    pf = plot(binCenter, fitvals);
    pf.LineStyle = '--';
    pf.LineWidth = 3;
    pf.Color = 'r';
    moment_table.emp_mode(regi) = binCenter(fitvals == max(fitvals));
    moment_table.n(regi)=length(fitvals);
%     axis square
end
saveas(gcf, './figs_HFS/fb-ibsr-long.png')