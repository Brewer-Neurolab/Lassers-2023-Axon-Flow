%% IBSR and BD ANOVA from Inputs

clear
clc
close all
cd("D:\Brewer lab data\HFS\")
% 11 SD
% nostim=load('D:\Brewer lab data\HFS\No Stim\23-Nov-2021_A\spike_burst_dyn_table_stim.mat');
% % hfs5=load('D:\Brewer lab data\HFS\Theta Stim\03-Mar-2022_B\spike_burst_dyn_table_stim.mat');
% hfs5=load('D:\Brewer lab data\HFS\Theta Stim\10-May-2022_A\spike_burst_dyn_table_stim.mat');
% hfs40=load('D:\Brewer lab data\HFS\HFS Stim\24-Nov-2021_A\spike_burst_dyn_table_stim.mat');

% 5SD
nostim=load('D:\Brewer lab data\HFS\No Stim\02-Oct-2022_A\spike_burst_dyn_table_stim.mat');
hfs5=load('D:\Brewer lab data\HFS\Theta Stim\27-Sep-2022_A\spike_burst_dyn_table_stim.mat');
hfs40=load('D:\Brewer lab data\HFS\HFS Stim\03-Oct-2022_A\spike_burst_dyn_table_stim.mat');

%% IBSR

%% Feedforward

regions=["EC-DG","DG-CA3","CA3-CA1","CA1-EC"];
stims=["No Stim", "5HFS", "40HFS"];

IBSR_ff_p=[];

binEdge = logspace(0,3,75);

for i=1:length(regions)
    
    nostim_ff_ibsr=cell2mat(nostim.spike_burst_dyn_table_stim.IntraBurstSpikeRate...
        (nostim.spike_burst_dyn_table_stim.regi==i & nostim.spike_burst_dyn_table_stim.if_ff==1));
    hfs5_ff_ibsr=cell2mat(hfs5.spike_burst_dyn_table_stim.IntraBurstSpikeRate...
        (hfs5.spike_burst_dyn_table_stim.regi==i & hfs5.spike_burst_dyn_table_stim.if_ff==1));
    hfs40_ff_ibsr=cell2mat(hfs40.spike_burst_dyn_table_stim.IntraBurstSpikeRate...
        (hfs40.spike_burst_dyn_table_stim.regi==i & hfs40.spike_burst_dyn_table_stim.if_ff==1));
    
    figure
    %     histogram(log10(nostim_ff_ibsr),100)
    histogram((nostim_ff_ibsr),(binEdge))
    ylabel("Counts")
    xlabel("IBSR (Hz)")
    ax=gca;
    ax.XScale="log";
    ax.FontSize=18;
    title(strcat("Feedforward IBSR ",regions(i)))
    saveas(gcf,strcat(".\three stim share figs\IBSR_Hist_ff_NoStim_",regions(i),".png"))
    
%     histogram(log10(hfs5_ff_ibsr),100)
    histogram((hfs5_ff_ibsr),(binEdge))
    ylabel("Counts")
    xlabel("IBSR (Hz)")
    ax=gca;
    ax.XScale="log";
    ax.FontSize=18;
    title(strcat("Feedforward IBSR ",regions(i)))
    saveas(gcf,strcat(".\three stim share figs\IBSR_Hist_ff_5HFS_",regions(i),".png"))
    
%     histogram(log10(hfs40_ff_ibsr),100)
    histogram((hfs40_ff_ibsr),(binEdge))
    ylabel("Counts")
    xlabel("IBSR (Hz)")
    ax=gca;
    ax.XScale="log";
    ax.FontSize=18;
    title(strcat("Feedforward IBSR ",regions(i)))
    saveas(gcf,strcat(".\three stim share figs\IBSR_Hist_ff_40HFS_",regions(i),".png"))
    
    % Log Space bins
    
    hist_vec=log10([nostim_ff_ibsr;hfs5_ff_ibsr;hfs40_ff_ibsr]);
    stim_group=[repmat(stims(1),length(nostim_ff_ibsr),1);repmat(stims(2),length(hfs5_ff_ibsr),1);...
        repmat(stims(3),length(hfs40_ff_ibsr),1)];
    
    figure
%     histogram(log10(nostim_ff_ibsr),100)
    histogram(log10(nostim_ff_ibsr),log10(binEdge))
    ylabel("Counts")
    xlabel("Log(IBSR)")
    title(strcat("Feedforward IBSR ",regions(i)))
    saveas(gcf,strcat(".\three stim share figs\log_IBSR_Hist_ff_NoStim_",regions(i),".png"))
    
%     histogram(log10(hfs5_ff_ibsr),100)
    histogram(log10(hfs5_ff_ibsr),log10(binEdge))
    ylabel("Counts")
    xlabel("Log(IBSR)")
    title(strcat("Feedforward IBSR ",regions(i)))
    saveas(gcf,strcat(".\three stim share figs\log_IBSR_Hist_ff_5HFS_",regions(i),".png"))
    
%     histogram(log10(hfs40_ff_ibsr),100)
    histogram(log10(hfs40_ff_ibsr),log10(binEdge))
    ylabel("Counts")
    xlabel("Log(IBSR)")
    title(strcat("Feedforward IBSR ",regions(i)))
    saveas(gcf,strcat(".\three stim share figs\log_IBSR_Hist_ff_40HFS_",regions(i),".png"))
    
    
%     figure
    [~,~,stats]=anova1(hist_vec,stim_group);
    saveas(gcf,strcat(".\three stim share figs\IBSR_Box_ff_",regions(i),".png"))
    figure
    [ff_c,ff_means]=multcompare(stats);
    
    diff1=((10^ff_means(2,1)-10^ff_means(1,1))/10^ff_means(1,1))*100;
    diff2=((10^ff_means(3,1)-10^ff_means(1,1))/10^ff_means(1,1))*100;
    diff3=((10^ff_means(3,1)-10^ff_means(2,1))/10^ff_means(2,1))*100;
    IBSR_ff_diff{i}=[diff1;diff2;diff3];
    
    IBSR_ff_p{i}=ff_c;
    IBSR_ff_m{i}=ff_means;
    
    
    figure
    hold on
    
    b=bar(categorical(stims,stims),ff_means(:,1),'FaceColor','flat');
    b.CData=[[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250]];
    e=errorbar(ff_means(:,1),ff_means(:,2),'k');
    e.LineStyle='none';
    hold off
    
    ylim([1.5,3])
    ylabel("Log(IBSR)")
    xlabel("Stimulation")
    title("Feedforward IBSR "+regions(i))
    
    saveas(gcf,strcat(".\three stim share figs\IBSR_Bar_ff_",regions(i),".png"))
    
    %close all
end

figure
hold on

means=[]; err=[];
for i=1:length(IBSR_ff_m)
    means=[means;IBSR_ff_m{i}(:,1)'];
    err=[err;IBSR_ff_m{i}(:,2)'];
end

b=bar(categorical(regions,regions),means,'FaceColor','flat');
% b.CData=[[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250]];

[ngroups,nbars] = size(means);
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end

e=errorbar(x',means,err,'k','LineStyle','none');
hold off

ylim([1.5,2.5])
ylabel("Log(IBSR)")
xlabel("Subregion")
title("Feedforward IBSR ")

saveas(gcf,strcat(".\three stim share figs\IBSR_Bar_ff.png"))
%% Feedback

regions=["DG-EC","CA3-DG","CA1-CA3","EC-CA1"];
stims=["No Stim", "5HFS", "40HFS"];

IBSR_fb_p=[];

% binEdge = logspace(0,4,100);

for i=1:length(regions)
    nostim_fb_ibsr=cell2mat(nostim.spike_burst_dyn_table_stim.IntraBurstSpikeRate...
        (nostim.spike_burst_dyn_table_stim.regi==i & nostim.spike_burst_dyn_table_stim.if_ff==0));
    hfs5_fb_ibsr=cell2mat(hfs5.spike_burst_dyn_table_stim.IntraBurstSpikeRate...
        (hfs5.spike_burst_dyn_table_stim.regi==i & hfs5.spike_burst_dyn_table_stim.if_ff==0));
    hfs40_fb_ibsr=cell2mat(hfs40.spike_burst_dyn_table_stim.IntraBurstSpikeRate...
        (hfs40.spike_burst_dyn_table_stim.regi==i & hfs40.spike_burst_dyn_table_stim.if_ff==0));
    
    figure
    
    histogram((nostim_fb_ibsr),(binEdge))
%     histogram(log10(nostim_fb_ibsr),100)
    ylabel("Counts")
    xlabel("IBSR (Hz)")
    ax=gca;
    ax.XScale="log";
    ax.FontSize=18;
    title(strcat("Feedback IBSR ",regions(i)))
    saveas(gcf,strcat(".\three stim share figs\IBSR_Hist_fb_NoStim_",regions(i),".png"))
    
    histogram((hfs5_fb_ibsr),(binEdge))
%     histogram(log10(hfs5_fb_ibsr),100)
    ylabel("Counts")
    xlabel("IBSR (Hz)")
    ax=gca;
    ax.XScale="log";
    ax.FontSize=18;
    title(strcat("Feedback IBSR ",regions(i)))
    saveas(gcf,strcat(".\three stim share figs\IBSR_Hist_fb_5HFS_",regions(i),".png"))
    
    histogram((hfs40_fb_ibsr),(binEdge))
%     histogram(log10(hfs40_fb_ibsr),100)
    ylabel("Counts")
    xlabel("IBSR (Hz)")
    ax=gca;
    ax.XScale="log";
    ax.FontSize=18;
    title(strcat("Feedback IBSR ",regions(i)))
    saveas(gcf,strcat(".\three stim share figs\IBSR_Hist_fb_40HFS_",regions(i),".png"))
    
    % log space bins
    
    hist_vec=log10([nostim_fb_ibsr;hfs5_fb_ibsr;hfs40_fb_ibsr]);
    stim_group=[repmat(stims(1),length(nostim_fb_ibsr),1);repmat(stims(2),length(hfs5_fb_ibsr),1);...
        repmat(stims(3),length(hfs40_fb_ibsr),1)];
    
    figure
    
    histogram(log10(nostim_fb_ibsr),log10(binEdge))
%     histogram(log10(nostim_fb_ibsr),100)
    ylabel("Counts")
    xlabel("Log(IBSR)")
    title(strcat("Feedback IBSR ",regions(i)))
    saveas(gcf,strcat(".\three stim share figs\log_IBSR_Hist_fb_NoStim_",regions(i),".png"))
    
    histogram(log10(hfs5_fb_ibsr),log10(binEdge))
%     histogram(log10(hfs5_fb_ibsr),100)
    ylabel("Counts")
    xlabel("Log(IBSR)")
    title(strcat("Feedback IBSR ",regions(i)))
    saveas(gcf,strcat(".\three stim share figs\log_IBSR_Hist_fb_5HFS_",regions(i),".png"))
    
    histogram(log10(hfs40_fb_ibsr),log10(binEdge))
%     histogram(log10(hfs40_fb_ibsr),100)
    ylabel("Counts")
    xlabel("Log(IBSR)")
    title(strcat("Feedback IBSR ",regions(i)))
    saveas(gcf,strcat(".\three stim share figs\log_IBSR_Hist_fb_40HFS_",regions(i),".png"))
    
%     figure
    [~,~,stats]=anova1(hist_vec,stim_group);
    saveas(gcf,strcat(".\three stim share figs\IBSR_Box_fb_",regions(i),".png"))
    figure
    [ff_c,ff_means]=multcompare(stats);
    
    diff1=((10^ff_means(2,1)-10^ff_means(1,1))/10^ff_means(1,1))*100;
    diff2=((10^ff_means(3,1)-10^ff_means(1,1))/10^ff_means(1,1))*100;
    diff3=((10^ff_means(3,1)-10^ff_means(2,1))/10^ff_means(2,1))*100;
    IBSR_fb_diff{i}=[diff1;diff2;diff3];
    
    IBSR_fb_p{i}=ff_c;
    IBSR_fb_m{i}=ff_means;
    
    figure
    hold on
    
    b=bar(categorical(stims,stims),ff_means(:,1),'FaceColor','flat');
    b.CData=[[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250]];
    e=errorbar(ff_means(:,1),ff_means(:,2),'k');
    e.LineStyle='none';
    hold off
    
    ylim([1.5,3])
    ylabel("Log(IBSR)")
    xlabel("Stimulation")
    title("Feedback IBSR "+regions(i))
    
    saveas(gcf,strcat(".\three stim share figs\IBSR_Bar_fb_",regions(i),".png"))
    
    %close all
end
figure
hold on

means=[]; err=[];
for i=1:length(IBSR_fb_m)
    means=[means;IBSR_fb_m{i}(:,1)'];
    err=[err;IBSR_fb_m{i}(:,2)'];
end

b=bar(categorical(regions,regions),means,'FaceColor','flat');
% b.CData=[[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250]];

[ngroups,nbars] = size(means);
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end

e=errorbar(x',means,err,'k','LineStyle','none');
hold off

ylim([1.5,2.5])
ylabel("Log(IBSR)")
xlabel("Subregion")
title("Feedback IBSR ")

saveas(gcf,strcat(".\three stim share figs\IBSR_Bar_fb.png"))

%% IBSR Well
cd("D:\Brewer lab data\HFS\")
nostim=load('D:\Brewer lab data\HFS\No Stim\Wells_5SD_500maxSD\well_spike_dynamics_table_nostim.mat');
hfs5=load('D:\Brewer lab data\HFS\Theta Stim\Wells_5SD_500maxSD\well_spike_dynamics_table_theta.mat');
hfs40=load('D:\Brewer lab data\HFS\HFS Stim\Wells_5SD_500maxSD\well_spike_dynamics_table_hfs.mat');

regions=["EC","DG","CA3","CA1"];
stims=["No Stim", "5HFS", "40HFS"];

IBSR_well_p=[];

binEdge = logspace(0,4,150);

for i=1:length(regions)
    nostim_fb_ibsr=cell2mat(nostim.well_spike_dynamics_table.IntraBurstSpikeRate...
        (nostim.well_spike_dynamics_table.regi==i));
    hfs5_fb_ibsr=cell2mat(hfs5.well_spike_dynamics_table.IntraBurstSpikeRate...
        (hfs5.well_spike_dynamics_table.regi==i));
    hfs40_fb_ibsr=cell2mat(hfs40.well_spike_dynamics_table.IntraBurstSpikeRate...
        (hfs40.well_spike_dynamics_table.regi==i));
    
    figure
    
    histogram((nostim_fb_ibsr),(binEdge))
%     histogram(log10(nostim_fb_ibsr),100)
    ylabel("Counts")
    xlabel("IBSR (Hz)")
    ax=gca;
    ax.XScale="log";
    ax.FontSize=18;
    title(strcat("Well IBSR ",regions(i)))
    saveas(gcf,strcat(".\three stim share figs\IBSR_Hist_well_NoStim_",regions(i),".png"))
    
    histogram((hfs5_fb_ibsr),(binEdge))
%     histogram(log10(hfs5_fb_ibsr),100)
    ylabel("Counts")
    xlabel("IBSR (Hz)")
    ax=gca;
    ax.XScale="log";
    ax.FontSize=18;
    title(strcat("Well IBSR ",regions(i)))
    saveas(gcf,strcat(".\three stim share figs\IBSR_Hist_well_5HFS_",regions(i),".png"))
    
    histogram((hfs40_fb_ibsr),(binEdge))
%     histogram(log10(hfs40_fb_ibsr),100)
    ylabel("Counts")
    xlabel("IBSR (Hz)")
    ax=gca;
    ax.XScale="log";
    ax.FontSize=18;
    title(strcat("Well IBSR ",regions(i)))
    saveas(gcf,strcat(".\three stim share figs\IBSR_Hist_well_40HFS_",regions(i),".png"))
    
    % log bins
    
    hist_vec=log10([nostim_fb_ibsr;hfs5_fb_ibsr;hfs40_fb_ibsr]);
    stim_group=[repmat(stims(1),length(nostim_fb_ibsr),1);repmat(stims(2),length(hfs5_fb_ibsr),1);...
        repmat(stims(3),length(hfs40_fb_ibsr),1)];
    
    figure
    
    histogram(log10(nostim_fb_ibsr),log10(binEdge))
%     histogram(log10(nostim_fb_ibsr),100)
    ylabel("Counts")
    xlabel("Log(IBSR)")
    title(strcat("Well IBSR ",regions(i)))
    saveas(gcf,strcat(".\three stim share figs\log_IBSR_Hist_well_NoStim_",regions(i),".png"))
    
    histogram(log10(hfs5_fb_ibsr),log10(binEdge))
%     histogram(log10(hfs5_fb_ibsr),100)
    ylabel("Counts")
    xlabel("Log(IBSR)")
    title(strcat("Well IBSR ",regions(i)))
    saveas(gcf,strcat(".\three stim share figs\log_IBSR_Hist_well_5HFS_",regions(i),".png"))
    
    histogram(log10(hfs40_fb_ibsr),log10(binEdge))
%     histogram(log10(hfs40_fb_ibsr),100)
    ylabel("Counts")
    xlabel("Log(IBSR)")
    title(strcat("Well IBSR ",regions(i)))
    saveas(gcf,strcat(".\three stim share figs\log_IBSR_Hist_well_40HFS_",regions(i),".png"))
    
%     figure
    [~,~,stats]=anova1(hist_vec,stim_group);
    saveas(gcf,strcat(".\three stim share figs\IBSR_Box_well_",regions(i),".png"))
    figure
    [ff_c,ff_means]=multcompare(stats);
    
    diff1=((10^ff_means(2,1)-10^ff_means(1,1))/10^ff_means(1,1))*100;
    diff2=((10^ff_means(3,1)-10^ff_means(1,1))/10^ff_means(1,1))*100;
    diff3=((10^ff_means(3,1)-10^ff_means(2,1))/10^ff_means(2,1))*100;
    IBSR_well_diff{i}=[diff1;diff2;diff3];
    
    IBSR_well_p{i}=ff_c;
    IBSR_well_m{i}=ff_means;
    
    figure
    hold on
    
    b=bar(categorical(stims,stims),ff_means(:,1),'FaceColor','flat');
    b.CData=[[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250]];
    e=errorbar(ff_means(:,1),ff_means(:,2),'k');
    e.LineStyle='none';
    hold off
    
    ylim([1.5,3])
    ylabel("Log(IBSR)")
    xlabel("Stimulation")
    title("Well IBSR "+regions(i))
    
    saveas(gcf,strcat(".\three stim share figs\IBSR_Bar_Well_",regions(i),".png"))
    
    %close all
end
figure
hold on

means=[]; err=[];
for i=1:length(IBSR_well_m)
    means=[means;IBSR_well_m{i}(:,1)'];
    err=[err;IBSR_well_m{i}(:,2)'];
end

b=bar(categorical(regions,regions),means,'FaceColor','flat');
% b.CData=[[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250]];

[ngroups,nbars] = size(means);
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end

e=errorbar(x',means,err,'k','LineStyle','none');
hold off

ylim([1.5,2.5])
ylabel("Log(IBSR)")
xlabel("Subregion")
title("Well IBSR ")

saveas(gcf,strcat(".\three stim share figs\IBSR_Bar_well.png"))

%% BD
% 11SD
cd("D:\Brewer lab data\HFS\")
% nostim=load('D:\Brewer lab data\HFS\No Stim\23-Nov-2021_A\spike_burst_dyn_table_stim.mat');
% % hfs5=load('D:\Brewer lab data\HFS\Theta Stim\03-Mar-2022_B\spike_burst_dyn_table_stim.mat');
% hfs5=load('D:\Brewer lab data\HFS\Theta Stim\10-May-2022_A\spike_burst_dyn_table_stim.mat');
% hfs40=load('D:\Brewer lab data\HFS\HFS Stim\24-Nov-2021_A\spike_burst_dyn_table_stim.mat');

%5SD
nostim=load('D:\Brewer lab data\HFS\No Stim\02-Oct-2022_A\spike_burst_dyn_table_stim.mat');
hfs5=load('D:\Brewer lab data\HFS\Theta Stim\27-Sep-2022_A\spike_burst_dyn_table_stim.mat');
hfs40=load('D:\Brewer lab data\HFS\HFS Stim\03-Oct-2022_A\spike_burst_dyn_table_stim.mat');
%% Feedforward

regions=["EC-DG","DG-CA3","CA3-CA1","CA1-EC"];
stims=["No Stim", "5HFS", "40HFS"];

BD_ff_p=[];

binEdge = logspace(0,4,100);

for i=1:length(regions)
    nostim_ff_bd=cell2mat(nostim.spike_burst_dyn_table_stim.BurstDuration...
        (nostim.spike_burst_dyn_table_stim.regi==i & nostim.spike_burst_dyn_table_stim.if_ff==1));
    hfs5_ff_bd=cell2mat(hfs5.spike_burst_dyn_table_stim.BurstDuration...
        (hfs5.spike_burst_dyn_table_stim.regi==i & hfs5.spike_burst_dyn_table_stim.if_ff==1));
    hfs40_ff_bd=cell2mat(hfs40.spike_burst_dyn_table_stim.BurstDuration...
        (hfs40.spike_burst_dyn_table_stim.regi==i & hfs40.spike_burst_dyn_table_stim.if_ff==1));
    
    figure
    
    histogram((nostim_ff_bd),(binEdge))
%     histogram(log10(nostim_ff_bd),100)
    ylabel("Counts")
    xlabel("BD (ms)")
    ax=gca;
    ax.XScale="log";
    ax.FontSize=18;
    title(strcat("Feedforward BD ",regions(i)))
    saveas(gcf,strcat(".\three stim share figs\BD_Hist_ff_NoStim_",regions(i),".png"))
    
    histogram((hfs5_ff_bd),(binEdge))
%     histogram(log10(hfs5_ff_bd),100)
    ylabel("Counts")
    xlabel("BD (ms)")
    ax=gca;
    ax.XScale="log";
    ax.FontSize=18;
    title(strcat("Feedforward BD ",regions(i)))
    saveas(gcf,strcat(".\three stim share figs\BD_Hist_ff_5HFS_",regions(i),".png"))
    
    histogram((hfs40_ff_bd),(binEdge))
%     histogram(log10(hfs40_ff_bd),100)
    ylabel("Counts")
    xlabel("BD (ms)")
    ax=gca;
    ax.XScale="log";
    ax.FontSize=18;
    title(strcat("Feedforward BD ",regions(i)))
    saveas(gcf,strcat(".\three stim share figs\BD_Hist_ff_40HFS_",regions(i),".png"))
    
    % log bins
    
    hist_vec=log10([nostim_ff_bd;hfs5_ff_bd;hfs40_ff_bd]);
    stim_group=[repmat(stims(1),length(nostim_ff_bd),1);repmat(stims(2),length(hfs5_ff_bd),1);...
        repmat(stims(3),length(hfs40_ff_bd),1)];
    
    figure
    
    histogram(log10(nostim_ff_bd),log10(binEdge))
%     histogram(log10(nostim_ff_bd),100)
    ylabel("Counts")
    xlabel("Log(BD)")
    title(strcat("Feedforward BD ",regions(i)))
    saveas(gcf,strcat(".\three stim share figs\log_BD_Hist_ff_NoStim_",regions(i),".png"))
    
    histogram(log10(hfs5_ff_bd),log10(binEdge))
%     histogram(log10(hfs5_ff_bd),100)
    ylabel("Counts")
    xlabel("Log(BD)")
    title(strcat("Feedforward BD ",regions(i)))
    saveas(gcf,strcat(".\three stim share figs\log_BD_Hist_ff_5HFS_",regions(i),".png"))
    
    histogram(log10(hfs40_ff_bd),log10(binEdge))
%     histogram(log10(hfs40_ff_bd),100)
    ylabel("Counts")
    xlabel("Log(BD)")
    title(strcat("Feedforward BD ",regions(i)))
    saveas(gcf,strcat(".\three stim share figs\log_BD_Hist_ff_40HFS_",regions(i),".png"))
    
%     figure
    [~,~,stats]=anova1(hist_vec,stim_group);
    saveas(gcf,strcat(".\three stim share figs\BD_Box_ff_",regions(i),".png"))
    figure
    [ff_c,ff_means]=multcompare(stats);
    
    diff1=((10^ff_means(2,1)-10^ff_means(1,1))/10^ff_means(1,1))*100;
    diff2=((10^ff_means(3,1)-10^ff_means(1,1))/10^ff_means(1,1))*100;
    diff3=((10^ff_means(3,1)-10^ff_means(2,1))/10^ff_means(2,1))*100;
    BD_ff_diff{i}=[diff1;diff2;diff3];
    
    BD_ff_p{i}=ff_c;
    BD_ff_m{i}=ff_means;
    
    figure
    hold on
    
    b=bar(categorical(stims,stims),ff_means(:,1),'FaceColor','flat');
    b.CData=[[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250]];
    e=errorbar(ff_means(:,1),ff_means(:,2),'k');
    e.LineStyle='none';
    hold off
    
    ylim([1.5,3])
    ylabel("Log(BD)")
    xlabel("Stimulation")
    title("Feedforward BD "+regions(i))
    
    saveas(gcf,strcat(".\three stim share figs\BD_Bar_ff_",regions(i),".png"))
    
    %close all
end

figure
hold on

means=[]; err=[];
for i=1:length(BD_ff_m)
    means=[means;BD_ff_m{i}(:,1)'];
    err=[err;BD_ff_m{i}(:,2)'];
end

b=bar(categorical(regions,regions),means,'FaceColor','flat');
% b.CData=[[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250]];

[ngroups,nbars] = size(means);
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end

e=errorbar(x',means,err,'k','LineStyle','none');
hold off

ylim([1.5,2.5])
ylabel("Log(Burst Duration)")
xlabel("Subregion")
title("Feedforward Burst Duration ")

saveas(gcf,strcat(".\three stim share figs\BD_Bar_ff.png"))
%% Feedback

regions=["DG-EC","CA3-DG","CA1-CA3","EC-CA1"];
stims=["No Stim", "5HFS", "40HFS"];

BD_fb_p=[];

% binEdge = logspace(0,4,100);

for i=1:length(regions)
    nostim_fb_bd=cell2mat(nostim.spike_burst_dyn_table_stim.BurstDuration...
        (nostim.spike_burst_dyn_table_stim.regi==i & nostim.spike_burst_dyn_table_stim.if_ff==0));
    hfs5_fb_bd=cell2mat(hfs5.spike_burst_dyn_table_stim.BurstDuration...
        (hfs5.spike_burst_dyn_table_stim.regi==i & hfs5.spike_burst_dyn_table_stim.if_ff==0));
    hfs40_fb_bd=cell2mat(hfs40.spike_burst_dyn_table_stim.BurstDuration...
        (hfs40.spike_burst_dyn_table_stim.regi==i & hfs40.spike_burst_dyn_table_stim.if_ff==0));
    
    histogram((nostim_fb_bd),(binEdge))
%     histogram(log10(nostim_fb_bd),100)
    ylabel("Counts")
    xlabel("BD (ms)")
    ax=gca;
    ax.XScale="log";
    ax.FontSize=18;
    title(strcat("Feedback BD ",regions(i)))
    saveas(gcf,strcat(".\three stim share figs\BD_Hist_fb_NoStim_",regions(i),".png"))
    
    histogram((hfs5_fb_bd),(binEdge))
%     histogram(log10(hfs5_fb_bd),100)
    ylabel("Counts")
    xlabel("BD (ms)")
    ax=gca;
    ax.XScale="log";
    ax.FontSize=18;
    title(strcat("Feedback BD ",regions(i)))
    saveas(gcf,strcat(".\three stim share figs\BD_Hist_fb_5HFS_",regions(i),".png"))
    
    histogram((hfs40_fb_bd),(binEdge))
%     histogram(log10(hfs40_fb_bd),100)
    ylabel("Counts")
    xlabel("BD (ms)")
    ax=gca;
    ax.XScale="log";
    ax.FontSize=18;
    title(strcat("Feedback BD ",regions(i)))
    saveas(gcf,strcat(".\three stim share figs\BD_Hist_fb_40HFS_",regions(i),".png"))
    
    %log bin
    
    hist_vec=log10([nostim_fb_bd;hfs5_fb_bd;hfs40_fb_bd]);
    stim_group=[repmat(stims(1),length(nostim_fb_bd),1);repmat(stims(2),length(hfs5_fb_bd),1);...
        repmat(stims(3),length(hfs40_fb_bd),1)];
    
    histogram(log10(nostim_fb_bd),log10(binEdge))
%     histogram(log10(nostim_fb_bd),100)
    ylabel("Counts")
    xlabel("Log(BD)")
    title(strcat("Feedback BD ",regions(i)))
    saveas(gcf,strcat(".\three stim share figs\log_BD_Hist_fb_NoStim_",regions(i),".png"))
    
    histogram(log10(hfs5_fb_bd),log10(binEdge))
%     histogram(log10(hfs5_fb_bd),100)
    ylabel("Counts")
    xlabel("Log(BD)")
    title(strcat("Feedback BD ",regions(i)))
    saveas(gcf,strcat(".\three stim share figs\log_BD_Hist_fb_5HFS_",regions(i),".png"))
    
    histogram(log10(hfs40_fb_bd),log10(binEdge))
%     histogram(log10(hfs40_fb_bd),100)
    ylabel("Counts")
    xlabel("Log(BD)")
    title(strcat("Feedback BD ",regions(i)))
    saveas(gcf,strcat(".\three stim share figs\log_BD_Hist_fb_40HFS_",regions(i),".png"))
    
%     figure
    [~,~,stats]=anova1(hist_vec,stim_group);
    saveas(gcf,strcat(".\three stim share figs\BD_Box_fb_",regions(i),".png"))
    figure
    [ff_c,ff_means]=multcompare(stats);
    
    diff1=((10^ff_means(2,1)-10^ff_means(1,1))/10^ff_means(1,1))*100;
    diff2=((10^ff_means(3,1)-10^ff_means(1,1))/10^ff_means(1,1))*100;
    diff3=((10^ff_means(3,1)-10^ff_means(2,1))/10^ff_means(2,1))*100;
    BD_fb_diff{i}=[diff1;diff2;diff3];
    
    BD_fb_p{i}=ff_c;
    BD_fb_m{i}=ff_means;
    
    figure
    hold on
    
    b=bar(categorical(stims,stims),ff_means(:,1),'FaceColor','flat');
    b.CData=[[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250]];
    e=errorbar(ff_means(:,1),ff_means(:,2),'k');
    e.LineStyle='none';
    hold off
    
    ylim([1.5,3])
    ylabel("Log(BD)")
    xlabel("Stimulation")
    title("Feedback BD "+regions(i))
    
    saveas(gcf,strcat(".\three stim share figs\BD_Bar_fb_",regions(i),".png"))
    
    %close all
end

figure
hold on

means=[]; err=[];
for i=1:length(BD_fb_m)
    means=[means;BD_fb_m{i}(:,1)'];
    err=[err;BD_fb_m{i}(:,2)'];
end

b=bar(categorical(regions,regions),means,'FaceColor','flat');
% b.CData=[[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250]];

[ngroups,nbars] = size(means);
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end

e=errorbar(x',means,err,'k','LineStyle','none');
hold off

ylim([1.5,2.5])
ylabel("Log(Burst Duration)")
xlabel("Subregion")
title("Feedback Burst Duration ")

saveas(gcf,strcat(".\three stim share figs\BD_Bar_fb.png"))
%% BD Well
cd("D:\Brewer lab data\HFS\")
nostim=load('D:\Brewer lab data\HFS\No Stim\Wells_5SD_500maxSD\well_spike_dynamics_table_nostim.mat');
hfs5=load('D:\Brewer lab data\HFS\Theta Stim\Wells_5SD_500maxSD\well_spike_dynamics_table_theta.mat');
hfs40=load('D:\Brewer lab data\HFS\HFS Stim\Wells_5SD_500maxSD\well_spike_dynamics_table_hfs.mat');

regions=["EC","DG","CA3","CA1"];
stims=["No Stim", "5HFS", "40HFS"];

BD_well_p=[];

binEdge = logspace(0,4,200);

for i=1:length(regions)
    nostim_ff_bd=cell2mat(nostim.well_spike_dynamics_table.BurstDuration...
        (nostim.well_spike_dynamics_table.regi==i));
    hfs5_ff_bd=cell2mat(hfs5.well_spike_dynamics_table.BurstDuration...
        (hfs5.well_spike_dynamics_table.regi==i));
    hfs40_ff_bd=cell2mat(hfs40.well_spike_dynamics_table.BurstDuration...
        (hfs40.well_spike_dynamics_table.regi==i));
    
    histogram((nostim_ff_bd),(binEdge))
%     histogram(log10(nostim_ff_bd),100)
    ylabel("Counts")
    xlabel("BD (ms)")
    ax=gca;
    ax.XScale="log";
    ax.FontSize=18;
    title(strcat("Well BD ",regions(i)))
    saveas(gcf,strcat(".\three stim share figs\BD_Hist_well_NoStim_",regions(i),".png"))
    
    histogram((hfs5_ff_bd),(binEdge))
%     histogram(log10(hfs5_ff_bd),100)
    ylabel("Counts")
    xlabel("BD (ms)")
    ax=gca;
    ax.XScale="log";
    ax.FontSize=18;
    title(strcat("Well BD ",regions(i)))
    saveas(gcf,strcat(".\three stim share figs\BD_Hist_well_5HFS_",regions(i),".png"))
    
    histogram((hfs40_ff_bd),(binEdge))
%     histogram(log10(hfs40_ff_bd),100)
    ylabel("Counts")
    xlabel("BD (ms)")
    ax=gca;
    ax.XScale="log";
    ax.FontSize=18;
    title(strcat("Well BD ",regions(i)))
    saveas(gcf,strcat(".\three stim share figs\BD_Hist_well_40HFS_",regions(i),".png"))
    
    %log bin
    
    hist_vec=log10([nostim_ff_bd;hfs5_ff_bd;hfs40_ff_bd]);
    stim_group=[repmat(stims(1),length(nostim_ff_bd),1);repmat(stims(2),length(hfs5_ff_bd),1);...
        repmat(stims(3),length(hfs40_ff_bd),1)];
    
    histogram(log10(nostim_ff_bd),log10(binEdge))
%     histogram(log10(nostim_ff_bd),100)
    ylabel("Counts")
    xlabel("Log(BD)")
    title(strcat("Well BD ",regions(i)))
    saveas(gcf,strcat(".\three stim share figs\log_BD_Hist_well_NoStim_",regions(i),".png"))
    
    histogram(log10(hfs5_ff_bd),log10(binEdge))
%     histogram(log10(hfs5_ff_bd),100)
    ylabel("Counts")
    xlabel("Log(BD)")
    title(strcat("Well BD ",regions(i)))
    saveas(gcf,strcat(".\three stim share figs\log_BD_Hist_well_5HFS_",regions(i),".png"))
    
    histogram(log10(hfs40_ff_bd),log10(binEdge))
%     histogram(log10(hfs40_ff_bd),100)
    ylabel("Counts")
    xlabel("Log(BD)")
    title(strcat("Well BD ",regions(i)))
    saveas(gcf,strcat(".\three stim share figs\log_BD_Hist_well_40HFS_",regions(i),".png"))
    
%     figure
    [~,~,stats]=anova1(hist_vec,stim_group);
    saveas(gcf,strcat(".\three stim share figs\BD_Box_well_",regions(i),".png"))
    figure
    [ff_c,ff_means]=multcompare(stats);
    
    diff1=((10^ff_means(2,1)-10^ff_means(1,1))/10^ff_means(1,1))*100;
    diff2=((10^ff_means(3,1)-10^ff_means(1,1))/10^ff_means(1,1))*100;
    diff3=((10^ff_means(3,1)-10^ff_means(2,1))/10^ff_means(2,1))*100;
    BD_well_diff{i}=[diff1;diff2;diff3];
    
    BD_well_p{i}=ff_c;
    BD_well_m{i}=ff_means;
    
    figure
    hold on
    
    b=bar(categorical(stims,stims),ff_means(:,1),'FaceColor','flat');
    b.CData=[[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250]];
    e=errorbar(ff_means(:,1),ff_means(:,2),'k');
    e.LineStyle='none';
    hold off
    
    ylim([1.5,3])
    ylabel("Log(BD)")
    xlabel("Stimulation")
    title("Well BD "+regions(i))
    
    saveas(gcf,strcat(".\three stim share figs\BD_Bar_well_",regions(i),".png"))
    
    %close all
end

figure
hold on

means=[]; err=[];
for i=1:length(BD_well_m)
    means=[means;BD_well_m{i}(:,1)'];
    err=[err;BD_well_m{i}(:,2)'];
end

b=bar(categorical(regions,regions),means,'FaceColor','flat');
% b.CData=[[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250]];

[ngroups,nbars] = size(means);
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end

e=errorbar(x',means,err,'k','LineStyle','none');
hold off

ylim([1.5,2.5])
ylabel("Log(Burst Duration)")
xlabel("Subregion")
title("Well Burst Duration ")

saveas(gcf,strcat(".\three stim share figs\BD_Bar_well.png"))

