%% Loads and Plots Medians and SE for all IBSR and BD distributions
% Sam Lassers 
% Plots median IBSR and BD from loaded variables. Also creates distribution
% based off of the ln transformed data centroids.

%% Load Data
clear
clc

cd("D:\Brewer lab data\HFS\")

%Well
nostim_well_bd=load("nostim_well_BD_5SD_500maxSD.mat");
nostim_well_IBSR=load("nostim_well_IBSR_5SD_500maxSD.mat");
theta_well_bd=load("theta_well_BD_5SD_500maxSD.mat");
theta_well_IBSR=load("theta_well_IBSR_5SD_500maxSD.mat");
hfs_well_bd=load("HFS_well_BD_5SD_500maxSD.mat");
hfs_well_IBSR=load("HFS_well_IBSR_5SD_500maxSD.mat");

%Tunnel
nostim_ff_bd=load("nostim_ff_BD.mat");
nostim_ff_IBSR=load("nostim_ff_IBSR.mat");
nostim_fb_bd=load("nostim_fb_BD.mat");
nostim_fb_IBSR=load("nostim_fb_IBSR.mat");
theta_ff_bd=load("Theta_ff_BD.mat");
theta_ff_IBSR=load("Theta_ff_IBSR.mat");
theta_fb_bd=load("Theta_fb_BD.mat");
theta_fb_IBSR=load("Theta_fb_IBSR.mat");
hfs_ff_bd=load("HFS_ff_bd.mat");
hfs_ff_IBSR=load("HFS_ff_IBSR.mat");
hfs_fb_bd=load("HFS_fb_bd.mat");
hfs_fb_IBSR=load("HFS_fb_IBSR.mat");

ff_labels={'EC-DG','DG-CA3','CA3-CA1','CA1-EC'};
ff_cats=categorical(ff_labels);
ff_cats=reordercats(ff_cats,string(ff_labels));

fb_labels={'DG-EC','CA3-DG','CA1-CA3','EC-CA1'};
fb_cats=categorical(fb_labels);
fb_cats=reordercats(fb_cats,string(fb_labels));

well_labels={'EC','DG','CA3','CA1'};
well_cats=categorical(well_labels);
well_cats=reordercats(well_cats,string(well_labels));

ln_sd=@(sd,mean) sqrt(log((sd^2/(mean))+1));
%% Feedforward Burst Duration
figure
medians=[];
err=[];
for i=1:4
    transformed_median=[log(nostim_ff_bd.moment_table.median(i)),log(theta_ff_bd.moment_table.median(i)),...
        log(hfs_ff_bd.moment_table.median(i))];
    medians=[medians;transformed_median];
    transformed_sd=[(nostim_ff_bd.moment_table.lnsd(i))...
        ,(theta_ff_bd.moment_table.lnsd(i)),...
        (hfs_ff_bd.moment_table.lnsd(i))];
    
    transformed_sd(1)=transformed_sd(1)/sqrt(nostim_ff_bd.moment_table.n_fit(i));
    transformed_sd(2)=transformed_sd(2)/sqrt(theta_ff_bd.moment_table.n_fit(i));
    transformed_sd(3)=transformed_sd(3)/sqrt(hfs_ff_bd.moment_table.n_fit(i));
    
    err=[err;transformed_sd];
end

b=bar(medians,'grouped');
hold on
[ngroups,nbars]=size(medians);
x=nan(nbars,ngroups);
for i=1:nbars
    x(i,:)=b(i).XEndPoints;
end
errorbar(x',medians,err,'k','linestyle','none')

ylim([3,6])
% title("Feedforward Burst Durtation")
ylabel("ln(Burst Duration (ms) )")
xticklabels(ff_cats)

ax=gca;
ax.FontSize=16;

hold off

saveas(gcf,".\three stim share figs\median_ff_bd.png")
%% Feedback Burst Duration
figure
medians=[];
err=[];
for i=1:4
    transformed_median=[log(nostim_fb_bd.moment_table.median(i)),log(theta_fb_bd.moment_table.median(i)),...
        log(hfs_fb_bd.moment_table.median(i))];
    medians=[medians;transformed_median];
    transformed_sd=[(nostim_fb_bd.moment_table.lnsd(i))...
        ,(theta_fb_bd.moment_table.lnsd(i)),...
        (hfs_fb_bd.moment_table.lnsd(i))];
    
    transformed_sd(1)=transformed_sd(1)/sqrt(nostim_fb_bd.moment_table.n_fit(i));
    transformed_sd(2)=transformed_sd(2)/sqrt(theta_fb_bd.moment_table.n_fit(i));
    transformed_sd(3)=transformed_sd(3)/sqrt(hfs_fb_bd.moment_table.n_fit(i));
    
    err=[err;transformed_sd];
end

b=bar(medians,'grouped');
hold on
[ngroups,nbars]=size(medians);
x=nan(nbars,ngroups);
for i=1:nbars
    x(i,:)=b(i).XEndPoints;
end
errorbar(x',medians,err,'k','linestyle','none')

ylim([3,6])
% title("Feedback Burst Durtation")
ylabel("ln(Burst Duration (ms) )")
xticklabels(fb_cats)

ax=gca;
ax.FontSize=16;

hold off
saveas(gcf,".\three stim share figs\median_fb_bd.png")
%% Well Burst Duration
figure
medians=[];
err=[];
for i=1:4
    transformed_median=[log(nostim_well_bd.moment_table.median(i)),log(theta_well_bd.moment_table.median(i)),...
        log(hfs_well_bd.moment_table.median(i))];
    medians=[medians;transformed_median];
    transformed_sd=[(nostim_well_bd.moment_table.lnsd(i))...
        ,(theta_well_bd.moment_table.lnsd(i)),...
        (hfs_well_bd.moment_table.lnsd(i))];
    
    transformed_sd(1)=transformed_sd(1)/sqrt(nostim_well_bd.moment_table.n_fit(i));
    transformed_sd(2)=transformed_sd(2)/sqrt(theta_well_bd.moment_table.n_fit(i));
    transformed_sd(3)=transformed_sd(3)/sqrt(hfs_well_bd.moment_table.n_fit(i));
    
    err=[err;transformed_sd];
end

b=bar(medians,'grouped');
hold on
[ngroups,nbars]=size(medians);
x=nan(nbars,ngroups);
for i=1:nbars
    x(i,:)=b(i).XEndPoints;
end
errorbar(x',medians,err,'k','linestyle','none')

ylim([3,6])
% title("Well Burst Durtation")
ylabel("ln(Burst Duration (ms) )")
xticklabels(well_cats)

ax=gca;
ax.FontSize=16;

hold off
saveas(gcf,".\three stim share figs\median_well_bd.png")
%% Feedforward IBSR
figure
medians=[];
err=[];
for i=1:4
    transformed_median=[log(nostim_ff_IBSR.moment_table.median(i)),log(theta_ff_IBSR.moment_table.median(i)),...
        log(hfs_ff_IBSR.moment_table.median(i))];
    medians=[medians;transformed_median];
    transformed_sd=[(nostim_ff_IBSR.moment_table.lnsd(i))...
        ,(theta_ff_IBSR.moment_table.lnsd(i)),...
        (hfs_ff_IBSR.moment_table.lnsd(i))];
    
    transformed_sd(1)=transformed_sd(1)/sqrt(nostim_ff_IBSR.moment_table.n_fit(i));
    transformed_sd(2)=transformed_sd(2)/sqrt(theta_ff_IBSR.moment_table.n_fit(i));
    transformed_sd(3)=transformed_sd(3)/sqrt(hfs_ff_IBSR.moment_table.n_fit(i));
    
    err=[err;transformed_sd];
end

b=bar(medians,'grouped');
hold on
[ngroups,nbars]=size(medians);
x=nan(nbars,ngroups);
for i=1:nbars
    x(i,:)=b(i).XEndPoints;
end
errorbar(x',medians,err,'k','linestyle','none')

ylim([3,6])
% title("Feedforward IBSR")
ylabel("ln(IBSR (Hz) )")
xticklabels(ff_cats)

ax=gca;
ax.FontSize=16;

hold off
saveas(gcf,".\three stim share figs\median_ff_ibsr.png")
%% Feedback IBSR
figure
medians=[];
err=[];
for i=1:4
    transformed_median=[log(nostim_fb_IBSR.moment_table.median(i)),log(theta_fb_IBSR.moment_table.median(i)),...
        log(hfs_fb_IBSR.moment_table.median(i))];
    medians=[medians;transformed_median];
    transformed_sd=[(nostim_fb_IBSR.moment_table.lnsd(i))...
        ,(theta_fb_IBSR.moment_table.lnsd(i)),...
        (hfs_fb_IBSR.moment_table.lnsd(i))];
    
    transformed_sd(1)=transformed_sd(1)/sqrt(nostim_fb_IBSR.moment_table.n_fit(i));
    transformed_sd(2)=transformed_sd(2)/sqrt(theta_fb_IBSR.moment_table.n_fit(i));
    transformed_sd(3)=transformed_sd(3)/sqrt(hfs_fb_IBSR.moment_table.n_fit(i));
    
    err=[err;transformed_sd];
end

b=bar(medians,'grouped');
hold on
[ngroups,nbars]=size(medians);
x=nan(nbars,ngroups);
for i=1:nbars
    x(i,:)=b(i).XEndPoints;
end
errorbar(x',medians,err,'k','linestyle','none')

ylim([3,6])
% title("Feedback IBSR")
ylabel("ln(IBSR (Hz) )")
xticklabels(fb_cats)

ax=gca;
ax.FontSize=16;

hold off
saveas(gcf,".\three stim share figs\median_fb_ibsr.png")
%% Well ISBR
figure
medians=[];
err=[];
for i=1:4
    transformed_median=[log(nostim_well_IBSR.moment_table.median(i)),log(theta_well_IBSR.moment_table.median(i)),...
        log(hfs_well_IBSR.moment_table.median(i))];
    medians=[medians;transformed_median];
    transformed_sd=[(nostim_well_IBSR.moment_table.lnsd(i))...
        ,(theta_well_IBSR.moment_table.lnsd(i)),...
        (hfs_well_IBSR.moment_table.lnsd(i))];
    
    transformed_sd(1)=transformed_sd(1)/sqrt(nostim_well_IBSR.moment_table.n_fit(i));
    transformed_sd(2)=transformed_sd(2)/sqrt(theta_well_IBSR.moment_table.n_fit(i));
    transformed_sd(3)=transformed_sd(3)/sqrt(hfs_well_IBSR.moment_table.n_fit(i));
    
    err=[err;transformed_sd];
end

b=bar(medians,'grouped');
hold on
[ngroups,nbars]=size(medians);
x=nan(nbars,ngroups);
for i=1:nbars
    x(i,:)=b(i).XEndPoints;
end
errorbar(x',medians,err,'k','linestyle','none')

ylim([3,6])
% title("Well IBSR")
ylabel("ln(IBSR (Hz) )")
xticklabels(well_cats)

ax=gca;
ax.FontSize=16;

hold off
saveas(gcf,".\three stim share figs\median_well_ibsr.png")
%% BD Well Create distribution

for i=1:4
nostim_dist=normrnd(nostim_well_bd.moment_table.lnmu(i),nostim_well_bd.moment_table.lnsd(i),[nostim_well_bd.moment_table.n_fit(i),1]);
hfs5_dist=normrnd(theta_well_bd.moment_table.lnmu(i),theta_well_bd.moment_table.lnsd(i),[theta_well_bd.moment_table.n_fit(i),1]);
hfs40_dist=normrnd(hfs_well_bd.moment_table.lnmu(i),hfs_well_bd.moment_table.lnsd(i),[hfs_well_bd.moment_table.n_fit(i),1]);

y=[nostim_dist;hfs5_dist;hfs40_dist];
g=[repmat("nostim",length(nostim_dist),1);repmat("5HFS",length(hfs5_dist),1);repmat("40HFS",length(hfs40_dist),1)];

[~,~,stats]=anova1(y,g);
[c,m]=multcompare(stats);

bd_well_c{i}=c;
bd_well_m{i}=m;
end
%% BD FF dist
for i=1:4
nostim_dist=normrnd(nostim_ff_bd.moment_table.lnmu(i),nostim_ff_bd.moment_table.lnsd(i),[nostim_ff_bd.moment_table.n_fit(i),1]);
hfs5_dist=normrnd(theta_ff_bd.moment_table.lnmu(i),theta_ff_bd.moment_table.lnsd(i),[theta_ff_bd.moment_table.n_fit(i),1]);
hfs40_dist=normrnd(hfs_ff_bd.moment_table.lnmu(i),hfs_ff_bd.moment_table.lnsd(i),[hfs_ff_bd.moment_table.n_fit(i),1]);

y=[nostim_dist;hfs5_dist;hfs40_dist];
g=[repmat("nostim",length(nostim_dist),1);repmat("5HFS",length(hfs5_dist),1);repmat("40HFS",length(hfs40_dist),1)];

[~,~,stats]=anova1(y,g);
[c,m]=multcompare(stats);

bd_ff_c{i}=c;
bd_ff_m{i}=m;
end
%% BD FB dist
for i=1:4
nostim_dist=normrnd(nostim_fb_bd.moment_table.lnmu(i),nostim_fb_bd.moment_table.lnsd(i),[nostim_fb_bd.moment_table.n_fit(i),1]);
hfs5_dist=normrnd(theta_fb_bd.moment_table.lnmu(i),theta_fb_bd.moment_table.lnsd(i),[theta_fb_bd.moment_table.n_fit(i),1]);
hfs40_dist=normrnd(hfs_fb_bd.moment_table.lnmu(i),hfs_fb_bd.moment_table.lnsd(i),[hfs_fb_bd.moment_table.n_fit(i),1]);

y=[nostim_dist;hfs5_dist;hfs40_dist];
g=[repmat("nostim",length(nostim_dist),1);repmat("5HFS",length(hfs5_dist),1);repmat("40HFS",length(hfs40_dist),1)];

[~,~,stats]=anova1(y,g);
[c,m]=multcompare(stats);

bd_fb_c{i}=c;
bd_fb_m{i}=m;
end
%% IBSR well Dist
for i=1:4
nostim_dist=normrnd(nostim_well_IBSR.moment_table.lnmu(i),nostim_well_IBSR.moment_table.lnsd(i),[nostim_well_IBSR.moment_table.n_fit(i),1]);
hfs5_dist=normrnd(theta_well_IBSR.moment_table.lnmu(i),theta_well_IBSR.moment_table.lnsd(i),[theta_well_IBSR.moment_table.n_fit(i),1]);
hfs40_dist=normrnd(hfs_well_IBSR.moment_table.lnmu(i),hfs_well_IBSR.moment_table.lnsd(i),[hfs_well_IBSR.moment_table.n_fit(i),1]);

y=[nostim_dist;hfs5_dist;hfs40_dist];
g=[repmat("nostim",length(nostim_dist),1);repmat("5HFS",length(hfs5_dist),1);repmat("40HFS",length(hfs40_dist),1)];

[~,~,stats]=anova1(y,g);
[c,m]=multcompare(stats);

ibsr_well_c{i}=c;
ibsr_well_m{i}=m;
end
%% IBSR FF dist
for i=1:4
nostim_dist=normrnd(nostim_ff_IBSR.moment_table.lnmu(i),nostim_ff_IBSR.moment_table.lnsd(i),[nostim_ff_IBSR.moment_table.n_fit(i),1]);
hfs5_dist=normrnd(theta_ff_IBSR.moment_table.lnmu(i),theta_ff_IBSR.moment_table.lnsd(i),[theta_ff_IBSR.moment_table.n_fit(i),1]);
hfs40_dist=normrnd(hfs_ff_IBSR.moment_table.lnmu(i),hfs_ff_IBSR.moment_table.lnsd(i),[hfs_ff_IBSR.moment_table.n_fit(i),1]);

y=[nostim_dist;hfs5_dist;hfs40_dist];
g=[repmat("nostim",length(nostim_dist),1);repmat("5HFS",length(hfs5_dist),1);repmat("40HFS",length(hfs40_dist),1)];

[~,~,stats]=anova1(y,g);
[c,m]=multcompare(stats);

ibsr_ff_c{i}=c;
ibsr_ff_m{i}=m;
end
%% IBSR FB dist
for i=1:4
nostim_dist=normrnd(nostim_fb_IBSR.moment_table.lnmu(i),nostim_fb_IBSR.moment_table.lnsd(i),[nostim_fb_IBSR.moment_table.n_fit(i),1]);
hfs5_dist=normrnd(theta_fb_IBSR.moment_table.lnmu(i),theta_fb_IBSR.moment_table.lnsd(i),[theta_fb_IBSR.moment_table.n_fit(i),1]);
hfs40_dist=normrnd(hfs_fb_IBSR.moment_table.lnmu(i),hfs_fb_IBSR.moment_table.lnsd(i),[hfs_fb_IBSR.moment_table.n_fit(i),1]);

y=[nostim_dist;hfs5_dist;hfs40_dist];
g=[repmat("nostim",length(nostim_dist),1);repmat("5HFS",length(hfs5_dist),1);repmat("40HFS",length(hfs40_dist),1)];

[~,~,stats]=anova1(y,g);
[c,m]=multcompare(stats);

ibsr_fb_c{i}=c;
ibsr_fb_m{i}=m;
end

%% BD Well Diff

for i=1:4
    diff_well_bd{i}=[100*(diff([nostim_well_bd.moment_table.median(i),theta_well_bd.moment_table.median(i)])/nostim_well_bd.moment_table.median(i));...
        100*(diff([nostim_well_bd.moment_table.median(i),hfs_well_bd.moment_table.median(i)])/nostim_well_bd.moment_table.median(i));...
        100*(diff([theta_well_bd.moment_table.median(i),hfs_well_bd.moment_table.median(i)])/theta_well_bd.moment_table.median(i))];
end

%% BD FF Diff
for i=1:4
    diff_ff_bd{i}=[100*(diff([nostim_ff_bd.moment_table.median(i),theta_ff_bd.moment_table.median(i)])/nostim_ff_bd.moment_table.median(i));...
        100*(diff([nostim_ff_bd.moment_table.median(i),hfs_ff_bd.moment_table.median(i)])/nostim_ff_bd.moment_table.median(i));...
        100*(diff([theta_ff_bd.moment_table.median(i),hfs_ff_bd.moment_table.median(i)])/theta_ff_bd.moment_table.median(i))];
end
%% BD FB Diff
for i=1:4
    diff_fb_bd{i}=[100*(diff([nostim_fb_bd.moment_table.median(i),theta_fb_bd.moment_table.median(i)])/nostim_fb_bd.moment_table.median(i));...
        100*(diff([nostim_fb_bd.moment_table.median(i),hfs_fb_bd.moment_table.median(i)])/nostim_fb_bd.moment_table.median(i));...
        100*(diff([theta_fb_bd.moment_table.median(i),hfs_fb_bd.moment_table.median(i)])/theta_fb_bd.moment_table.median(i))];
end
%% IBSR Well Diff
for i=1:4
    diff_well_IBSR{i}=[100*(diff([nostim_well_IBSR.moment_table.median(i),theta_well_IBSR.moment_table.median(i)])/nostim_well_IBSR.moment_table.median(i));...
        100*(diff([nostim_well_IBSR.moment_table.median(i),hfs_well_IBSR.moment_table.median(i)])/nostim_well_IBSR.moment_table.median(i));...
        100*(diff([theta_well_IBSR.moment_table.median(i),hfs_well_IBSR.moment_table.median(i)])/theta_well_IBSR.moment_table.median(i))];
end
%% IBSR FF Diff
for i=1:4
    diff_ff_IBSR{i}=[100*(diff([nostim_ff_IBSR.moment_table.median(i),theta_ff_IBSR.moment_table.median(i)])/nostim_ff_IBSR.moment_table.median(i));...
        100*(diff([nostim_ff_IBSR.moment_table.median(i),hfs_ff_IBSR.moment_table.median(i)])/nostim_ff_IBSR.moment_table.median(i));...
        100*(diff([theta_ff_IBSR.moment_table.median(i),hfs_ff_IBSR.moment_table.median(i)])/theta_ff_IBSR.moment_table.median(i))];
end
%% IBSR Fb Diff
for i=1:4
    diff_fb_IBSR{i}=[100*(diff([nostim_fb_IBSR.moment_table.median(i),theta_fb_IBSR.moment_table.median(i)])/nostim_fb_IBSR.moment_table.median(i));...
        100*(diff([nostim_fb_IBSR.moment_table.median(i),hfs_fb_IBSR.moment_table.median(i)])/nostim_fb_IBSR.moment_table.median(i));...
        100*(diff([theta_fb_IBSR.moment_table.median(i),hfs_fb_IBSR.moment_table.median(i)])/theta_fb_IBSR.moment_table.median(i))];
end