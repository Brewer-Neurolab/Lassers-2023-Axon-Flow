function [anova_subregion,anova_ff,anova_fb]=network_stats_211110(ff_table,fb_table)

%get number of subregions
subregions=unique(ff_table{1}.Subregion,'stable');

for i=1:length(subregions)
    current_subregion=subregions(i);
    num_axons=0;
    axons_per_tunnel=[];
    conduction_times=[];
    for j=1:length(ff_table)
        table=ff_table{j}(ff_table{j}.Subregion==current_subregion,:);
        num_axons=num_axons+length([table.("Conduction Time"){:}]);
        for k=1:length([table.("Conduction Time")])
            axons_per_tunnel=[axons_per_tunnel,length([table.("Conduction Time"){k}])];
        end
        axons_per_tunnel=axons_per_tunnel(axons_per_tunnel~=0);
        conduction_times=[conduction_times,[table.("Conduction Time"){:}]];
    end
    
    ff_axons(i)=num_axons;
    ff_axons_per_tunnel{i}=axons_per_tunnel;
    ff_conduction_times{i}=conduction_times;
    
    num_axons=0;
    conduction_times=[];
    for j=1:length(fb_table)
        table=fb_table{j}(fb_table{j}.Subregion==current_subregion,:);
        num_axons=num_axons+length([table.("Conduction Time"){:}]);
        for k=1:length([table.("Conduction Time")])
            axons_per_tunnel=[axons_per_tunnel,length([table.("Conduction Time"){k}])];
        end
        axons_per_tunnel=axons_per_tunnel(axons_per_tunnel~=0);
        conduction_times=[conduction_times,[table.("Conduction Time"){:}]];
    end
    
    fb_axons(i)=num_axons;
    fb_axons_per_tunnel{i}=axons_per_tunnel;
    fb_conduction_times{i}=conduction_times;
    
end

%pie charts
figure
for k=1:length(subregions)
    L={strcat('n=',string(ff_axons(k))),strcat('n=',string(fb_axons(k)))};
    ax(k)=subplot(1,length(subregions),k);
    subplot(ax(k))
    H=pie([ff_axons(k),fb_axons(k)]);
    title(subregions(k))
    colormap([1,0,0;0,0,1])
    T = H(strcmpi(get(H,'Type'),'text'));
    P = cell2mat(get(T,'Position'));
    set(T,{'Position'},num2cell(P*0.5,2))
    P(1,1)=P(1,1)*1.5;
    text(P(:,1),P(:,2),L(:),'FontSize',15)
    axplot=gca;
    axplot.FontSize=20;
end
set(gcf,'Units','normalized','OuterPosition',[0 0 1 1])

%count bar graphs
figure
title('Axons Per Tunnel')
xlabel('Number of axons per tunnel')
ylabel('Counts')
loops=1;
for k=1:2:2*length(subregions)
    ax1(k)=subplot(1,2*length(subregions),k);
    ax2(k)=subplot(1,2*length(subregions),k+1);
    subplot(ax1(k))
    [GC,GR]=groupcounts(ff_axons_per_tunnel{loops}');
    x=(1:max(GR))';
    y=[];
    y(ismember(x,GR))=GC;
    y(~ismember(x,GR))=0;
    bar(x,y,'r')
    title(subregions(loops)+' FF')
    %xlabel('Number of axons per tunnel')
    %ylabel('Counts')
    
    ax=gca;
    ax.FontSize=10;
    
    subplot(ax2(k))
    [GC,GR]=groupcounts(fb_axons_per_tunnel{loops}');
    bar(GR,GC,'b')
    title(subregions(loops)+' FB')
    %xlabel('Number of axons per tunnel')
    %ylabel('Counts')
    
    ax=gca;
    ax.FontSize=10;
    
    loops=loops+1;
end
set(gcf,'Units','normalized','OuterPosition',[0 0 1 1])
[~,h1]=suplabel('Number of Axons Per Tunnel','t');
[~,h2]=suplabel('Number of Axons');
[~,h3]=suplabel('Count','y');

set(h1,'Fontsize',16)
set(h2,'Fontsize',16)
set(h3,'Fontsize',16)

%link axes for bar graph 
linkaxes([ax1,ax2],'y')

%convert conduction times to conduction velocity
for i=1:length(subregions)
    ff_conduction_times{i}=(1./(ff_conduction_times{i}/1000)).*(200/1e6);
    fb_conduction_times{i}=(1./(ff_conduction_times{i}/1000)).*(200/1e6);
end

% %compute 2way anova for direction
for k=1:length(subregions)
    %[p(k),tbl(k),stats{k}]=anova2([ff_conduction_times{k}',fb_conduction_times{k}']);
    y=[ff_conduction_times{k},fb_conduction_times{k}]';
    g=[repmat("Feedforward",1,length(ff_conduction_times{k})),...
        repmat("Feedback",1,length(fb_conduction_times{k}))]';
    [pval,anovatable,anovastats]=anovan(y,{g});
    p{k}=pval;
    tbl{k}=anovatable;
    stats{k}=anovastats;
end
anova_subregion={p,tbl,stats};

%compute 2way anova for ff and fb between all subregions
% y=[ff_conduction_times{:}]';
% g=[];
% for i=1:length(subregions)
%     g=repmat(subregions(i),1,length(ff_conduction_times{i}));
%     varnames{i}=subregions(i);
% end
% [ff_p,ff_tbl,ff_stats]=anovan([ff_conduction_times{:}]',{g'});
% 
% y=[fb_conduction_times{:}]';
% g=[];
% for i=1:length(subregions)
%     g=[g,repmat(subregions(i),1,length(fb_conduction_times{i}))];
% end
% [fb_p,fb_tbl,fb_stats]=anovan([fb_conduction_times{:}]',{g'});
% 
% anova_ff={ff_p,ff_tbl,ff_stats};
% anova_fb={fb_p,fb_tbl,fb_stats};

paired_permute=nchoosek(subregions,2);
paired_permute_idx=nchoosek([1:length(subregions)],2);
%feedforward
for i=1:length(paired_permute(:,1))
    y=[ff_conduction_times{paired_permute_idx(i,1)},ff_conduction_times{paired_permute_idx(i,2)}]';
    g=[repmat(subregions(paired_permute_idx(i,1)),1,length(ff_conduction_times{paired_permute_idx(i,1)}))...
        repmat(subregions(paired_permute_idx(i,2)),1,length(ff_conduction_times{paired_permute_idx(i,2)}))]';
    comparing=strcat(subregions(paired_permute_idx(i,1)),'|',subregions(paired_permute_idx(i,2)));
    [p,tbl,stats]=anovan(y,{g},'model',2,'varnames',comparing);
    anova_ff(i).name=comparing;
    anova_ff(i).p=p;
    anova_ff(i).tbl=tbl;
    anova_ff(i).stats=stats;
end
%feedback
for i=1:length(paired_permute(:,1))
    y=[fb_conduction_times{paired_permute_idx(i,1)},fb_conduction_times{paired_permute_idx(i,2)}]';
    g=[repmat(subregions(paired_permute_idx(i,1)),1,length(fb_conduction_times{paired_permute_idx(i,1)}))...
        repmat(subregions(paired_permute_idx(i,2)),1,length(fb_conduction_times{paired_permute_idx(i,2)}))]';
    comparing=strcat(subregions(paired_permute_idx(i,1)),'|',subregions(paired_permute_idx(i,2)));
    [p,tbl,stats]=anovan(y,{g},'model',2,'varnames',comparing);
    anova_fb(i).name=comparing;
    anova_fb(i).p=p;
    anova_fb(i).tbl=tbl;
    anova_fb(i).stats=stats;
end

%Compute two way anova
% g1=[];
% g2=[];
% g3=[];
% g4=[];
% for i=1:length(subregions)
%     g1=[g1,repmat(subregions(i),1,length(ff_conduction_times{i}))];
%     g2=[g2,repmat(subregions(i),1,length(fb_conduction_times{i}))];
%     g3=[g3,repmat("Feedforward",1,length(ff_conduction_times{i}))];
%     g4=[g4,repmat("Feedback",1,length(fb_conduction_times{i}))];
% end
% 
% group_subregion=[g1,g2]';
% group_direction=[g3,g4]';
% y=[[ff_conduction_times{:}],fb_conduction_times{:}]';
% [p,tbl,stats]=anovan(y,{group_subregion, group_direction},'model',2,'varnames',{'Subregioin','Direction'});

%bar([mean([ff_conduction_times{:}]);mean([fb_conduction_times{:}])],'r','b')
vals1=[];
vals2=[];
vals3=[];
vals4=[];
for i=1:length(subregions)
    % 200um, conduction times in ms
    vals1=[vals1;(mean(ff_conduction_times{i}))];
    vals2=[vals2;(mean(fb_conduction_times{i}))];
    vals3=[vals3;std(ff_conduction_times{i})/sqrt(length(ff_conduction_times{i}))];
    vals4=[vals4;std(fb_conduction_times{i})/sqrt(length(fb_conduction_times{i}))];
end

vals=[vals1,vals2];
error=[vals3,vals4];
figure
hold on
b=bar(vals);
b(1).FaceColor='r';
b(2).FaceColor='b';
ngroups = size(vals, 1);
nbars = size(vals, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, vals(:,i), error(:,i), '.k');
end
hold off
labels=[];
for i=1:length(subregions)
    split_name=strsplit(subregions(i),{'-'});
    name=strcat(subregions(i),'|',strcat(split_name(2),'-',split_name(1)));
    labels{i}=name;
end
xticks([1:length(labels)])
xticklabels(labels)
set(gcf,'Units','normalized','OuterPosition',[0 0 1 1])
title('Conduction Velocity')
ylabel('Conduction Velocity (m/s)')
ax=gca;
ax.FontSize=20;
end