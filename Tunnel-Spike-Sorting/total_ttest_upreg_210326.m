function [h, p, ci, stats, n, stderr1, stderr2]=total_ttest_upreg_210326(table1, table2)

% Assumes tables 1 and 2 are already organized
% input the metadata from the spontaneous and stimulated datasets into this function

% if all(~ismember({'Electrode Pairs','Spike Area','Conduction Time','Unit Pairs','Direction'}, table1.Properties.VariableNames))...
%         && all(~ismember({'Electrode Pairs','Spike Area','Conduction Time','Unit Pairs','Direction'}, table2.Properties.VariableNames))
%     error('This table does not have the appropriate variables: Electrode Pairs,Spike Area,Conduction Time,Unit Pairs,Direction.')
% end

% Modify t1 and t2 to only have values in stim that go up
for i=1:length(table1)
    [rwtb,~]=size(table1{i});
    for j=1:rwtb
        stim_greater=table1{i}.("Spike Area"){j}<=table2{i}.("Spike Area"){j};
        table1{i}.("Spike Area"){j}=table1{i}.("Spike Area"){j}(stim_greater)./300;
        table1{i}.("Conduction Time"){j}=table1{i}.("Conduction Time"){j}(stim_greater);
        table1{i}.("Unit Pairs"){j}=table1{i}.("Unit Pairs"){j}(stim_greater);
        table1{i}.Direction{j}=table1{i}.Direction{j}(stim_greater);
        
        table2{i}.("Spike Area"){j}=table2{i}.("Spike Area"){j}(stim_greater)./300;
        table2{i}.("Conduction Time"){j}=table2{i}.("Conduction Time"){j}(stim_greater);
        table2{i}.("Unit Pairs"){j}=table2{i}.("Unit Pairs"){j}(stim_greater);
        table2{i}.Direction{j}=table2{i}.Direction{j}(stim_greater);
        
    end
    
end


sub_conx=unique(table1{1}.Subregion,'stable');

%Number of elements in each column of bar graph
n=zeros(1,length(sub_conx));
for i=1:length(sub_conx)
    
    %assumes number of n's are equal for both tables
    for j=1:length(table1)
        n(i)=sum([n(i),length([table1{j}.("Spike Area"){table1{j}.Subregion==sub_conx(i)}])]);
    end
    
end

n={sub_conx;n};

%zero values seed
t1_spike_all=zeros(1,length(sub_conx));
t2_spike_all=zeros(1,length(sub_conx));

for i=1:length(sub_conx)
%     if length([all_table1.("Spike Area"){all_table1.Subregion==sub_conx(i)}])~=length([table2.("Spike Area"){table2.Subregion==sub_conx(i)}])
%         error('The tables conain a different number of spike area elements. Try zero padding the data.')
%     end
    
    t1_spike_vector=[];
    for j=1:length(table1)
        t1_spike_all(i)=sum([t1_spike_all(i),[table1{j}.("Spike Area"){table1{j}.Subregion==sub_conx(i)}]]);
        t1_spike_vector=[t1_spike_vector,[table1{j}.("Spike Area"){table1{j}.Subregion==sub_conx(i)}]];
        
    end
    if length(t1_spike_vector)<3
        t1_spike_all(i)=0;
        t1_spike_vector=0;
    end
    spike_mean_1(i)=mean(t1_spike_vector);
    
    t2_spike_vector=[];
    for j=1:length(table2)
        t2_spike_all(i)=sum([t2_spike_all(i),[table2{j}.("Spike Area"){table2{j}.Subregion==sub_conx(i)}]]);
        t2_spike_vector=[t2_spike_vector,[table2{j}.("Spike Area"){table2{j}.Subregion==sub_conx(i)}]];
    end
    if length(t2_spike_vector)<3
        t2_spike_all(i)=0;
        t2_spike_vector=0;
    end
    spike_mean_2(i)=mean(t2_spike_vector);
    
    if ~isempty(t1_spike_vector)&&~isempty(t2_spike_vector)
        [h(i) p(i) ci(i,:) stats{i}]=ttest(t1_spike_vector,t2_spike_vector);
        tempstderr1=std(t1_spike_vector)/sqrt(length(t1_spike_vector));
        stderr1(i,:)=[-tempstderr1,tempstderr1];
        tempstderr2=std(t2_spike_vector)/sqrt(length(t2_spike_vector));
        stderr2(i,:)=[-tempstderr2,tempstderr2];
    end
end

%bar_axis=categorical(cellstr(sub_conx));
%bar_axis=reordercats(bar_axis,cellstr(sub_conx));
bar_axis=[1,2,3,4];
means=[spike_mean_1;spike_mean_2];

%low_weight=0.25;
fontsize=15;

figure
hold on
%bar(bar_axis,spike_mean_1)
%bar(bar_axis,spike_mean_2,low_weight)
b=bar(bar_axis,means);

% Hard code distances
er1=errorbar(bar_axis-.15,spike_mean_1,stderr1(:,1),stderr1(:,2));
er1.Color=[0 0 0];
er1.LineStyle='none';

er2=errorbar(bar_axis+.15,spike_mean_2,stderr2(:,1),stderr2(:,2));
er2.Color=[0 0 0];
er2.LineStyle='none';

ylabel('Hz')
xlabel('Region')
if all([table1{5}.Direction{1:end}])
    title("Upregulation: Feed Forward",'FontSize',fontsize)
elseif ~all([table1{5}.Direction{1:end}])
    title("Upregulation: Feed Back",'FontSize',fontsize)
else
    title("Upregulation: CHECK DIRECTION",'FontSize',fontsize)
end

ax=gca;
ax.XAxis.FontSize=fontsize;
ax.YAxis.FontSize=fontsize;
xtips1=b(1).XEndPoints;
ytips1=b(1).YEndPoints;
xtips2=b(2).XEndPoints;
ytips2=b(2).YEndPoints;
labels=string(repelem("n=", length(sub_conx)))+n{2};
text(xtips1,ytips1,labels,'HorizontalAlignment','center','VerticalAlignment','bottom')
text(xtips2,ytips2,labels,'HorizontalAlignment','center','VerticalAlignment','bottom')
space_ax=repelem(" ", length(sub_conx));
xnames=[sub_conx';space_ax];
xnames=xnames(:)';
xnames(end)=[];
%ax.XTickLabel=xnames;
%set(gca, 'XTickLabel',{xnames});
xticklabels(xnames)

hold off

end