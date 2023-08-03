function ttest_to_bar_210226(table1,table2,ttable)
% Assumes tables 1 and 2 are already organized
% IMPORTANT NOTE: UP AND DOWN REGULATION WERE SWITCHED DURING CREATION,
% NAMES IN SAVED PICTURES SWITCHED BUT VARIABLES HAVE NOT BEEN

if all(~ismember({'Electrode Pairs','Spike Area','Conduction Time','Unit Pairs','Direction'}, table1.Properties.VariableNames))...
        && all(~ismember({'Electrode Pairs','Spike Area','Conduction Time','Unit Pairs','Direction'}, table2.Properties.VariableNames))
    error('This table does not have the appropriate variables: Electrode Pairs,Spike Area,Conduction Time,Unit Pairs,Direction.')
end

sub_conx=unique(table1.Subregion);
ci_low_two_tail=[];
ci_low_right=[];
ci_low_left=[];

ci_hig_two_tail=[];
ci_hig_right=[];
ci_hig_left=[];


for i=1:length(sub_conx)
    if length([table1.("Spike Area"){table1.Subregion==sub_conx(i)}])~=length([table2.("Spike Area"){table2.Subregion==sub_conx(i)}])
        error('The tables conain a different number of spike area elements. Try zero padding the data.')
    end
    if length([table1.("Spike Area"){table1.Subregion==sub_conx(i)}])~=length([table2.("Spike Area"){table2.Subregion==sub_conx(i)}])
        error('The tables conain a different number of spike area elements. Try zero padding the data.')
    end
    
    range=table1.Subregion==sub_conx(i);
    
    t1_spikes(i)=sum([table1.("Spike Area"){range}]);
    t2_spikes(i)=sum([table2.("Spike Area"){range}]);
    
    if isempty(t1_spikes(i))
        t1_spikes(i)=0;
    end
    
    if isempty(t2_spikes(i))
        t2_spikes(i)=0;
    end
    
    %greater_spikes=[greater_spikes,max([t1_spikes(i),t2_spikes(i)])];
    
    ci_low_two_tail=[ci_low_two_tail,ttable.ci{1}(i,1)];
    ci_low_right=[ci_low_right,ttable.ci{2}(i,1)];
    ci_low_left=[ci_low_left,ttable.ci{3}(i,1)];
    
    ci_hig_two_tail=[ci_hig_two_tail,ttable.ci{1}(i,2)];
    ci_hig_right=[ci_hig_right,ttable.ci{2}(i,2)];
    ci_hig_left=[ci_hig_left,ttable.ci{3}(i,2)];
end

bar_axis=categorical(cellstr(sub_conx));

low_weight=0.25;

figure
hold on
bar(bar_axis,t1_spikes)
bar(bar_axis,t2_spikes,low_weight)
er=errorbar(bar_axis,t1_spikes,ci_low_two_tail,ci_hig_two_tail);
er.Color=[0 0 0];
er.LineStyle='none';
ylabel('Spikes')
xlabel('Region')
title({strcat(table1.("Recording Name"){1}(1:34), ' '), strcat(ttable.direction(1), ' T-Test Results Two Tail')})
hold off
legend('Spontaneous','Stimulated')

two_tail_bar=gcf;
saveas(two_tail_bar,fullfile('D:\Brewer lab data\ttest files',strcat(table1.("Recording Name"){1}(1:34), ' ', ttable.direction(1), ' T-Test Results Two Tail')),'png')

figure
hold on
bar(bar_axis,t1_spikes)
bar(bar_axis,t2_spikes,low_weight)
er=errorbar(bar_axis,t1_spikes,ci_low_right,ci_hig_right);
er.Color=[0 0 0];
er.LineStyle='none';
ylabel('Spikes')
xlabel('Region')
title({strcat(table1.("Recording Name"){1}(1:34), ' '), strcat(ttable.direction(1), ' T-Test Results Downregulation')})
hold off
legend('Spontaneous','Stimulated')

upreg_bar=gcf;
saveas(upreg_bar,fullfile('D:\Brewer lab data\ttest files',strcat(table1.("Recording Name"){1}(1:34), ' ', ttable.direction(1), ' T-Test Results Downregulation')),'png')

figure
hold on
bar(bar_axis,t1_spikes)
bar(bar_axis,t2_spikes,low_weight)
er=errorbar(bar_axis,t1_spikes,ci_low_left,ci_hig_left);
er.Color=[0 0 0];
er.LineStyle='none';
ylabel('Spikes')
xlabel('Region')
title({strcat(table1.("Recording Name"){1}(1:34), ' '), strcat(ttable.direction(1), ' T-Test Results Upregulation')})
hold off
legend('Spontaneous','Stimulated')

downreg_bar=gcf;
saveas(downreg_bar,fullfile('D:\Brewer lab data\ttest files',strcat(table1.("Recording Name"){1}(1:34), ' ', ttable.direction(1), ' T-Test Results Upregulation')),'png')
end