function [h, p, stats, n, stderr_spon, stderr_stim]=no_pad_wilcoxon_210701(table1, table2)

%if conduction time isnt present in both, drop time and spikes etc
for i=1:length(table1)
    [rwtb,~]=size(table1{i});
    for j=1:rwtb
        
        time_similar_1=zeros(1,length(table1{i}.("Conduction Time"){j}));
        for k=1:length(table2{i}.("Conduction Time"){j})
             time_similar_1=or(time_similar_1,table1{i}.("Conduction Time"){j}==table2{i}.("Conduction Time"){j}(k)+.05);
             time_similar_1=or(time_similar_1,table1{i}.("Conduction Time"){j}==table2{i}.("Conduction Time"){j}(k)-.05);
             time_similar_1=or(time_similar_1,table1{i}.("Conduction Time"){j}==table2{i}.("Conduction Time"){j}(k));
%             time_similar_1=or(time_similar_1,ismember(string(table2{i}.("Conduction Time"){j}(k)+.05),string(table1{i}.("Conduction Time"){j})));
%             time_similar_1=or(time_similar_1,ismember(string(table2{i}.("Conduction Time"){j}(k)-.05),string(table1{i}.("Conduction Time"){j})));
%             time_similar_1=or(time_similar_1,ismember(string(table2{i}.("Conduction Time"){j}(k)),string(table1{i}.("Conduction Time"){j})));
        end
        
        time_similar_2=zeros(1,length(table2{i}.("Conduction Time"){j}));
        for k=1:length(table1{i}.("Conduction Time"){j})
             time_similar_2=or(time_similar_2,table2{i}.("Conduction Time"){j}==table1{i}.("Conduction Time"){j}(k)+.05);
             time_similar_2=or(time_similar_2,table2{i}.("Conduction Time"){j}==table1{i}.("Conduction Time"){j}(k)-.05);
             time_similar_2=or(time_similar_2,table2{i}.("Conduction Time"){j}==table1{i}.("Conduction Time"){j}(k));
%             time_similar_2=or(time_similar_2,ismember(string(table1{i}.("Conduction Time"){j}(k)+.05),string(table2{i}.("Conduction Time"){j})));
%             time_similar_2=or(time_similar_2,ismember(string(table1{i}.("Conduction Time"){j}(k)-.05),string(table2{i}.("Conduction Time"){j})));
%             time_similar_2=or(time_similar_2,ismember(string(table1{i}.("Conduction Time"){j}(k)),string(table2{i}.("Conduction Time"){j})));
        end
        
        if ~isempty(time_similar_1) && any(time_similar_1)
            table1{i}.("Spike Area"){j}=table1{i}.("Spike Area"){j}(time_similar_1);
            table1{i}.("Conduction Time"){j}=table1{i}.("Conduction Time"){j}(time_similar_1);
            table1{i}.("Unit Pairs"){j}=table1{i}.("Unit Pairs"){j}(time_similar_1);
            table1{i}.Direction{j}=table1{i}.Direction{j}(time_similar_1);
        elseif ~any(time_similar_1)
            table1{i}.("Spike Area"){j}=[];
            table1{i}.("Conduction Time"){j}=[];
            table1{i}.("Unit Pairs"){j}=[];
            table1{i}.Direction{j}=[];
        end
        
        if ~isempty(time_similar_2) && any(time_similar_2)
            table2{i}.("Spike Area"){j}=table2{i}.("Spike Area"){j}(time_similar_2);
            table2{i}.("Conduction Time"){j}=table2{i}.("Conduction Time"){j}(time_similar_2);
            table2{i}.("Unit Pairs"){j}=table2{i}.("Unit Pairs"){j}(time_similar_2);
            table2{i}.Direction{j}=table2{i}.Direction{j}(time_similar_2);
        elseif ~any(time_similar_2)
            table2{i}.("Spike Area"){j}=[];
            table2{i}.("Conduction Time"){j}=[];
            table2{i}.("Unit Pairs"){j}=[];
            table2{i}.Direction{j}=[];
        end
        
    end
    
end

%% Compute ttest on summation of all similar datasets across all recordings

[total_h, total_p, total_stats,total_n, total_stderr_spon, total_stderr_stim]=total_wilcox_210630(table1,table2);
%%
[total_h_up, total_p_up, total_stats_up,total_n_up, total_stderr_up_spon, total_stderr_up_stim]=total_wilcox_upreg_210630(table1,table2);
%%
[total_h_down, total_p_down, total_stats_down,total_n_down, total_stderr_down_spon, total_stderr_down_stim]=total_wilcox_downreg_210630(table1,table2);

%%
h={total_h total_h_up total_h_down};
p={total_p total_p_up total_p_down};
stats={total_stats total_stats_up total_stats_down};
n={total_n total_n_up total_n_down};
stderr_spon={total_stderr_spon total_stderr_up_spon total_stderr_down_spon};
stderr_stim={total_stderr_stim total_stderr_up_stim total_stderr_down_stim};
end