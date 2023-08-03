function [pre_stim_vals,post_stim_vals]=pre_post_time_comp_nopad_210511(table1,table2)
% Tables 1 and 2 should be pre stim and post stim in that order

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

pre_stim_vals=[];
post_stim_vals=[];
for i=1:length(table1)
    pre_stim_vals=[pre_stim_vals,[table1{i}.("Spike Area"){1:end}]];
    post_stim_vals=[post_stim_vals,[table2{i}.("Spike Area"){1:end}]];
end

%300 is the number of seconds in the recording we are using, so as a note:
%HARD CODED
% pre_stim_vals=pre_stim_vals./300;
% post_stim_vals=post_stim_vals./300;

%Changed to no division, calculates in per seconds earlier
pre_stim_vals=pre_stim_vals;
post_stim_vals=post_stim_vals;