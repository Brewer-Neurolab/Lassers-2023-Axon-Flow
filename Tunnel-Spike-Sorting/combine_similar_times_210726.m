function [combined_table]=combine_similar_times_210726(table)

%7/26/21 changed from 0.05 change in conduction time to 0.10

if all(~ismember({'Subregion','Electrode Pairs','Spike Area','Conduction Time','Unit Pairs','Direction'}, table.Properties.VariableNames))
    error('This table does not have the appropriate variables: Electrode Pairs,Spike Area,Conduction Time,Unit Pairs,Direction.')
end

[rwtb, ~]=size(table);

for i=1:rwtb
    
    range_to_sum={};
    if isempty(table.("Conduction Time"){i})||length(table.("Conduction Time"){i})==1
        continue
    end
    
    subset_sum=[];
    range_ticker=1;
    for j=length(table.("Conduction Time"){i}):-1:2
        if isempty(subset_sum)
            subset_sum=j;
        end
        %changed to tolerance comparison because MATLAB schenanigans
        %7/19/21
        if (abs((table.("Conduction Time"){i}(j)-0.1)-table.("Conduction Time"){i}(j-1))<0.00001)...
                ||(abs(table.("Conduction Time"){i}(j)-table.("Conduction Time"){i}(j-1))<0.00001)
            if j==2
                range_to_sum{range_ticker}=[subset_sum,1];
            else
                subset_sum=[subset_sum,j-1];
            end
        else
            if isempty(subset_sum)
                subset_sum=j;
            end
            if j==2
                range_to_sum{range_ticker}=subset_sum;
                range_ticker=range_ticker+1;
                range_to_sum{range_ticker}=1;
            else
                range_to_sum{range_ticker}=subset_sum;
                subset_sum=[];
                range_ticker=range_ticker+1;
            end
        end
        

    end
    
    temp_area=[];
    temp_time=[];
    temp_dir=[];
    temp_up=[];
    temp_down=[];
    
    for j=length(range_to_sum):-1:1
        temp_area=[temp_area, sum(table.("Spike Area"){i}(range_to_sum{j}))];
        [~,max_spike_idx]=max(table.("Spike Area"){i}(range_to_sum{j}));
        spikes_in_range=table.("Conduction Time"){i}(range_to_sum{j});
        temp_time=[temp_time, spikes_in_range(max_spike_idx)];
        temp_dir=[temp_dir,table.Direction{i}(1)];
        if strmatch('up_times',table.Properties.VariableNames)~=0
            temp_up=[temp_up,{[table.up_times{i}{range_to_sum{j}}]}];
            temp_down=[temp_down,{[table.down_times{i}{range_to_sum{j}}]}];
        end
    end
    
    table.("Spike Area"){i}=temp_area;
    table.("Conduction Time"){i}=temp_time;
    table.Direction{i}=temp_dir;
    table.("Unit Pairs"){i}=table.("Unit Pairs"){i}(1:length(table.("Conduction Time"){i}));
    if strmatch('up_times',table.Properties.VariableNames)~=0
        table.up_times{i}=temp_up;
        table.down_times{i}=temp_down;
    end
    
end

combined_table=table;

end
