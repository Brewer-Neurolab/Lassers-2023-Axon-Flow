function [zeropad_table1, zeropad_table2]=zero_pad_table_210804(table1,table2)

% Zero pads two tables with identical fields
% Assumes Combine has been run first so the two vectors are of equal length
% and the spikes are ordered acending 

if all(~ismember({'Subregion','Electrode Pairs','Spike Area','Conduction Time','Unit Pairs','Direction'}, table1.Properties.VariableNames))...
        && all(~ismember({'Subregion','Electrode Pairs','Spike Area','Conduction Time','Unit Pairs','Direction'}, table2.Properties.VariableNames))
    error('This table does not have the appropriate variables: Electrode Pairs,Spike Area,Conduction Time,Unit Pairs,Direction.')
end

[rwtb1, ~]=size(table1);
[rwtb2, ~]=size(table2);

if ~isequal(rwtb1,rwtb2)
    error("Rows across tables are not equal.")
end

for i=1:rwtb1
    % Check if number of elements in area, time, and direction are the same
%     if ~isequal(numel(table1.("Spike Area"){i}),numel(table1.("Conduction Time"){i}),numel(table1.Direction{i}))...
%             || ~isequal(numel(table2.("Spike Area"){i}),numel(table2.("Conduction Time"){i}),numel(table2.Direction{i}))
%         error('Table does not have same number of elements across area, time and direction')
%     end
    
    if isempty(table1.("Conduction Time"){i}) && ~isempty(table2.("Conduction Time"){i})
        table1.("Conduction Time"){i}=table2.("Conduction Time"){i};
        table1.("Spike Area"){i}=zeros(length(table2.("Conduction Time"){i}),1);
        table1.Direction{i}=all([table1.Direction{1:end}]);
        table1.("Unit Pairs"){i}=strcat(table1.("Electrode Pairs"){i},'-u1');
        if strmatch('up_times',table1.Properties.VariableNames)~=0 && strmatch('up_times',table2.Properties.VariableNames)~=0
            table1.up_times{i}=cell(length(table2.("Conduction Time"){i}),1);
            table1.down_times{i}=cell(length(table2.("Conduction Time"){i}),1);
        end

    elseif ~isempty(table1.("Conduction Time"){i}) && isempty(table2.("Conduction Time"){i})
        table2.("Conduction Time"){i}=table1.("Conduction Time"){i};
        table2.("Spike Area"){i}=zeros(length(table1.("Conduction Time"){i}),1);
        table2.Direction{i}=all([table2.Direction{1:end}]);
        table2.("Unit Pairs"){i}=strcat(table2.("Electrode Pairs"){i},'-u1');
        if strmatch('up_times',table1.Properties.VariableNames)~=0 && strmatch('up_times',table2.Properties.VariableNames)~=0
            table2.up_times{i}=cell(length(table2.("Conduction Time"){i}),1);
            table2.down_times{i}=cell(length(table1.("Conduction Time"){i}),1);
        end
    end
       

    temp_times_1=table1.("Conduction Time"){i};
    temp_times_2=table2.("Conduction Time"){i};
    
    %items to add to the first times
    %convert all items to strings to avoid being difficult
    to_cat_1=[];
    for j=1:length(temp_times_2)
        upbound=temp_times_2(j)+.08;
        downbound=temp_times_2(j)-.08;
        comp=or(abs(upbound-temp_times_1)<0.00001,abs(downbound-temp_times_1)<0.00001);
        comp=or(comp,abs(temp_times_2(j)-temp_times_1)<0.00001);
        if ~any(comp)
            to_cat_1=[to_cat_1,temp_times_2(j)];
        end
    end
    
    to_cat_2=[];
    for j=1:length(temp_times_1)
        upbound=temp_times_1(j)+.08;
        downbound=temp_times_1(j)-.08;
        comp=or(abs(upbound-temp_times_2)<0.00001,abs(downbound-temp_times_2)<0.00001);
        comp=or(comp,abs(temp_times_1(j)-temp_times_2)<0.00001);
        if ~any(comp)
            to_cat_2=[to_cat_2,temp_times_1(j)];
        end
    end
    
    temp_times_1=sort([temp_times_1,to_cat_1]);
    temp_times_2=sort([temp_times_2,to_cat_2]);
    
%     original_times_locs_1={};
%     for j=1:length(temp_times_1)
%         original_times_locs_1{j}=(temp_times_1(j)==table1.("Conduction Time"){i});
%     end
%     
%     original_times_locs_2={};
%     for j=1:length(temp_times_2)
%         original_times_locs_2{j}=(temp_times_2(j)==table2.("Conduction Time"){i});
%     end
    

    original_times_locs_1=zeros(1,length(temp_times_1));
    original_times_locs_2=zeros(1,length(temp_times_2));
    
    for j=1:length(table1.("Conduction Time"){i})
        original_times_locs_1=or(original_times_locs_1,(temp_times_1==table1.("Conduction Time"){i}(j)));
    end
    
    for j=1:length(table2.("Conduction Time"){i})
        original_times_locs_2=or(original_times_locs_2,(temp_times_2==table2.("Conduction Time"){i}(j)));
    end
    
    newSpikes1=zeros(1,length(temp_times_1));
    newSpikes2=zeros(1,length(temp_times_2));
    
    if ~isempty(original_times_locs_1)
        newSpikes1(original_times_locs_1)=table1.("Spike Area"){i};
    end
    
    if ~isempty(original_times_locs_2)
        newSpikes2(original_times_locs_2)=table2.("Spike Area"){i};
    end
    
    if strmatch('up_times',table1.Properties.VariableNames)~=0 && strmatch('up_times',table2.Properties.VariableNames)~=0
        newUp1=cell(length(temp_times_1),1);
        newUp2=cell(length(temp_times_2),1);

        if ~isempty(original_times_locs_1)
            newUp1(original_times_locs_1)=table1.up_times{i};
        end

        if ~isempty(original_times_locs_2)
            newUp2(original_times_locs_2)=table2.up_times{i};
        end

        newDown1=cell(length(temp_times_1),1);
        newDown2=cell(length(temp_times_2),1);

        if ~isempty(original_times_locs_1)
            newDown1(original_times_locs_1)=table1.down_times{i};
        end

        if ~isempty(original_times_locs_2)
            newDown2(original_times_locs_2)=table2.down_times{i};
        end
        table1.up_times{i}=newUp1;
        table2.up_times{i}=newUp2;
        table1.down_times{i}=newDown1;
        table2.down_times{i}=newDown2;
        
        if ~isempty(table1.up_times{i})
        
    end
    
    table1.("Conduction Time"){i}=temp_times_1;
    table2.("Conduction Time"){i}=temp_times_2;
    table1.("Spike Area"){i}=newSpikes1;
    table2.("Spike Area"){i}=newSpikes2;
    
end

for i=1:rwtb1
    %Re-evaluates unit pairs
    table1.("Unit Pairs"){i}=[];
    table2.("Unit Pairs"){i}=[];

    for j=1:length(table1.("Spike Area"){i})
        table1.("Unit Pairs"){i}=[table1.("Unit Pairs"){i};strcat(table1.("Electrode Pairs"){i},"-u",string(j))];
    end

    for j=1:length(table2.("Spike Area"){i})
        table2.("Unit Pairs"){i}=[table2.("Unit Pairs"){i};strcat(table2.("Electrode Pairs"){i},"-u",string(j))];
    end
end

direction=all([[table1.Direction{1:end}],[table2.Direction{1:end}]]);
for i=1:rwtb1
   %Re-evaluates directions 
    if direction==1
        table1.Direction{i}=ones(1,length(table1.("Conduction Time"){i}));
        table2.Direction{i}=ones(1,length(table2.("Conduction Time"){i}));
    elseif direction==0
        table1.Direction{i}=zeros(1,length(table1.("Conduction Time"){i}));
        table2.Direction{i}=zeros(1,length(table2.("Conduction Time"){i}));
    end
    
end

zeropad_table1=table1;
zeropad_table2=table2;

end