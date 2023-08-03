function [allregion_unit_matched]=peak_train_assembly_210914(table_ff, table_fb, time_idx)

for rec=1:length(table_ff)
    
    [rwtb,~]=size(table_ff{rec});
    
    times_table=table('Size',[rwtb,10],'VariableTypes',{'string','string','double','double','cell','cell','cell','cell','cell','cell'},...
        'VariableNames',{'Subregion','Electrode Pairs','regi','chani','up_ff','down_ff','up_fb','down_fb','ff_cdt','fb_cdt'});
    times_table.Subregion=table_ff{rec}.Subregion;
    times_table.("Electrode Pairs")=table_ff{rec}.("Electrode Pairs");
    
    subregions=unique(table_ff{rec}.Subregion,'stable');
    
    for rows=1:rwtb
        for i=1:length(subregions)
            if strcmpi(times_table.Subregion(rows),subregions(i))
                times_table.regi(rows)=i;
            end
            region_locs=find(strcmpi(times_table.Subregion,subregions(i)));
            for chid=1:sum(strcmpi(times_table.Subregion,subregions(i)))
                times_table.chani(region_locs(chid))=chid;
            end
        end
    end
    
    %feedforward
    for i=1:(rwtb)
        
        for j=1:length(table_ff{rec}.up_times{i})
            up=ismember(time_idx{rec},int64(table_ff{rec}.up_times{i}{j}.*1000));
            down=ismember(time_idx{rec},int64(table_ff{rec}.down_times{i}{j}.*1000));
            times_table.up_ff{i}{j}=up;
            times_table.down_ff{i}{j}=down;
            times_table.ff_cdt{i}{j}=table_ff{rec}.("Conduction Time"){i}(j);
        end
    end

    [rwtb,~]=size(table_fb{rec});
    
    %feedback
    for i=1:(rwtb)
            
        for j=1:length(table_fb{rec}.up_times{i})
            up=ismember(time_idx{rec},int64(table_fb{rec}.up_times{i}{j}.*1000));
            down=ismember(time_idx{rec},int64(table_fb{rec}.down_times{i}{j}.*1000));
            times_table.up_fb{i}{j}=up;
            times_table.down_fb{i}{j}=down;
            times_table.fb_cdt{i}{j}=table_fb{rec}.("Conduction Time"){i}(j);
        end
    end
    allregion_unit_matched{rec}=times_table;
end

end