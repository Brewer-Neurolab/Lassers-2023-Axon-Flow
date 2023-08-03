function [allregion_unit_matched]=directory_peak_train_assembly_220103(well_folder,time_idx,direction)

cd(well_folder)
current_folder=dir;
folders=current_folder([current_folder.isdir]);
folders=folders(3:end); %removes relative pathing
recordings=string({folders.name});

for rec=1:length(recordings)
    
    cd(strcat(well_folder,'\',recordings(rec)))
    
    if direction(rec)
        load(strcat(well_folder,'\','matching_table_wells_CW.mat'))
    elseif ~direction(rec)
        load(strcat(well_folder,'\','matching_table_wells_CCW.mat'))
    else
        error('Direction not specified correctly. -SBL')
    end
    
    rwtb=size(matching_table,1);
    
    times_table=table('Size',[rwtb,4],'VariableTypes',{'string','string','double','cell'},...
        'VariableNames',{'Subregion','Electrode','regi','Spike Train'});
    times_table.Subregion=matching_table.subregion; %need to load matching table
    times_table.("Electrode")=matching_table.electrode;
    
    subregions=unique(matching_table.subregion,'stable');
    
    %ensure matching table cw and ccw have the same subregion order
    for rows=1:rwtb
        for i=1:length(subregions)
            if strcmpi(times_table.Subregion(rows),subregions(i))
                times_table.regi(rows)=i;
            end
        end
    end
    
    %feedforward
    for i=1:(rwtb)
        load(strcat(well_folder,'\',recordings(rec),'\',matching_table.electrode(i),'_spikes.mat'))
        times_table.("Spike Train"){i}=ismember(time_idx{rec},int64(index.*1000));
    end

    allregion_unit_matched{rec}=times_table;
end

end
