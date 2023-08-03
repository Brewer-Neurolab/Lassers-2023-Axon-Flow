function [spon_areas]=reorganize_NMI_outputs_210905(channel_pairs, spon_areas,if_cw,parent_dir)

spon_areas_organized=channel_pairs;
for i=1:length(spon_areas)
    
    if if_cw(i)
        channel_pairs=load(strcat(parent_dir,"\","matching_table_cw.mat"));
        channel_pairs=table2cell(channel_pairs.matching_table);
        spon_areas_organized=(channel_pairs);
    elseif ~if_cw(i)
        channel_pairs=load(strcat(parent_dir,"\","matching_table_ccw.mat"));
        channel_pairs=table2cell(channel_pairs.matching_table);
        spon_areas_organized=(channel_pairs);
    else
        error("Please lable files correctly as cw or ccw.")
    end
    
    for j=1:length(channel_pairs)
        units=[];
        conduction_times=[];
        units_names=[];
        direction=[];
        file_name=[];
        subregion=[];
        if isfield(spon_areas{i},'up_times')
            up_times=[];
        end
        if isfield(spon_areas{i},'down_times')
            down_times=[];
        end
        for k=1:length(spon_areas{i}.channel_pair)
            if contains(string(spon_areas{i}.channel_pair(k)),string(channel_pairs(j,2))) && spon_areas{i}.area(k) ~= -1
                units=[units; spon_areas{i}.area(k)];
                units_names=[units_names; spon_areas{i}.channel_pair(k)];
                %subregion=[subregion; spon_areas{i}.sub_region(k)];
                direction=[direction; spon_areas{i}.direction(k)];
                conduction_times=[conduction_times;spon_areas{i}.conduction_time(k)];
                if isfield(spon_areas{i},'up_times')
                    up_times=[up_times;spon_areas{i}.up_times(k)];
                end
                if isfield(spon_areas{i},'down_times')
                    down_times=[down_times;spon_areas{i}.down_times(k)];
                end
                
            end
            file_name=[spon_areas{i}.file_name];
            rec_fid=spon_areas{i}.fid;
        end
        spon_areas_organized{j,3}=units;
        %spon_areas_organized{j,3}=subregion;
        spon_areas_organized{j,4}=conduction_times;
        spon_areas_organized{j,5}=units_names;
        spon_areas_organized{j,6}=direction;
        spon_areas_organized{j,7}=file_name;
        spon_areas_organized{j,8}=rec_fid;
        if isfield(spon_areas{i},'up_times')
            spon_areas_organized{j,9}=up_times;
            has_up_spon(i)=1;
        else
            has_up_spon(i)=0;
        end
        if isfield(spon_areas{i},'down_times')
            spon_areas_organized{j,10}=down_times;
            has_down_spon(i)=1;
        else
            has_down_spon(i)=0;
        end
    end
    % Saves over variable
    spon_areas{i}=spon_areas_organized;
end

end