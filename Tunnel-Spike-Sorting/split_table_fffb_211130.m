function [spon_areas_table_ff,spon_areas_table_fb]=split_table_fffb_211130(spon_areas_table)

spon_areas_table_ff={};
spon_areas_table_fb={};
for i=1:length(spon_areas_table)
    if strmatch('up_times',spon_areas_table{i}.Properties.VariableNames)~=0
        spon_areas_table_ff{i}=table('Size',size(spon_areas_table{i}),'VariableTypes',{'string','string','cell','cell','cell','cell','cell','double','cell','cell'},'VariableNames',{'Subregion','Electrode Pairs','Spike Area','Conduction Time',...
        'Unit Pairs', 'Direction', 'Recording Name','fid','up_times','down_times'});
        spon_areas_table_fb{i}=table('Size',size(spon_areas_table{i}),'VariableTypes',{'string','string','cell','cell','cell','cell','cell','double','cell','cell'},'VariableNames',{'Subregion','Electrode Pairs','Spike Area','Conduction Time',...
        'Unit Pairs', 'Direction', 'Recording Name','fid','up_times','down_times'});
    else
        spon_areas_table_ff{i}=table('Size',size(spon_areas_table{i}),'VariableTypes',{'string','string','cell','cell','cell','cell','cell','double'},'VariableNames',{'Subregion','Electrode Pairs','Spike Area','Conduction Time',...
        'Unit Pairs', 'Direction', 'Recording Name','fid'});
        spon_areas_table_fb{i}=table('Size',size(spon_areas_table{i}),'VariableTypes',{'string','string','cell','cell','cell','cell','cell','double'},'VariableNames',{'Subregion','Electrode Pairs','Spike Area','Conduction Time',...
        'Unit Pairs', 'Direction', 'Recording Name','fid'});
    end

    for j=1:length(spon_areas_table{i}.("Direction"))
        if isempty(spon_areas_table{i}.("Direction"){j})
            spon_areas_table_ff{i}.("Subregion"){j}=spon_areas_table{i}.("Subregion"){j};
            spon_areas_table_fb{i}.("Subregion"){j}=spon_areas_table{i}.("Subregion"){j};
            spon_areas_table_ff{i}.("Electrode Pairs"){j}=spon_areas_table{i}.("Electrode Pairs"){j};
            spon_areas_table_fb{i}.("Electrode Pairs"){j}=spon_areas_table{i}.("Electrode Pairs"){j};
            spon_areas_table_ff{i}.("Spike Area"){j}=spon_areas_table{i}.("Spike Area"){j};
            spon_areas_table_fb{i}.("Spike Area"){j}=spon_areas_table{i}.("Spike Area"){j};
            spon_areas_table_ff{i}.("Conduction Time"){j}=spon_areas_table{i}.("Conduction Time"){j};
            spon_areas_table_fb{i}.("Conduction Time"){j}=spon_areas_table{i}.("Conduction Time"){j};
            spon_areas_table_ff{i}.("Unit Pairs"){j}=spon_areas_table{i}.("Unit Pairs"){j};
            spon_areas_table_fb{i}.("Unit Pairs"){j}=spon_areas_table{i}.("Unit Pairs"){j};
            spon_areas_table_ff{i}.("Direction"){j}=spon_areas_table{i}.("Direction"){j};
            spon_areas_table_fb{i}.("Direction"){j}=spon_areas_table{i}.("Direction"){j};
            spon_areas_table_ff{i}.("Recording Name"){j}=spon_areas_table{i}.("Recording Name"){j};
            spon_areas_table_fb{i}.("Recording Name"){j}=spon_areas_table{i}.("Recording Name"){j};
            spon_areas_table_ff{i}.("fid")(j)=spon_areas_table{i}.("fid")(j);
            spon_areas_table_fb{i}.("fid")(j)=spon_areas_table{i}.("fid")(j);
            if strmatch('up_times',spon_areas_table{i}.Properties.VariableNames)~=0
                spon_areas_table_ff{i}.("up_times"){j}=spon_areas_table{i}.("up_times"){j};
                spon_areas_table_fb{i}.("up_times"){j}=spon_areas_table{i}.("up_times"){j};
                spon_areas_table_ff{i}.("down_times"){j}=spon_areas_table{i}.("down_times"){j};
                spon_areas_table_fb{i}.("down_times"){j}=spon_areas_table{i}.("down_times"){j};
            end
%             spon_areas_table_ff{i}(j,["Electrode Pairs",'Spike Area','Conduction Time','Unit Pairs','Direction','Recording Name'])=...
%                 spon_areas_table{i}(j,["Electrode Pairs",'Spike Area','Conduction Time','Unit Pairs','Direction','Recording Name']);
%             spon_areas_table_fb{i}(j,["Electrode Pairs",'Spike Area','Conduction Time','Unit Pairs','Direction','Recording Name'])=...
%                 spon_areas_table{i}(j,["Electrode Pairs",'Spike Area','Conduction Time','Unit Pairs','Direction','Recording Name']);
            
            continue
        end
        spon_areas_table_ff{i}.("Subregion"){j}=spon_areas_table{i}.("Subregion"){j};
        spon_areas_table_fb{i}.("Subregion"){j}=spon_areas_table{i}.("Subregion"){j};
        spon_areas_table_ff{i}.("Electrode Pairs"){j}=spon_areas_table{i}.("Electrode Pairs"){j};
        spon_areas_table_fb{i}.("Electrode Pairs"){j}=spon_areas_table{i}.("Electrode Pairs"){j};
        spon_areas_table_ff{i}.("fid")(j)=spon_areas_table{i}.("fid")(j);
        spon_areas_table_fb{i}.("fid")(j)=spon_areas_table{i}.("fid")(j);
        for k=1:length(spon_areas_table{i}.("Direction"){j})
            if spon_areas_table{i}.("Direction"){j}(k)==1
                spon_areas_table_ff{i}.("Spike Area"){j}=[spon_areas_table_ff{i}.("Spike Area"){j},spon_areas_table{i}.("Spike Area"){j}(k)];
                spon_areas_table_ff{i}.("Conduction Time"){j}=[spon_areas_table_ff{i}.("Conduction Time"){j}, spon_areas_table{i}.("Conduction Time"){j}(k)];
                spon_areas_table_ff{i}.("Unit Pairs"){j}=[spon_areas_table_ff{i}.("Unit Pairs"){j},spon_areas_table{i}.("Unit Pairs"){j}(k)];
                spon_areas_table_ff{i}.("Direction"){j}=[spon_areas_table_ff{i}.("Direction"){j}, spon_areas_table{i}.("Direction"){j}(k)];
                if strmatch('up_times',spon_areas_table{i}.Properties.VariableNames)~=0
                    spon_areas_table_ff{i}.("up_times"){j}=[spon_areas_table_ff{i}.("up_times"){j},spon_areas_table{i}.("up_times"){j}(k)];
                    spon_areas_table_ff{i}.("down_times"){j}=[spon_areas_table_ff{i}.("down_times"){j},spon_areas_table{i}.("down_times"){j}(k)];
                end
                
            elseif spon_areas_table{i}.("Direction"){j}(k)==0
                spon_areas_table_fb{i}.("Spike Area"){j}=[spon_areas_table_fb{i}.("Spike Area"){j},spon_areas_table{i}.("Spike Area"){j}(k)];
                spon_areas_table_fb{i}.("Conduction Time"){j}=[spon_areas_table_fb{i}.("Conduction Time"){j}, spon_areas_table{i}.("Conduction Time"){j}(k)];
                spon_areas_table_fb{i}.("Unit Pairs"){j}=[spon_areas_table_fb{i}.("Unit Pairs"){j}, spon_areas_table{i}.("Unit Pairs"){j}(k)];
                spon_areas_table_fb{i}.("Direction"){j}=[spon_areas_table_fb{i}.("Direction"){j}, spon_areas_table{i}.("Direction"){j}(k)];
                if strmatch('up_times',spon_areas_table{i}.Properties.VariableNames)~=0
                    spon_areas_table_fb{i}.("up_times"){j}=[spon_areas_table_fb{i}.("up_times"){j},spon_areas_table{i}.("up_times"){j}(k)];
                    spon_areas_table_fb{i}.("down_times"){j}=[spon_areas_table_fb{i}.("down_times"){j},spon_areas_table{i}.("down_times"){j}(k)];
                end
                
            end

        end
        spon_areas_table_ff{i}.("Recording Name"){j}=spon_areas_table{i}.("Recording Name"){j};
        spon_areas_table_fb{i}.("Recording Name"){j}=spon_areas_table{i}.("Recording Name"){j};
    end
end

end