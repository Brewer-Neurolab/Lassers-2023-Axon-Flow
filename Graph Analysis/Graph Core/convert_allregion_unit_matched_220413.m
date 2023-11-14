function new_matched=convert_allregion_unit_matched_220413(allregion_unit_matched)

new_matched=[];
for i=1:length(allregion_unit_matched)
    recording_data=[];
    Subregion=unique(allregion_unit_matched{i}.Subregion,'stable');
    for j=1:length(Subregion)
        tunnels=cellstr(allregion_unit_matched{i}.("Electrode Pairs")(allregion_unit_matched{i}.Subregion==Subregion(j)));
        original_chan=allregion_unit_matched{i}.("Electrode Pairs")(allregion_unit_matched{i}.Subregion==Subregion(j));
        clusters=[];
        direction=[];
        ff=allregion_unit_matched{i}.up_ff(allregion_unit_matched{i}.Subregion==Subregion(j));
        fb=allregion_unit_matched{i}.up_fb(allregion_unit_matched{i}.Subregion==Subregion(j));
        for k=1:length(tunnels)
            clusters{k}=[[1:length(ff{k});1:length(ff{k})]';[1:length(fb{k});1:length(fb{k})]'];
            direction{k}=[ones(length(ff{k}),1);zeros(length(fb{k}),1)];
        end
        %alphabetized order
        tunnels_str=string((tunnels));
        first_tunnel_electrode=split(tunnels_str,{'-'});
%         [~,idx]=sort(first_tunnel_electrode(:,1));
        first_tunnel_electrode=cellstr(first_tunnel_electrode(:,1));
        R = cell2mat(regexp(first_tunnel_electrode ,'(?<Name>\D+)(?<Nums>\d+)','names'));
        tmp = sortrows([{R.Name}' num2cell(cellfun(@(x)str2double(x),{R.Nums}'))]);
        SortedText = strcat(tmp(:,1) ,cellfun(@(x) num2str(x), tmp(:,2),'unif',0));
        
        idx=[];
        for elec=1:length(tunnels)
            idx(elec)=find(contains(string(tunnels),string(SortedText{elec})));
        end
        
        clusters=clusters(idx);
        tunnels=tunnels(idx);
        direction=direction(idx);
        original_chan=original_chan(idx);
        
        recording_data{j}=[clusters',tunnels,direction',original_chan'];
        
    end
    new_matched{i}=recording_data;
end

end