function directional_pie_chart_210413(areas)
% Make sure input is a cell array where cells contain a singular table with
% fields: Subregion, Electrode Pairs, Spike Area, Conduction Time, Unit
% Pairs, Direction, Recording Name, and fid. This function from that
% creates pie charts of the direction of each channel.
subregion=unique(areas{1}.Subregion,'stable');
direction_count=table([0;0;0],[0;0;0],[0;0;0],[0;0;0],...
        'VariableNames',{'EC-DG','DG-CA3','CA3-CA1','CA1-EC'},...
        'RowNames',{'Forward','Backward','Both'});
for i=1:length(areas)
    for j=1:length(subregion)
        for k=1:sum(areas{i}.Subregion==subregion(j))
            subregion_locs=find(areas{i}.Subregion==subregion(j));
            directions=areas{i}.Direction{subregion_locs(k)};
            %Both directions
            if any(ismember(directions,0))...
                    && any(ismember(directions,1))
                direction_count.(subregion(j))(3)=direction_count.(subregion(j))(3)+1;
            %Feedforward
            elseif any(ismember(directions,1))
                direction_count.(subregion(j))(1)=direction_count.(subregion(j))(1)+1;
            %Feedback
            elseif any(ismember(directions,0))
                direction_count.(subregion(j))(2)=direction_count.(subregion(j))(2)+1;
            end
        end
    end
end

row_names=direction_count.Properties.RowNames;
for i=1:length(subregion)
    figure
    counts=direction_count.(subregion(i));
    p=pie(counts);
    pText = findobj(p,'Type','text');
    percentValues = get(pText,'String'); 
%     ffperc=direction_count.(subregion(i))(1)/sum(counts);
%     fbperc=direction_count.(subregion(i))(2)/sum(counts);
%     bothperc=direction_count.(subregion(i))(3)/sum(counts);
%    percentages={ffperc,fbperc,bothperc};
    colons={': ';': ';': '};
    nums={['('+string(counts(1))+') '];['('+string(counts(2))+') '];['('+string(counts(3))+') ']};
    
    for j=1:length(nums)
        nums{j}=convertStringsToChars(nums{j});
    end
    
    combinedtxt=strcat(row_names,colons,nums,percentValues);
    pText(1).String = combinedtxt(1);
    pText(2).String = combinedtxt(2);
    pText(3).String = combinedtxt(3);
    newColors=[...
        1, 0, 0;
        0, 0, 1;
        .5, .5, .5];
    ax=gca();
    ax.Colormap=newColors;
    title(subregion(i))
    for j=1:length(nums)
        t=p(2*j);
        t.FontSize=20;
    end
    
end
end