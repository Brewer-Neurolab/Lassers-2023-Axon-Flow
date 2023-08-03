function [h_test,p_test,ci_test,stats_test, percent_diff]=unpaired_ttest_spikes_211130(table1,table2,table1_duration,table2_duration)

%go by subregion
%ensure table 1 and 2 have the same number of subregions
subregions_table1=unique(table1{1}.Subregion,'stable');
subregions_table2=unique(table2{1}.Subregion,'stable');

if length(subregions_table1)==length(subregions_table2)
    if any(ismember(subregions_table1,subregions_table2)==0)
        error('tables do not contain the same subregions')
    end
    subregions=subregions_table1;
elseif length(subregions_table1)<length(subregions_table2)
    if any(ismember(subregions_table1,subregions_table2)==0)
        error('tables do not contain the same subregions')
    end
    subregions=subregions_table1;
elseif length(subregions_table1)>length(subregions_table2)
    if any(ismember(subregions_table2,subregions_table1)==0)
        error('tables do not contain the same subregions')
    end
    subregions=subregions_table2;
else
    error('Please check fields of input tables.')
end


for i=1:length(subregions)
    %table 1 area count
    table1_areas=[];
    for j=1:length(table1)
        temp_area=(table1{j}.("Spike Area")(table1{j}.Subregion==subregions(i)));
        temp_area=[temp_area{:}];
        table1_areas=[table1_areas,temp_area/table1_duration];
    end
    table1_areas_n(i)=length(table1_areas);
    table1_mean(i)=mean(table1_areas);
    
    table2_areas=[];
    %table 2 area count
    for j=1:length(table2)
        temp_area=(table2{j}.("Spike Area")(table2{j}.Subregion==subregions(i)));
        temp_area=[temp_area{:}];
        table2_areas=[table2_areas,temp_area/table2_duration];
    end
    table2_areas_n(i)=length(table2_areas);
    table2_mean(i)=mean(table2_areas);
    
    [h,p,ci,stats]=ttest2(table1_areas,table2_areas);
    h_test(i)=h;
    p_test(i)=p;
    ci_test(i,:)=ci;
    stats_test{i}=stats;
    std_err_1(i)=std(table1_areas)/sqrt(length(table1_areas));
    std_err_2(i)=std(table2_areas)/sqrt(length(table2_areas));
end


x=categorical(cellstr(subregions));
x=reordercats(x,subregions);
vals=[table1_mean;table2_mean];
figure
b=bar(x,vals);
hold on

err=[std_err_1;std_err_2];

[ngroups,nbars]=size(vals);
x = nan(nbars, ngroups);
for i = 1:ngroups
    x(:,i) = b(i).XEndPoints;
end
errorbar(x',vals,err,'k','linestyle','none');
hold off
% groupwidth = min(0.8, nbars/(nbars + 1.5));
% for i=1:nbars
%     x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
%     errorbar(x,vals(:,i),err(:,i),'k.')
% end
xlabel("Subregions")
ylabel("Spike Frequency (Hz)")

percent_diff=((vals(2,:)-vals(1,:))./vals(1,:))*100;
end