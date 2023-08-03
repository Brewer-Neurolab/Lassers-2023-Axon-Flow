function [areas,total_bins,overlap_bins, total_area, overlap_area]=...
    ISI_conduction_time_peak_finding_210909...
    (stim,stim_dir,times_struct_stim,stim_duration,parent_dir,if_cw)
% Testing peak finding capabilities across multiple units for elimination
% of zeros
% Fixed saving issue and added bin count for percent overlap bins 7/7/21

total_bins=0;
overlap_bins=0;
total_area=0;
overlap_area=0;

if if_cw
    load(strcat(parent_dir,'\','matching_table_cw'))
elseif ~if_cw
    load(strcat(parent_dir,'\','matching_table_ccw'))
else
    error("Please define if the array is cw or ccw with 1 and 0.")
end

status_stim=stim.status;
stim_tunnel=stim.current_Tunnel_Pair;
ISI=stim.ISI;

insert_peak_table=table(matching_table.tunnel_pairs);

%%
% find peaks for each
%stimulated conduction times
peaks_stim=cell(length(insert_peak_table.Var1),1);
heights_stim=cell(length(insert_peak_table.Var1),1);
stim_stat=cell(length(insert_peak_table.Var1),1);
[stim_rows,~]=size(status_stim);

peakprom_mult=0.12;
std_num=1.9;

for stat_idx=1:stim_rows
    figure
    hold on
    bincountvec=status_stim(stat_idx,:);
    bar(ISI,bincountvec)
    std_from_mean=std_num*std(bincountvec);
    yline(std_from_mean)

    peakprom=peakprom_mult*max(bincountvec);

    [max_pk_fb,max_fb_locs]=findpeaks(bincountvec(1:length(bincountvec)/2),ISI(1:length(ISI)/2),'MinPeakProminence',peakprom,'MinPeakHeight',std_from_mean,'MinPeakDistance',0.14);
    [max_pk_ff,max_ff_locs]=findpeaks(bincountvec((length(bincountvec)/2)+1:end),ISI((length(ISI)/2)+1:end),'MinPeakProminence',peakprom,'MinPeakHeight',std_from_mean,'MinPeakDistance',.14);
    
    plot(max_fb_locs,max_pk_fb,'x','MarkerSize',12,'MarkerEdgeColor','k','Color','k')
    plot(max_ff_locs,max_pk_ff,'x','MarkerSize',12,'MarkerEdgeColor','k','Color','k')
    num_peaks=length(max_fb_locs)+length(max_ff_locs);
    peak_locs=[max_fb_locs,max_ff_locs];
    peak_heights=[max_pk_fb,max_pk_ff];
    
    for j=1:length(insert_peak_table.Var1)
        if contains(stim_tunnel{stat_idx},insert_peak_table.Var1(j))
            peaks_stim{j}=[peaks_stim{j},{peak_locs}];
            heights_stim{j}=[heights_stim{j},{peak_heights}];
            stim_stat{j}=[stim_stat{j},{bincountvec}];
%             stim_locs_ff{j}=[stim_locs_ff{j},{max_ff_locs}];
%             stim_locs_fb{j}=[stim_locs_fb{j},{max_fb_locs}];
%             stim_pks_ff{j}=[stim_pks_ff{j},{max_pk_ff}];
%             stim_pks_fb{j}=[stim_pks_fb{j},{max_pk_fb}];
        end
    end
    hold off
    close all
end



insert_peak_table=table(insert_peak_table.Var1...
    ,peaks_stim,heights_stim,stim_stat,...
    'VariableNames',{'tunnel_pairs',...
    'stim_peaks','stim_heights','stim_stat'});



%% Calculate areas
% Stim Areas

stim_locs=[insert_peak_table.stim_peaks{1:end}];
peak_heights=[insert_peak_table.stim_heights{1:end}];
for i=1:length(stim_locs)
    stim_fb_locs{i}=stim_locs{i}(stim_locs{i}<0);
    stim_ff_locs{i}=stim_locs{i}(stim_locs{i}>0);
    stim_fb_heights{i}=peak_heights{i}(stim_locs{i}<0);
    stim_ff_heights{i}=peak_heights{i}(stim_locs{i}>0);
end

[rw_stim,~]=size(status_stim);
current_Tunnel_Pair=stim_tunnel;
areas{2}.file_name=stim.file_name;
areas{2}.area=[];
areas{2}.conduction_time=[];
areas{2}.direction=[];
areas{2}.num_peaks=[];
areas{2}.channel_pair=[];
areas{2}.up_times=[];
areas{2}.down_times=[];

for stat_idx=1:rw_stim
    bincountvec=status_stim(stat_idx,:);
    b=bar(ISI,bincountvec);
    std_from_mean=std_num*std(bincountvec);
    yline(std_from_mean)
    [conduction_time_temp,area_temp,direction_temp,num_peaks_temp,up_times_stim_temp,down_times_stim_temp]=...
        ISI_sum_area_nonpeakstest_210903(bincountvec,ISI,b,stim_fb_locs{stat_idx},stim_ff_locs{stat_idx},stim_fb_heights{stat_idx},stim_ff_heights{stat_idx},times_struct_stim(stat_idx),stim_duration);
    
    if conduction_time_temp~=0
        set(gca,'fontsize',13)
        xlabel('<- Feedback            Conduction time (ms)         Feedforward->');
        ylabel('Spike Hz');
        title(char(current_Tunnel_Pair{stat_idx})+" Area:"+string(area_temp))
        saveas(gcf,strcat(stim_dir,'\',char(current_Tunnel_Pair{stat_idx})+string(fix(max(area_temp)))),'png')

        mult_tun_pair=[];
        for len_tunnel_pair=1:length(area_temp)
            mult_tun_pair=[mult_tun_pair;current_Tunnel_Pair{stat_idx}];
        end

        areas{2}.channel_pair=[areas{2}.channel_pair;string(mult_tun_pair)];
        areas{2}.area=[areas{2}.area;area_temp'];
        areas{2}.conduction_time=[areas{2}.conduction_time;abs(conduction_time_temp')];
        areas{2}.direction=[areas{2}.direction;direction_temp'];
        areas{2}.num_peaks=[areas{2}.num_peaks;num_peaks_temp'];
        %indeterminate=indeterminate+indeterminate_temp;
        areas{2}.up_times=[areas{2}.up_times;up_times_stim_temp'];
        areas{2}.down_times=[areas{2}.down_times;down_times_stim_temp'];
    else

        areas{2}.channel_pair=[areas{2}.channel_pair;string(current_Tunnel_Pair{stat_idx})];
        areas{2}.area=[areas{1}.area;-1];
        areas{2}.conduction_time=[areas{2}.conduction_time;conduction_time_temp];
        areas{2}.direction=[areas{2}.direction;direction_temp];
        areas{2}.up_times=[areas{2}.up_times;{[]}];
        areas{2}.down_times=[areas{2}.down_times;{[]}];
    end
end
end