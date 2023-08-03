function [areas,total_bins,overlap_bins, total_area, overlap_area]=...
    ISI_interunit_peak_finding_211108...
    (spon,stim,spon_dir,stim_dir,times_struct_spon,times_struct_stim,spon_duration,...
    stim_duration,parent_dir,if_cw)
% Testing peak finding capabilities across multiple units for elimination
% of zeros
% Ensure status and tunnel pairs are loaded spon and post stim
% Fixed saving issue and added bin count for percent overlap bins 7/7/21
% Fixed error saving spikes in dead time zone 11/4/21 due to smaller nmi
% threshold previously

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

status_spon=spon.status;
status_stim=stim.status;
spon_tunnel=spon.current_Tunnel_Pair;
stim_tunnel=stim.current_Tunnel_Pair;
ISI=spon.ISI;
time_param_spon=times_struct_spon;
time_param_stim=times_struct_stim;

%Need to convert time_param_xxxx back to times struct 11/8/21

insert_peak_table=table(matching_table.tunnel_pairs);

%%
% find peaks for each

%spontaneous conduction times
peaks_spon=cell(length(insert_peak_table.Var1),1);
heights_spon=cell(length(insert_peak_table.Var1),1);
spon_stat=cell(length(insert_peak_table.Var1),1);
% spon_locs_ff=cell(length(insert_peak_table.Var1),1);
% spon_locs_fb=cell(length(insert_peak_table.Var1),1);
% spon_pks_ff=cell(length(insert_peak_table.Var1),1);
% spon_pks_fb=cell(length(insert_peak_table.Var1),1);
peakprom_mult=0.12;
std_num=1.9;
[spon_rows,~]=size(status_spon);

min_difference_1=0.2;
min_difference_2=0.12;
%max_difference=1.0;
bin_width=diff(ISI);
bin_width=bin_width(1);
%ref_per=1.5;

% IMPORTANT NOTE: HARD CODE UNSURE IF THIS WORKS ACROSS ALL VARIATIONS,
% Original preserves original conduction time hist, first two bins on
% either side of 0 are taken out, edges correspond to minimum difference at
% my current resolution
for stat_idx=1:spon_rows
    figure
    hold on
    original=status_spon(stat_idx,:);
    bincountvec=status_spon(stat_idx,:);
    bincountvec(and(ISI<((2*bin_width)+(0.5*bin_width)),ISI>-((2*bin_width)+(0.5*bin_width))))=0;
    bar(ISI,bincountvec)
    std_from_mean=std_num*std(bincountvec);
    xticks(ISI)
    yline(std_from_mean)

    peakprom=peakprom_mult*max(bincountvec);
    
    %adjust min peak distance for two bins
    [max_pk_fb,max_fb_locs]=findpeaks(bincountvec(1:length(bincountvec)/2),ISI(1:length(ISI)/2),'MinPeakProminence',peakprom,'MinPeakHeight',std_from_mean,'MinPeakDistance',0.23);
    [max_pk_ff,max_ff_locs]=findpeaks(bincountvec((length(bincountvec)/2)+1:end),ISI((length(ISI)/2)+1:end),'MinPeakProminence',peakprom,'MinPeakHeight',std_from_mean,'MinPeakDistance',.23);
    
    
    if any(ismember(round(max_fb_locs,5),round(-(3*bin_width),5)))
        bincountvec(and(ISI<0,ISI>-((2*bin_width)+(0.5*bin_width))))=original(and(ISI<0,ISI>-((2*bin_width)+(0.5*bin_width))));
        if original(round(ISI,5)==round(-(2*bin_width),5))>original(round(ISI,5)==round(-(3*bin_width),5))
            max_pk_fb(ismember(round(max_fb_locs,5),round(-(3*bin_width),5)))=[];
            max_fb_locs(ismember(round(max_fb_locs,5),round(-(3*bin_width),5)))=[];
            delete_ISI=and(time_param_spon(stat_idx).ISI>-min_difference_1,time_param_spon(stat_idx).ISI<0);
            time_param_spon(stat_idx).ISI(delete_ISI)=NaN;
            time_param_spon(stat_idx).up(delete_ISI)=NaN;
            time_param_spon(stat_idx).down(delete_ISI)=NaN;
            bincountvec(and(ISI<0,ISI>-((2*bin_width)+(0.5*bin_width))))=0;
            status_spon(stat_idx,:)=bincountvec;
        else
            delete_ISI=and(time_param_spon(stat_idx).ISI>-min_difference_2,time_param_spon(stat_idx).ISI<0);
            time_param_spon(stat_idx).ISI(delete_ISI)=NaN;
            time_param_spon(stat_idx).up(delete_ISI)=NaN;
            time_param_spon(stat_idx).down(delete_ISI)=NaN;
            bincountvec(and(ISI<0,ISI>-((1*bin_width)+(0.5*bin_width))))=0;
            status_spon(stat_idx,:)=bincountvec;
        end
    else
        delete_ISI=and(time_param_spon(stat_idx).ISI>-min_difference_1,time_param_spon(stat_idx).ISI<0);
        time_param_spon(stat_idx).ISI(delete_ISI)=NaN;
        time_param_spon(stat_idx).up(delete_ISI)=NaN;
        time_param_spon(stat_idx).down(delete_ISI)=NaN;
        bincountvec(and(ISI<0,ISI>-((2*bin_width)+(0.5*bin_width))))=0;
        status_spon(stat_idx,:)=bincountvec;
    end
    if any(ismember(round(max_ff_locs,5),round((3*bin_width),5)))
        bincountvec(and(ISI>0,ISI<((2*bin_width)+(0.5*bin_width))))=original(and(ISI>0,ISI<((2*bin_width)+(0.5*bin_width))));
        if original(round(ISI,5)==round((2*bin_width),5))>original(round(ISI,5)==round((3*bin_width),5))
            max_pk_ff(ismember(round(max_ff_locs,5),round((3*bin_width),5)))=[];
            max_ff_locs(ismember(round(max_ff_locs,5),round((3*bin_width),5)))=[];
            delete_ISI=and(time_param_spon(stat_idx).ISI<min_difference_1,time_param_spon(stat_idx).ISI>0);
            time_param_spon(stat_idx).ISI(delete_ISI)=NaN;
            time_param_spon(stat_idx).up(delete_ISI)=NaN;
            time_param_spon(stat_idx).down(delete_ISI)=NaN;
            bincountvec(and(ISI<0,ISI>-((2*bin_width)+(0.5*bin_width))))=0;
            status_spon(stat_idx,:)=bincountvec;
        else
            delete_ISI=and(time_param_spon(stat_idx).ISI<min_difference_2,time_param_spon(stat_idx).ISI>0);
            time_param_spon(stat_idx).ISI(delete_ISI)=NaN;
            time_param_spon(stat_idx).up(delete_ISI)=NaN;
            time_param_spon(stat_idx).down(delete_ISI)=NaN;
            bincountvec(and(ISI>0,ISI<((1*bin_width)+(0.5*bin_width))))=0;
            status_spon(stat_idx,:)=bincountvec;
        end
    else
        delete_ISI=and(time_param_spon(stat_idx).ISI<min_difference_1,time_param_spon(stat_idx).ISI>0);
        time_param_spon(stat_idx).ISI(delete_ISI)=NaN;
        time_param_spon(stat_idx).up(delete_ISI)=NaN;
        time_param_spon(stat_idx).down(delete_ISI)=NaN;
        bincountvec(and(ISI>0,ISI<((2*bin_width)+(0.5*bin_width))))=0;
        status_spon(stat_idx,:)=bincountvec;
    end
    
    spon_pk_fb{stat_idx}=max_pk_fb;
    spon_pk_ff{stat_idx}=max_pk_ff;
    
    plot(max_fb_locs,max_pk_fb,'x','MarkerSize',12,'MarkerEdgeColor','k','Color','k')
    plot(max_ff_locs,max_pk_ff,'x','MarkerSize',12,'MarkerEdgeColor','k','Color','k')
    num_peaks=length(max_fb_locs)+length(max_ff_locs);
    peak_locs=[max_fb_locs,max_ff_locs];
    peak_heights=[max_pk_fb,max_pk_ff];
    
    for j=1:length(insert_peak_table.Var1)
        if contains(spon_tunnel{stat_idx},insert_peak_table.Var1(j))
            peaks_spon{j}=[peaks_spon{j},{peak_locs}];
            heights_spon{j}=[heights_spon{j},{peak_heights}];
            spon_stat{j}=[spon_stat{j},{bincountvec}];
%             spon_locs_ff{j}=[spon_locs_ff{j},{max_ff_locs}];
%             spon_locs_fb{j}=[spon_locs_fb{j},{max_fb_locs}];
%             spon_pks_ff{j}=[spon_pks_ff{j},{max_pk_ff}];
%             spon_pks_fb{j}=[spon_pks_fb{j},{max_pk_fb}];
        end
    end
    hold off
    close all
end

%return back to times struct
times_struct_spon=time_param_spon;

%stimulated conduction times
peaks_stim=cell(length(insert_peak_table.Var1),1);
heights_stim=cell(length(insert_peak_table.Var1),1);
stim_stat=cell(length(insert_peak_table.Var1),1);
% stim_locs_ff=cell(length(insert_peak_table.Var1),1);
% stim_locs_fb=cell(length(insert_peak_table.Var1),1);
% stim_pks_ff=cell(length(insert_peak_table.Var1),1);
% stim_pks_fb=cell(length(insert_peak_table.Var1),1);
[stim_rows,~]=size(status_stim);

for stat_idx=1:stim_rows
    figure
    hold on
    original=status_stim(stat_idx,:);
    bincountvec=status_stim(stat_idx,:);
    bincountvec(and(ISI<((2*bin_width)+(0.5*bin_width)),ISI>-((2*bin_width)+(0.5*bin_width))))=0;
    bar(ISI,bincountvec)
    xticks(ISI)
    std_from_mean=std_num*std(bincountvec);
    yline(std_from_mean)

    peakprom=peakprom_mult*max(bincountvec);

    [max_pk_fb,max_fb_locs]=findpeaks(bincountvec(1:length(bincountvec)/2),ISI(1:length(ISI)/2),'MinPeakProminence',peakprom,'MinPeakHeight',std_from_mean,'MinPeakDistance',0.24);
    [max_pk_ff,max_ff_locs]=findpeaks(bincountvec((length(bincountvec)/2)+1:end),ISI((length(ISI)/2)+1:end),'MinPeakProminence',peakprom,'MinPeakHeight',std_from_mean,'MinPeakDistance',.24);
    
    if any(ismember(round(max_fb_locs,5),round(-(3*bin_width),5)))
        bincountvec(and(ISI<0,ISI>-((2*bin_width)+(0.5*bin_width))))=original(and(ISI<0,ISI>-((2*bin_width)+(0.5*bin_width))));
        if original(round(ISI,5)==round(-(2*bin_width),5))>original(round(ISI,5)==round(-(3*bin_width),5))
            max_pk_fb(ismember(round(max_fb_locs,5),round(-(3*bin_width),5)))=[];
            max_fb_locs(ismember(round(max_fb_locs,5),round(-(3*bin_width),5)))=[];
            delete_ISI=and(time_param_stim(stat_idx).ISI>-min_difference_1,time_param_stim(stat_idx).ISI<0);
            time_param_stim(stat_idx).ISI(delete_ISI)=NaN;
            time_param_stim(stat_idx).up(delete_ISI)=NaN;
            time_param_stim(stat_idx).down(delete_ISI)=NaN;
            bincountvec(and(ISI<0,ISI>-((2*bin_width)+(0.5*bin_width))))=0;
            status_stim(stat_idx,:)=bincountvec;
        else
            delete_ISI=and(time_param_stim(stat_idx).ISI>-min_difference_2,time_param_stim(stat_idx).ISI<0);
            time_param_stim(stat_idx).ISI(delete_ISI)=NaN;
            time_param_stim(stat_idx).up(delete_ISI)=NaN;
            time_param_stim(stat_idx).down(delete_ISI)=NaN;
            bincountvec(and(ISI<0,ISI>-((1*bin_width)+(0.5*bin_width))))=0;
            status_stim(stat_idx,:)=bincountvec;
        end
    else
        delete_ISI=and(time_param_stim(stat_idx).ISI>-min_difference_1,time_param_stim(stat_idx).ISI<0);
        time_param_stim(stat_idx).ISI(delete_ISI)=NaN;
        time_param_stim(stat_idx).up(delete_ISI)=NaN;
        time_param_stim(stat_idx).down(delete_ISI)=NaN;
        bincountvec(and(ISI<0,ISI>-((2*bin_width)+(0.5*bin_width))))=0;
        status_stim(stat_idx,:)=bincountvec;
    end
    if any(ismember(round(max_ff_locs,5),round((3*bin_width),5)))
        bincountvec(and(ISI>0,ISI<((2*bin_width)+(0.5*bin_width))))=original(and(ISI>0,ISI<((2*bin_width)+(0.5*bin_width))));
        if original(round(ISI,5)==round((2*bin_width),5))>original(round(ISI,5)==round((3*bin_width),5))
            max_pk_ff(ismember(round(max_ff_locs,5),round((3*bin_width),5)))=[];
            max_ff_locs(ismember(round(max_ff_locs,5),round((3*bin_width),5)))=[];
            delete_ISI=and(time_param_stim(stat_idx).ISI<min_difference_1,time_param_stim(stat_idx).ISI>0);
            time_param_stim(stat_idx).ISI(delete_ISI)=NaN;
            time_param_stim(stat_idx).up(delete_ISI)=NaN;
            time_param_stim(stat_idx).down(delete_ISI)=NaN;
            bincountvec(and(ISI<0,ISI>-((2*bin_width)+(0.5*bin_width))))=0;
            status_stim(stat_idx,:)=bincountvec;
        else
            delete_ISI=and(time_param_stim(stat_idx).ISI<min_difference_2,time_param_stim(stat_idx).ISI>0);
            time_param_stim(stat_idx).ISI(delete_ISI)=NaN;
            time_param_stim(stat_idx).up(delete_ISI)=NaN;
            time_param_stim(stat_idx).down(delete_ISI)=NaN;
            bincountvec(and(ISI>0,ISI<((1*bin_width)+(0.5*bin_width))))=0;
            status_stim(stat_idx,:)=bincountvec;
        end
    else
        delete_ISI=and(time_param_stim(stat_idx).ISI<min_difference_1,time_param_stim(stat_idx).ISI>0);
        time_param_stim(stat_idx).ISI(delete_ISI)=NaN;
        time_param_stim(stat_idx).up(delete_ISI)=NaN;
        time_param_stim(stat_idx).down(delete_ISI)=NaN;
        bincountvec(and(ISI<0,ISI>-((2*bin_width)+(0.5*bin_width))))=0;
        status_stim(stat_idx,:)=bincountvec;
    end
    
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

%return back to times struct
times_struct_stim=time_param_stim;

insert_peak_table=table(insert_peak_table.Var1,peaks_spon,heights_spon,spon_stat...
    ,peaks_stim,heights_stim,stim_stat,...
    'VariableNames',{'tunnel_pairs','spon_peaks','spon_heights','spon_stat',...
    'stim_peaks','stim_heights','stim_stat'});


%% Place Axons V4
%changed form 0.05 to 0.1 6/3/21

for tunnel_idx=1:length(insert_peak_table.tunnel_pairs)
    spon_peaks=cell2mat(insert_peak_table.spon_peaks{tunnel_idx});
    spon_heights=cell2mat(insert_peak_table.spon_heights{tunnel_idx});
    stim_peaks=cell2mat(insert_peak_table.stim_peaks{tunnel_idx});
    stim_heights=cell2mat(insert_peak_table.stim_heights{tunnel_idx});
    
    %sorting
    [spon_peaks,spon_I]=sort(spon_peaks);
    [stim_peaks,stim_I]=sort(stim_peaks);
    spon_heights=spon_heights(spon_I);
    stim_heights=stim_heights(stim_I);
    
    %spon side
    for spons=1:length(insert_peak_table.spon_peaks{tunnel_idx})
        for stims_to_check=1:length(stim_peaks)
            
            pos_peaks_dif=insert_peak_table.spon_peaks{tunnel_idx}{spons}(insert_peak_table.spon_peaks{tunnel_idx}{spons}>0);
            neg_peaks_dif=insert_peak_table.spon_peaks{tunnel_idx}{spons}(insert_peak_table.spon_peaks{tunnel_idx}{spons}<0);
            
            if ~any(abs(pos_peaks_dif-stim_peaks(stims_to_check))<=0.16001) && stim_peaks(stims_to_check)>0
                insert_peak_table.spon_peaks{tunnel_idx}{spons}=[insert_peak_table.spon_peaks{tunnel_idx}{spons},stim_peaks(stims_to_check)];
                insert_peak_table.spon_heights{tunnel_idx}{spons}=...
                    [insert_peak_table.spon_heights{tunnel_idx}{spons},...
                    (insert_peak_table.spon_stat{tunnel_idx}{spons}(ISI==stim_peaks(stims_to_check)))];
            end
            
            if ~any(abs(neg_peaks_dif-stim_peaks(stims_to_check))<=0.16001) && stim_peaks(stims_to_check)<0
                insert_peak_table.spon_peaks{tunnel_idx}{spons}=[insert_peak_table.spon_peaks{tunnel_idx}{spons},stim_peaks(stims_to_check)];
                insert_peak_table.spon_heights{tunnel_idx}{spons}=...
                    [insert_peak_table.spon_heights{tunnel_idx}{spons},...
                    (insert_peak_table.spon_stat{tunnel_idx}{spons}(ISI==stim_peaks(stims_to_check)))];
            end
            
            
        end
    end
    
    spon_peaks=cell2mat(insert_peak_table.spon_peaks{tunnel_idx});
    spon_heights=cell2mat(insert_peak_table.spon_heights{tunnel_idx});
    stim_peaks=cell2mat(insert_peak_table.stim_peaks{tunnel_idx});
    stim_heights=cell2mat(insert_peak_table.stim_heights{tunnel_idx});
    
    %sorting
    [spon_peaks,spon_I]=sort(spon_peaks);
    [stim_peaks,stim_I]=sort(stim_peaks);
    spon_heights=spon_heights(spon_I);
    stim_heights=stim_heights(stim_I);
    
    %Stim side
    for stims=1:length(insert_peak_table.stim_peaks{tunnel_idx})
        for spons_to_check=1:length(spon_peaks)
            
            pos_peaks_dif=insert_peak_table.stim_peaks{tunnel_idx}{stims}(insert_peak_table.stim_peaks{tunnel_idx}{stims}>0);
            neg_peaks_dif=insert_peak_table.stim_peaks{tunnel_idx}{stims}(insert_peak_table.stim_peaks{tunnel_idx}{stims}<0);

            if ~any(abs(pos_peaks_dif-spon_peaks(spons_to_check))<=0.16001) && spon_peaks(spons_to_check)>0
                insert_peak_table.stim_peaks{tunnel_idx}{stims}=[insert_peak_table.stim_peaks{tunnel_idx}{stims},spon_peaks(spons_to_check)];
                insert_peak_table.stim_heights{tunnel_idx}{stims}=...
                    [insert_peak_table.stim_heights{tunnel_idx}{stims},...
                    (insert_peak_table.stim_stat{tunnel_idx}{stims}(ISI==spon_peaks(spons_to_check)))];
            end
            
            if ~any(abs(neg_peaks_dif-spon_peaks(spons_to_check))<=0.16001) && spon_peaks(spons_to_check)<0
                insert_peak_table.stim_peaks{tunnel_idx}{stims}=[insert_peak_table.stim_peaks{tunnel_idx}{stims},spon_peaks(spons_to_check)];
                insert_peak_table.stim_heights{tunnel_idx}{stims}=...
                    [insert_peak_table.stim_heights{tunnel_idx}{stims},...
                    (insert_peak_table.stim_stat{tunnel_idx}{stims}(ISI==spon_peaks(spons_to_check)))];
            end
            
            
        end
    end
    
end

%% Calculate areas
spon_locs=[insert_peak_table.spon_peaks{1:end}];
peak_heights=[insert_peak_table.spon_heights{1:end}];
for i=1:length(spon_locs)
    spon_fb_locs{i}=spon_locs{i}(spon_locs{i}<0);
    spon_ff_locs{i}=spon_locs{i}(spon_locs{i}>0);
    spon_fb_heights{i}=peak_heights{i}(spon_locs{i}<0);
    spon_ff_heights{i}=peak_heights{i}(spon_locs{i}>0);
end

[rw_spon,~]=size(status_spon);
current_Tunnel_Pair=spon_tunnel;
areas{1}.file_name=spon.file_name;
areas{1}.area=[];
areas{1}.conduction_time=[];
areas{1}.direction=[];
areas{1}.num_peaks=[];
areas{1}.channel_pair=[];
areas{1}.up_times=[];
areas{1}.down_times=[];

for stat_idx=1:(rw_spon)
    bincountvec=status_spon(stat_idx,:);
    %bincountvec(and(ISI<((2*bin_width)+(0.5*bin_width)),ISI>-((2*bin_width)+(0.5*bin_width))))=0;
    bar(ISI,bincountvec)
    b=bar(ISI,bincountvec);
    std_from_mean=std_num*std(bincountvec);
    yline(std_from_mean)
%     [conduction_time_temp,area_temp,direction_temp,num_peaks_temp,up_times_spon_temp,down_times_spon_temp]=...
%         ISI_sum_area_nonpeakstest_210722(bincountvec,ISI,b,spon_fb_locs{stat_idx},spon_ff_locs{stat_idx},spon_fb_heights{stat_idx},spon_ff_heights{stat_idx}, times_struct_spon(stat_idx),spon_duration);
    [conduction_time_temp,area_temp,direction_temp,num_peaks_temp,up_times_spon_temp,down_times_spon_temp]=...
        ISI_sum_area_nonpeakstest_211025(bincountvec,ISI,b,spon_fb_locs{stat_idx},spon_ff_locs{stat_idx},spon_fb_heights{stat_idx},spon_ff_heights{stat_idx}, times_struct_spon(stat_idx),spon_duration);
    
    if conduction_time_temp~=0
        set(gca,'fontsize',13)
        xlabel('<- Feedback            Conduction time (ms)         Feedforward->');
        %ylabel('Spike Hz');
        ylabel('Spike Count')
        title(char(current_Tunnel_Pair{stat_idx})+" Total Area: "+string(sum(bincountvec))+" Area:"+string(area_temp))
        xticks(ISI)
        saveas(gcf,strcat(spon_dir,'\',char(current_Tunnel_Pair{stat_idx})+string(fix(max(area_temp)))),'png')

        mult_tun_pair=[];
        for len_tunnel_pair=1:length(area_temp)
            mult_tun_pair=[mult_tun_pair;current_Tunnel_Pair{stat_idx}];
        end

        areas{1}.channel_pair=[areas{1}.channel_pair;string(mult_tun_pair)];
        areas{1}.area=[areas{1}.area;area_temp'];
        areas{1}.conduction_time=[areas{1}.conduction_time;abs(conduction_time_temp')];
        areas{1}.direction=[areas{1}.direction;direction_temp'];
        areas{1}.num_peaks=[areas{1}.num_peaks;num_peaks_temp'];
        %indeterminate=indeterminate+indeterminate_temp;
        
        %remove times less than min diff
%         for axons=1:length(up_times_spon_temp)
%             %esure each value is not shared
%             if length(up_times_spon_temp{axons})>length(down_times_spon_temp{axons})
%                 [~,LocB]=ismembertol(up_times_spon_temp{axons},down_times_spon_temp{axons},1,'DataScale',1);
%                 up=unique(LocB);
%                 up_times_spon_temp{axons}=up_times_spon_temp{axons}(up);
%             elseif length(up_times_spon_temp{axons})<length(down_times_spon_temp{axons})
%                 [~,LocB]=ismembertol(down_times_spon_temp{axons},up_times_spon_temp{axons},1,'DataScale',1);
%                 down=unique(LocB);
%                 down_times_spon_temp{axons}=down_times_spon_temp{axons}(down);
%             end
%             
%             time_diff=abs([up_times_spon_temp{axons}-down_times_spon_temp{axons}]);
%             less_than_min=time_diff<min_difference;
%             up_times_spon_temp{axons}(less_than_min)=[];
%             down_times_spon_temp{axons}(less_than_min)=[];
%         end
        
        areas{1}.up_times=[areas{1}.up_times;up_times_spon_temp'];
        areas{1}.down_times=[areas{1}.down_times;down_times_spon_temp'];
    else

        areas{1}.channel_pair=[areas{1}.channel_pair;string(current_Tunnel_Pair{stat_idx})];
        areas{1}.area=[areas{1}.area;-1];
        areas{1}.conduction_time=[areas{1}.conduction_time;conduction_time_temp];
        areas{1}.direction=[areas{1}.direction;direction_temp];
        areas{1}.up_times=[areas{1}.up_times;{[]}];
        areas{1}.down_times=[areas{1}.down_times;{[]}];
    end
end

%% Stim Areas

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
    %bincountvec(and(ISI<((2*bin_width)+(0.5*bin_width)),ISI>-((2*bin_width)+(0.5*bin_width))))=0;
    b=bar(ISI,bincountvec);
    std_from_mean=std_num*std(bincountvec);
    yline(std_from_mean)
    [conduction_time_temp,area_temp,direction_temp,num_peaks_temp,up_times_stim_temp,down_times_stim_temp]=...
        ISI_sum_area_nonpeakstest_211025(bincountvec,ISI,b,stim_fb_locs{stat_idx},stim_ff_locs{stat_idx},stim_fb_heights{stat_idx},stim_ff_heights{stat_idx},times_struct_stim(stat_idx),stim_duration);
    
    if conduction_time_temp~=0
        set(gca,'fontsize',13)
        xlabel('<- Feedback            Conduction time (ms)         Feedforward->');
        %ylabel('Spike Hz');
        ylabel('Spike Count')
        title(char(current_Tunnel_Pair{stat_idx})+" Total Area: "+string(sum(bincountvec))+" Area:"+string(area_temp))
        xticks(ISI)
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
        %remove times less than min diff
%         for axons=1:length(up_times_stim_temp)
%             %esure each value is not shared
%             if length(up_times_stim_temp{axons})>length(down_times_stim_temp{axons})
%                 [~,LocB]=ismembertol(up_times_stim_temp{axons},down_times_stim_temp{axons},1,'DataScale',1);
%                 up=unique(LocB);
%                 up_times_stim_temp{axons}=up_times_stim_temp{axons}(up);
%             elseif length(up_times_stim_temp{axons})<length(down_times_stim_temp{axons})
%                 [~,LocB]=ismembertol(down_times_stim_temp{axons},up_times_stim_temp{axons},1,'DataScale',1);
%                 down=unique(LocB);
%                 down_times_stim_temp{axons}=down_times_stim_temp{axons}(down);
%             end
%             
%             time_diff=abs([up_times_stim_temp{axons}-down_times_stim_temp{axons}]);
%             less_than_min=time_diff<min_difference;
%             up_times_stim_temp{axons}(less_than_min)=[];
%             down_times_stim_temp{axons}(less_than_min)=[];
%         end
%         
        areas{2}.up_times=[areas{2}.up_times;up_times_stim_temp'];
        areas{2}.down_times=[areas{2}.down_times;down_times_stim_temp'];
     else

        areas{2}.channel_pair=[areas{2}.channel_pair;string(current_Tunnel_Pair{stat_idx})];
        areas{2}.area=[areas{2}.area;-1];
        areas{2}.conduction_time=[areas{2}.conduction_time;conduction_time_temp];
        areas{2}.direction=[areas{2}.direction;direction_temp];
        areas{2}.up_times=[areas{2}.up_times;{[]}];
        areas{2}.down_times=[areas{2}.down_times;{[]}];
    end
end
end