function [conduction_time,max_area,direction,indeterminate,num_peaks]=sum_area_210504(bincountvec, bincount_ISI,bar_graph, indeterminate,peakprom_mult,std_num)

% Find feed forward and feedback peaks
% Assumes even number of elements in bincount vec and bin count ISI
% Make sure this function goes after get_spike_delay_hist to ensure that a
% figure is currently open
% This function also colors the histogram to show which bins are being
% summed

%Added multiple peak detection 5/4/21

hold on
std_from_mean=std_num*std(bincountvec);
yline(std_from_mean)

peakprom=peakprom_mult*max(bincountvec);

[max_pk_fb,max_fb_locs]=findpeaks(bincountvec(1:length(bincountvec)/2),bincount_ISI(1:length(bincount_ISI)/2),'MinPeakProminence',peakprom,'MinPeakHeight',std_from_mean,'MinPeakDistance',0.14);
[max_pk_ff,max_ff_locs]=findpeaks(bincountvec((length(bincountvec)/2)+1:end),bincount_ISI((length(bincount_ISI)/2)+1:end),'MinPeakProminence',peakprom,'MinPeakHeight',std_from_mean,'MinPeakDistance',.14);
plot(max_fb_locs,max_pk_fb,'x','MarkerSize',12,'MarkerEdgeColor','k','Color','k')
plot(max_ff_locs,max_pk_ff,'x','MarkerSize',12,'MarkerEdgeColor','k','Color','k')
num_peaks=length(max_fb_locs)+length(max_ff_locs);

% Choose viable peaks based on percentage spikes captures
hist_width_to_sum_ff=2;
hist_width_to_sum_fb=2;
% 
% if length(max_fb_locs)>1
%     hist_width_to_sum_fb=1;
% end
% 
% if length(max_ff_locs)>1
%     hist_width_to_sum_ff=1;
% end


fb_pks=[];
fb_locs=[];
fb_areas=[];
left_bound_fb=[];
right_bound_fb=[];
lfb_bound=[];
rfb_bound=[];
fb_cdt=[];
for i=1:length(max_pk_fb)
        % Left most
    if find(bincount_ISI==max_fb_locs(i))-hist_width_to_sum_fb>0
        left_bound_fb=[left_bound_fb,find(bincount_ISI==max_fb_locs(i))-hist_width_to_sum_fb];
    else
        left_bound_fb=[left_bound_fb,1];
    end

    % Right most
    size_vec=length(bincountvec)/2;
    if find(bincount_ISI==max_fb_locs(i))+hist_width_to_sum_fb<size_vec
        right_bound_fb=[right_bound_fb,find(bincount_ISI==max_fb_locs(i))+hist_width_to_sum_fb];
    else
        right_bound_fb=[right_bound_fb,size_vec];
    end
    
    % Get area
    area=sum(bincountvec(left_bound_fb(end):right_bound_fb(end)));
    total_area=sum(bincountvec(1:length(bincountvec)/2));
    
    %Check area is at least 80% of spikes
    %if area>=total_area*.5
        fb_pks=[fb_pks,find(bincount_ISI==max_fb_locs(i))];
        fb_locs=[fb_locs,find(bincount_ISI==max_fb_locs(i))];
        fb_cdt=[fb_cdt,bincount_ISI(bincount_ISI==max_fb_locs(i))];
        fb_areas=[fb_areas,area];
        lfb_bound=[lfb_bound,(left_bound_fb(i))];
        rfb_bound=[rfb_bound,(right_bound_fb(i))];
    %end
    
end

ff_pks=[];
ff_locs=[];
ff_areas=[];
left_bound_ff=[];
right_bound_ff=[];
lff_bound=[];
rff_bound=[];
ff_cdt=[]; %conduction_time
for i=1:length(max_pk_ff)
    size_vec=length(bincountvec);
        % Left most
    if find(bincount_ISI==max_ff_locs(i))-hist_width_to_sum_ff>(size_vec/2)+1
        left_bound_ff=[left_bound_ff,find(bincount_ISI==max_ff_locs(i))-hist_width_to_sum_ff];
    else
        left_bound_ff=[left_bound_ff,length(bincountvec/2)+1];
    end

    % Right most
    
    if find(bincount_ISI==max_ff_locs(i))+hist_width_to_sum_ff<size_vec
        right_bound_ff=[right_bound_ff,find(bincount_ISI==max_ff_locs(i))+hist_width_to_sum_ff];
    else
        right_bound_ff=[right_bound_ff,size_vec];
    end
    
    % Get area
    area=sum(bincountvec(left_bound_ff(end):right_bound_ff(end)));
    total_area=sum(bincountvec(size_vec/2+1:size_vec));
    
    %Check area is at least 80% of spikes
    %if area>=total_area*.5
        ff_pks=[ff_pks,max_pk_ff(i)];
        ff_locs=[ff_locs,find(bincount_ISI==max_ff_locs(i))];
        ff_areas=[ff_areas,area];
        ff_cdt=[ff_cdt,bincount_ISI(bincount_ISI==max_ff_locs(i))];
        lff_bound=[lff_bound,left_bound_ff(i)];
        rff_bound=[rff_bound,right_bound_ff(i)];
    %end
    
end



%Concatonate and sort peaks, choose top four
%fb first then ff
pks=[fb_pks,ff_pks];
locs=[fb_locs,ff_locs];
areas=[fb_areas,ff_areas];
cdt=[fb_cdt,ff_cdt];
l_bound=[lfb_bound,lff_bound];
r_bound=[rfb_bound,rff_bound];

%no peaks found case
if isempty(pks)
    max_area=-1;
    pks=0;
    direction=-1;
    conduction_time=0;
    indeterminate=indeterminate+1;
    return
end

if length(pks)>4
    [~,I]=sort(pks,'descend');
    pks=pks(sort(I(1:4)));
    locs=locs(sort(I(1:4)));
    areas=areas(sort(I(1:4)));
    cdt=cdt(sort(I(1:4)));
    l_bound=l_bound(sort(I(1:4)));
    r_bound=r_bound(sort(I(1:4)));
end

%bins to next peak
bins2peak=[];
for i=1:length(pks)
    if i==1 && length(pks)>1
        bins2peak(i,1)=Inf;
        bins2peak(i,2)=locs(i+1)-locs(i);
    elseif i==length(pks) && length(pks)>1
        bins2peak(i,2)=Inf;
        bins2peak(i,1)=locs(i)-locs(i-1);
    elseif length(pks)==1
        bins2peak(i,2)=Inf;
        bins2peak(i,1)=Inf;
    else
        bins2peak(i,1)=locs(i)-locs(i-1);
        bins2peak(i,2)=locs(i+1)-locs(i);
    end
end

colors=[1,0,0;.5,0,.5;0,1,0;1,1,0];
max_area=zeros(1,length(pks));
direction=max_area;
conduction_time=max_area;
for i=1:length(pks)
    
%     if i<length(pks) && length(pks)>1
%         if locs(i+1)-locs(i)<3
%             bin_locs_overlap=[l_bound(i+1):r_bound(i)];
%             area_overlap=sum(bincountvec(binc_locs_overlap));
%             half_over_1=round(area_overlap)/2;
%             half_over_2=area_overlap-half_over_1;
%             max_area(i)=sum(bincountvec(lbound(i):bin_locs_overlap(1)-1))+half_over_1;
%             max_area(i+1)=sum(bincountvec(bin_locs_overlap+1:locs(i+1)))+half_over_2;
%             if locs(i+2)-locs(i+1)<3
%                 
%             else
%                 max_area(i+1)=max_area(i+1)+sum(bincountvec(locs(i+1)+1:r_bound(i+1)));
%             end
    if i==1 && length(pks)>1
        if bins2peak(i,2)<5
            bin_locs_overlap=[l_bound(i+1):r_bound(i)];
            area_overlap=sum(bincountvec(bin_locs_overlap));
            half_over=round(area_overlap/2);
            max_area(i)=max_area(i)+sum(bincountvec(l_bound(i):bin_locs_overlap(1)-1))+half_over;
            for j=l_bound(i):bin_locs_overlap(1)
                bar_graph.CData(j,:)=[1 0 0];
            end
        else
            max_area(i)=areas(i);
            for j=l_bound(i):r_bound(i)
                bar_graph.CData(j,:)=[1 0 0];
            end
        end
        conduction_time(i)=cdt(i);
        if cdt(i)>0
            direction(i)=1;
        else
            direction(i)=0;
        end
    elseif i==length(pks) && length(pks)>1
        if bins2peak(i,1)<5
            bin_locs_overlap=[l_bound(i):r_bound(i-1)];
            area_overlap=sum(bincountvec(bin_locs_overlap));
            half_over=round(area_overlap/2);
            max_area(i)=max_area(i)+sum(bincountvec(bin_locs_overlap(end)+1:r_bound(i)))+half_over;
            bar_graph.FaceColor='flat';
            for j=bin_locs_overlap(end):r_bound(i)
                bar_graph.CData(j,:)=[1 0 0];
            end
            if cdt(i)>0
                direction(i)=1;
            else
                direction(i)=0;
            end
        else
            max_area(i)=areas(i);
            for j=l_bound(i):r_bound(i)
                bar_graph.CData(j,:)=[1 0 0];
            end
            if cdt(i)>0
                direction(i)=1;
            else
                direction(i)=0;
            end
        end
        conduction_time(i)=cdt(i);
        if cdt(i)>0
            direction(i)=1;
        else
            direction(i)=0;
        end
    elseif ~isempty(bins2peak)
        if bins2peak(i,2)<5 && bins2peak(i,1)<5
            %overlap on right
            bin_locs_overlap1=[l_bound(i+1):r_bound(i)];
            area_overlap=sum(bincountvec(bin_locs_overlap1));
            half_over=round(area_overlap/2);
            max_area(i)=max_area(i)+half_over;
            %overlap on left
            bin_locs_overlap2=[l_bound(i):r_bound(i-1)];
            area_overlap=sum(bincountvec(bin_locs_overlap2));
            half_over=round(area_overlap/2);
            max_area(i)=max_area(i)+sum(bincountvec(bin_locs_overlap2(end)+1:bin_locs_overlap1(end)-1))+half_over;
            for j=bin_locs_overlap2(end):bin_locs_overlap1(end)
                bar_graph.CData(j,:)=[1 0 0];
            end
        else
            max_area(i)=areas(i);
            for j=l_bound(i):r_bound(i)
                bar_graph.CData(j,:)=[1 0 0];
            end
        end
        conduction_time(i)=cdt(i);
        if cdt(i)>0
            direction(i)=1;
        else
            direction(i)=0;
        end
    end
end


hold off
end