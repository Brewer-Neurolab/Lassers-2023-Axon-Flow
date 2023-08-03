function [conduction_time,max_area,direction,num_peaks,up_times,down_times]=...
    ISI_sum_area_nonpeakstest_211025(bincountvec, bincount_ISI,bar_graph,max_fb_locs,max_ff_locs,max_pk_fb,max_pk_ff,times_struct,duration)

% Find feed forward and feedback peaks
% Assumes even number of elements in bincount vec and bin count ISI
% Make sure this function goes after get_spike_delay_hist to ensure that a
% figure is currently open
% This function also colors the histogram to show which bins are being
% summed

% Added multiple peak detection 5/4/21
% Fixed color error and distribution of spikes to ratio of summation of
% peaks, considers all peaks not just top 4 6/21/21
% Counts bins summed and bins overlapped 7/7/21
% Fixed saving spike times instead of conduction times 211025

hold on

total_bins=0;
overlap_bins=0;
total_spikes=0;
overlap_spikes=0;

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



%Concatonate and sort peaks
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
    up_times=[];
    down_times=[];
    hold off
    %indeterminate=indeterminate+1;
    return
end

%Sort peaks and corresponding areas
[locs,I]=sort(locs);
pks=pks(I);
areas=areas(I);
cdt=cdt(I);
l_bound=l_bound(I);
r_bound=r_bound(I);

% Only consider the tallest 4 peaks
% if length(pks)>4
%     [~,I]=sort(pks,'descend');
%     pks=pks(sort(I(1:4)));
%     locs=locs(sort(I(1:4)));
%     areas=areas(sort(I(1:4)));
%     cdt=cdt(sort(I(1:4)));
%     l_bound=l_bound(sort(I(1:4)));
%     r_bound=r_bound(sort(I(1:4)));
% end

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

% bins2peak not lined up with lbound and rbound, additionally negatives
% create errors check in morning

colors=[1,0,0;.5,0,.5;0,1,0;1,1,0];
max_area=zeros(1,length(pks));
direction=max_area;
conduction_time=max_area;
% Overlap no longer exists with ISI, variable names kept the same for
% convinience 

for i=1:length(pks)
    
    %Left most peak considered, right handed overlap
    if i==1 && length(pks)>1
        if bins2peak(i,2)<5 && mod(bins2peak(i,2)-1,2)==0
            bin_locs_overlap=[l_bound(i+1):r_bound(i)];
            isolate_bins=bin_locs_overlap(1:length(bin_locs_overlap)/2);
            isolate_spike_count=sum(bincountvec(isolate_bins));
            max_area(i)=max_area(i)+sum(bincountvec(l_bound(i):bin_locs_overlap(1)-1))+isolate_spike_count;
            total_spikes=total_spikes+max_area(i);
            bar_graph.FaceColor='flat';
            for j=l_bound(i):bin_locs_overlap(1)
                bar_graph.CData(j,:)=[1 0 0];
            end
            %flag spikes
            bins=[l_bound(i):bin_locs_overlap(1)-1,isolate_bins];
            if cdt(i)<0
                edge_bins=times_struct.edges([bins,bins(end)+1]);
            elseif cdt(i)>0
                edge_bins=times_struct.edges([bins(1)-1,bins]);
            else
                error('Conduction time not supported.')
            end
            
            ISI_bins=(and(times_struct.ISI>=edge_bins(1),times_struct.ISI<edge_bins(end)));
            up_times{i}=times_struct.up(ISI_bins);
            down_times{i}=times_struct.down(ISI_bins);
            
        elseif bins2peak(i,2)<5 && mod(bins2peak(i,2)-1,2)==1
            bin_locs_overlap=[l_bound(i+1):r_bound(i)];
            isolate_bins=bin_locs_overlap(1:length(bin_locs_overlap));
            half_ISI_bin=isolate_bins(ceil(length(isolate_bins)/2));
            
            if cdt(i)<0
                edge_bin=times_struct.edges(half_ISI_bin);
            elseif cdt(i)>0
                edge_bin=times_struct.edges(half_ISI_bin-1);
            end

            half_ISI=bincount_ISI(half_ISI_bin);
            
            %edge_bin=times_struct.ISI(half_ISI_bin);
            %half_ISI=(edge_bin+times_struct.ISI(half_ISI_bin+1))/2;
            ISI_bin_split=(and(times_struct.ISI>=edge_bin,times_struct.ISI<half_ISI));
            isolate_spike_count=sum(ISI_bin_split)/duration;
            max_area(i)=max_area(i)+sum(bincountvec(l_bound(i):bin_locs_overlap(1)-1))+isolate_spike_count;
            total_spikes=total_spikes+max_area(i);
            bar_graph.FaceColor='flat';
            for j=l_bound(i):bin_locs_overlap(1)
                bar_graph.CData(j,:)=[1 0 0];
            end
            %flag spikes
            bins=[l_bound(i):bin_locs_overlap(1)-1,half_ISI_bin];
            edge_bins=times_struct.edges(bins);
            ISI_bins=(and(times_struct.ISI>=edge_bins(1),times_struct.ISI<half_ISI));
            up_times{i}=times_struct.up(ISI_bins);
            down_times{i}=times_struct.down(ISI_bins);
            
        else
            max_area(i)=areas(i);
            bar_graph.FaceColor='flat';
            for j=l_bound(i):r_bound(i)
                bar_graph.CData(j,:)=[1 0 0];
            end
            
            %flag spikes
            bins=[l_bound(i):r_bound(i)];
            %edge_bins=times_struct.edges(bins);
            
            if cdt(i)<0
                edge_bins=times_struct.edges([bins,bins(end)+1]);
            elseif cdt(i)>0
                edge_bins=times_struct.edges([bins(1)-1,bins]);
            else
                error('Conduction time not supported.')
            end
            
            ISI_bins=(and(times_struct.ISI>=edge_bins(1),times_struct.ISI<=edge_bins(end)));
            up_times{i}=times_struct.up(ISI_bins);
            down_times{i}=times_struct.down(ISI_bins);
        end
        conduction_time(i)=cdt(i);
        if cdt(i)>0
            direction(i)=1;
        else
            direction(i)=0;
        end
        
    %Right most peak considered, left handed overlap    
    elseif i==length(pks) && length(pks)>1
        if bins2peak(i,1)<5 && mod(bins2peak(i,1)-1,2)==0
            bin_locs_overlap=[l_bound(i):r_bound(i-1)];
            isolate_bins=bin_locs_overlap((length(bin_locs_overlap)/2)+1:end);
            isolate_spike_count=sum(bincountvec(isolate_bins));
            max_area(i)=max_area(i)+sum(bincountvec(bin_locs_overlap(1)+1:r_bound(i)))+isolate_spike_count;
            total_spikes=total_spikes+max_area(i);
            bar_graph.FaceColor='flat';
            for j=bin_locs_overlap(end):r_bound(i)
                bar_graph.CData(j,:)=[1 0 0];
            end
            %flag spikes
            bins=[isolate_bins,bin_locs_overlap(end)+1:r_bound(i)];
            
            if cdt(i)>0
                edge_bins=times_struct.edges([bins(1)-1,bins]);
            elseif cdt(i)<0
                edge_bins=times_struct.edges([bins,bins(end)+1]);
            else
                error('Conduction time not supported')
            end
            
            ISI_bins=(and(times_struct.ISI>=edge_bins(1),times_struct.ISI<edge_bins(end)));
            up_times{i}=times_struct.up(ISI_bins);
            down_times{i}=times_struct.down(ISI_bins);
        elseif bins2peak(i,1)<5 && mod(bins2peak(i,1)-1,2)==1
            bin_locs_overlap=[l_bound(i):r_bound(i-1)];
            isolate_bins=bin_locs_overlap(ceil(length(bin_locs_overlap)/2):end);
            half_ISI_bin=isolate_bins(ceil(length(isolate_bins)/2));
            %edge_bin=times_struct.ISI(half_ISI_bin);
            
            if cdt(i)<0
                edge_bin=times_struct.edges(half_ISI_bin);
            elseif cdt(i)>0
                edge_bin=times_struct.edges(half_ISI_bin-1);
            end

            half_ISI=bincount_ISI(half_ISI_bin);
            
            %half_ISI=(edge_bin+times_struct.ISI(half_ISI_bin-1))/2;
            ISI_bin_split=(and(times_struct.ISI<=edge_bin,times_struct.ISI>=half_ISI));
            isolate_spike_count=sum(ISI_bin_split)/duration;
            max_area(i)=max_area(i)+sum(bincountvec(bin_locs_overlap(1)+1:r_bound(i)))+isolate_spike_count;
            total_spikes=total_spikes+max_area(i);
            bar_graph.FaceColor='flat';
            for j=bin_locs_overlap(end):r_bound(i)
                bar_graph.CData(j,:)=[1 0 0];
            end
            %flag spikes
            bins=[half_ISI_bin,bin_locs_overlap(end)+1:r_bound(i)];
            edge_bins=times_struct.edges(bins);
            ISI_bins=(and(times_struct.ISI<=edge_bins(end),times_struct.ISI>=half_ISI));
            up_times{i}=times_struct.up(ISI_bins);
            down_times{i}=times_struct.down(ISI_bins);
        else
            max_area(i)=areas(i);
            bar_graph.FaceColor='flat';
            for j=l_bound(i):r_bound(i)
                bar_graph.CData(j,:)=[1 0 0];
            end
            
            %flag spikes
            bins=[l_bound(i):r_bound(i)];
            %edge_bins=times_struct.edges(bins);
            
            if cdt(i)<0
                edge_bins=times_struct.edges([bins,bins(end)+1]);
            elseif cdt(i)>0
                edge_bins=times_struct.edges([bins(1)-1,bins]);
            else
                error('Conduction time not supported.')
            end
            
            ISI_bins=(and(times_struct.ISI>=edge_bins(1),times_struct.ISI<=edge_bins(end)));
            up_times{i}=times_struct.up(ISI_bins);
            down_times{i}=times_struct.down(ISI_bins);
        end
        conduction_time(i)=cdt(i);
        if cdt(i)>0
            direction(i)=1;
        else
            direction(i)=0;
        end
        
    % Both sides considered    
    elseif ~isempty(bins2peak)
        if bins2peak(i,2)<5 && bins2peak(i,1)<5 
            %overlap on right
            if bins2peak(i,2)<5 && mod(bins2peak(i,2)-1,2)==0
                bin_locs_overlap=[l_bound(i+1):r_bound(i)];
                isolate_bins=bin_locs_overlap(1:length(bin_locs_overlap)/2);
                isolate_spike_count=sum(bincountvec(isolate_bins));
                max_area(i)=max_area(i)+sum(bincountvec(l_bound(i):bin_locs_overlap(1)-1))+isolate_spike_count;
                total_spikes=total_spikes+max_area(i);
                bar_graph.FaceColor='flat';
                for j=l_bound(i):bin_locs_overlap(1)
                    bar_graph.CData(j,:)=[1 0 0];
                end
                %flag spikes
                bins=[l_bound(i):bin_locs_overlap(1)-1,isolate_bins];
                %edge_bins=times_struct.edges(bins);
                
                if cdt(i)<0
                    edge_bins=times_struct.edges([bins,bins(end)+1]);
                elseif cdt(i)>0
                    edge_bins=times_struct.edges([bins(1)-1,bins]);
                else
                    error('Conduction time not supported.')
                end
                
                ISI_bins=(and(times_struct.ISI>=edge_bins(1),times_struct.ISI<edge_bins(end)));
                up_times{i}=times_struct.up(ISI_bins);
                down_times{i}=times_struct.down(ISI_bins);
                rightmost=isolate_bins(end);
                
            elseif bins2peak(i,2)<5 && mod(bins2peak(i,2)-1,2)==1
                bin_locs_overlap=[l_bound(i+1):r_bound(i)];
                isolate_bins=bin_locs_overlap(1:length(bin_locs_overlap));
                half_ISI_bin=isolate_bins(ceil(length(isolate_bins)/2));
                %edge_bin=times_struct.ISI(half_ISI_bin);
                %half_ISI=(edge_bin+times_struct.ISI(half_ISI_bin+1))/2;
                
                if cdt(i)<0
                    edge_bin=times_struct.edges(half_ISI_bin);
                elseif cdt(i)>0
                    edge_bin=times_struct.edges(half_ISI_bin-1);
                end
                
                half_ISI=bincount_ISI(half_ISI_bin);
                
                ISI_bin_split=(and(times_struct.ISI>=edge_bin,times_struct.ISI<half_ISI));
                isolate_spike_count=sum(ISI_bin_split)/duration;
                max_area(i)=max_area(i)+sum(bincountvec(l_bound(i):bin_locs_overlap(1)-1))+isolate_spike_count;
                total_spikes=total_spikes+max_area(i);
                bar_graph.FaceColor='flat';
                for j=l_bound(i):bin_locs_overlap(1)
                    bar_graph.CData(j,:)=[1 0 0];
                end
                %flag spikes
                bins=[l_bound(i):bin_locs_overlap(1)-1,half_ISI_bin];
                edge_bins=times_struct.edges(bins);
                ISI_bins=(and(times_struct.ISI>=edge_bins(1),times_struct.ISI<half_ISI));
                up_times{i}=times_struct.up(ISI_bins);
                down_times{i}=times_struct.down(ISI_bins);
                rightmost=half_ISI_bin(end);
            end
            %overlap on left
            if bins2peak(i,1)<5 && mod(bins2peak(i,1)-1,2)==0
                bin_locs_overlap=[l_bound(i):r_bound(i-1)];
                isolate_bins=bin_locs_overlap((length(bin_locs_overlap)/2)+1:end);
                isolate_spike_count=sum(bincountvec(isolate_bins));
                max_area(i)=max_area(i)+isolate_spike_count;
                total_spikes=total_spikes+max_area(i);
                bar_graph.FaceColor='flat';
                for j=bin_locs_overlap(end):r_bound(i)
                    bar_graph.CData(j,:)=[1 0 0];
                end
                %flag spikes
                bins=[isolate_bins,bin_locs_overlap(end)+1:r_bound(i)];
                %edge_bins=times_struct.edges([bins(1)-1,bins]);
                
                if cdt(i)>0
                    edge_bins=times_struct.edges([bins(1)-1,bins]);
                elseif cdt(i)<0
                    edge_bins=times_struct.edges([bins,bins(end)+1]);
                else
                    error('Conduction time not supported')
                end
                
                ISI_bins=(and(times_struct.ISI>=edge_bins(1),times_struct.ISI<edge_bins(end)));
                up_times{i}=[up_times{i},times_struct.up(ISI_bins)];
                down_times{i}=[down_times{i},times_struct.down(ISI_bins)];
                leftmost=isolate_bins(1);
                
            elseif bins2peak(i,1)<5 && mod(bins2peak(i,1)-1,2)==1
                bin_locs_overlap=[l_bound(i-1):r_bound(i)];
                %isolate_bins=bin_locs_overlap(1:length(bin_locs_overlap));
                isolate_bins=bin_locs_overlap(ceil(length(bin_locs_overlap)/2):end);
                half_ISI_bin=isolate_bins(ceil(length(isolate_bins)/2));
                %edge_bin=times_struct.ISI(half_ISI_bin);
                %half_ISI=(edge_bin+times_struct.ISI(half_ISI_bin-1))/2;
                
                if cdt(i)<0
                    edge_bin=times_struct.edges(half_ISI_bin);
                elseif cdt(i)>0
                    edge_bin=times_struct.edges(half_ISI_bin-1);
                end
                
                half_ISI=bincount_ISI(half_ISI_bin);
                
                ISI_bin_split=(and(times_struct.ISI<=edge_bin,times_struct.ISI>=half_ISI));
                isolate_spike_count=sum(ISI_bin_split)/duration;
                max_area(i)=max_area(i)+isolate_spike_count;
                total_spikes=total_spikes+max_area(i);
                bar_graph.FaceColor='flat';
                for j=bin_locs_overlap(end):r_bound(i)
                    bar_graph.CData(j,:)=[1 0 0];
                end
                %flag spikes
                bins=[half_ISI_bin,bin_locs_overlap(end)+1:r_bound(i)];
                edge_bins=times_struct.edges(bins);
                ISI_bins=(and(times_struct.ISI<=edge_bins(end),times_struct.ISI>=half_ISI));
                up_times{i}=[up_times{i}, times_struct.up(ISI_bins)];
                down_times{i}=[down_times{i},times_struct.down(ISI_bins)];
                leftmost=half_ISI_bin(1);
            end
            
            up_times{i}=unique(up_times{i});
            down_times{i}=unique(down_times{i});
            
            
            for j=leftmost:rightmost
                bar_graph.CData(j,:)=[1 0 0];
            end
            
        %left overlap
        elseif bins2peak(i,1)<5 && ~bins2peak(i,2)<5
            if bins2peak(i,1)<5 && mod(bins2peak(i,1)-1,2)==0
                bin_locs_overlap=[l_bound(i):r_bound(i-1)];
                isolate_bins=bin_locs_overlap((length(bin_locs_overlap)/2)+1:end);
                isolate_spike_count=sum(bincountvec(isolate_bins));
                max_area(i)=max_area(i)+sum(bincountvec(bin_locs_overlap(1)+1:r_bound(i)))+isolate_spike_count;
                total_spikes=total_spikes+max_area(i);
                bar_graph.FaceColor='flat';
                for j=bin_locs_overlap(end):r_bound(i)
                    bar_graph.CData(j,:)=[1 0 0];
                end
                %flag spikes
                bins=[isolate_bins,bin_locs_overlap(end)+1:r_bound(i)];
                %edge_bins=times_struct.edges([bins(1)-1,bins]);
                
                if cdt(i)>0
                    edge_bins=times_struct.edges([bins(1)-1,bins]);
                elseif cdt(i)<0
                    edge_bins=times_struct.edges([bins,bins(end)+1]);
                else
                    error('Conduction time not supported')
                end
                
                ISI_bins=(and(times_struct.ISI>=edge_bins(1),times_struct.ISI<edge_bins(end)));
                up_times{i}=times_struct.up(ISI_bins);
                down_times{i}=times_struct.down(ISI_bins);
            elseif bins2peak(i,1)<5 && mod(bins2peak(i,1)-1,2)==1
                bin_locs_overlap=[l_bound(i):r_bound(i-1)];
                %isolate_bins=bin_locs_overlap(1:length(bin_locs_overlap));
                isolate_bins=bin_locs_overlap(ceil(length(bin_locs_overlap)/2):end);
                half_ISI_bin=isolate_bins(ceil(length(isolate_bins)/2));
                %edge_bin=times_struct.ISI(half_ISI_bin);
                %half_ISI=(edge_bin+times_struct.ISI(half_ISI_bin-1))/2;
                
                if cdt(i)<0
                    edge_bin=times_struct.edges(half_ISI_bin);
                elseif cdt(i)>0
                    edge_bin=times_struct.edges(half_ISI_bin-1);
                end
                
                half_ISI=bincount_ISI(half_ISI_bin);
                
                ISI_bin_split=(and(times_struct.ISI<=edge_bin,times_struct.ISI>=half_ISI));
                isolate_spike_count=sum(ISI_bin_split)/duration;
                max_area(i)=max_area(i)+sum(bincountvec(bin_locs_overlap(1)+1:r_bound(i)))+isolate_spike_count;
                total_spikes=total_spikes+max_area(i);
                bar_graph.FaceColor='flat';
                for j=bin_locs_overlap(end):r_bound(i)
                    bar_graph.CData(j,:)=[1 0 0];
                end
                %flag spikes
                bins=[half_ISI_bin,bin_locs_overlap(end)+1:r_bound(i)];
                edge_bins=times_struct.edges(bins);
                ISI_bins=(and(times_struct.ISI<=edge_bins(end),times_struct.ISI>=half_ISI));
                up_times{i}=times_struct.up(ISI_bins);
                down_times{i}=times_struct.down(ISI_bins);
            else
                max_area(i)=areas(i);
                bar_graph.FaceColor='flat';
                for j=l_bound(i):r_bound(i)
                    bar_graph.CData(j,:)=[1 0 0];
                end

                %flag spikes
                bins=[l_bound(i):r_bound(i)];
                %edge_bins=times_struct.edges(bins);
                
                if cdt(i)<0
                    edge_bins=times_struct.edges([bins,bins(end)+1]);
                elseif cdt(i)>0
                    edge_bins=times_struct.edges([bins(1)-1,bins]);
                else
                    error('Conduction time not supported.')
                end
                
                ISI_bins=(and(times_struct.ISI>=edge_bins(1),times_struct.ISI<=edge_bins(end)));
                up_times{i}=times_struct.up(ISI_bins);
                down_times{i}=times_struct.down(ISI_bins);
            end
            conduction_time(i)=cdt(i);
            if cdt(i)>0
                direction(i)=1;
            else
                direction(i)=0;
            end
        %right overlap
        elseif bins2peak(i,2)<5 && ~bins2peak(i,1)<5
            if bins2peak(i,2)<5 && mod(bins2peak(i,2)-1,2)==0
                bin_locs_overlap=[l_bound(i+1):r_bound(i)];
                isolate_bins=bin_locs_overlap(1:length(bin_locs_overlap)/2);
                isolate_spike_count=sum(bincountvec(isolate_bins));
                max_area(i)=max_area(i)+sum(bincountvec(l_bound(i):bin_locs_overlap(1)-1))+isolate_spike_count;
                total_spikes=total_spikes+max_area(i);
                bar_graph.FaceColor='flat';
                for j=l_bound(i):bin_locs_overlap(1)
                    bar_graph.CData(j,:)=[1 0 0];
                end
                %flag spikes
                bins=[l_bound(i):bin_locs_overlap(1)-1,isolate_bins];
                %edge_bins=times_struct.edges(bins);
                
                if cdt(i)<0
                    edge_bins=times_struct.edges([bins,bins(end)+1]);
                elseif cdt(i)>0
                    edge_bins=times_struct.edges([bins(1)-1,bins]);
                else
                    error('Conduction time not supported.')
                end
                
                ISI_bins=(and(times_struct.ISI>=edge_bins(1),times_struct.ISI<edge_bins(end)));
                up_times{i}=times_struct.up(ISI_bins);
                down_times{i}=times_struct.down(ISI_bins);

            elseif bins2peak(i,2)<5 && mod(bins2peak(i,2)-1,2)==1
                bin_locs_overlap=[l_bound(i+1):r_bound(i)];
                isolate_bins=bin_locs_overlap(1:length(bin_locs_overlap));
                half_ISI_bin=isolate_bins(ceil(length(isolate_bins)/2));
                %edge_bin=times_struct.ISI(half_ISI_bin);
                
                if cdt(i)<0
                    edge_bin=times_struct.edges(half_ISI_bin);
                elseif cdt(i)>0
                    edge_bin=times_struct.edges(half_ISI_bin-1);
                end
                
                half_ISI=bincount_ISI(half_ISI_bin);
                ISI_bin_split=(and(times_struct.ISI>=edge_bin,times_struct.ISI<half_ISI));
                isolate_spike_count=sum(ISI_bin_split)/duration;
                max_area(i)=max_area(i)+sum(bincountvec(l_bound(i):bin_locs_overlap(1)-1))+isolate_spike_count;
                total_spikes=total_spikes+max_area(i);
                bar_graph.FaceColor='flat';
                for j=l_bound(i):bin_locs_overlap(1)
                    bar_graph.CData(j,:)=[1 0 0];
                end
                %flag spikes
                bins=[l_bound(i):bin_locs_overlap(1)-1,half_ISI_bin];
                %edge_bins=times_struct.edges([bins(1)-1,bins]);
                edge_bins=times_struct.edges(bins);
                ISI_bins=(and(times_struct.ISI>=edge_bins(1),times_struct.ISI<half_ISI));
                up_times{i}=times_struct.up(ISI_bins);
                down_times{i}=times_struct.down(ISI_bins);

            else
                max_area(i)=areas(i);
                bar_graph.FaceColor='flat';
                for j=l_bound(i):r_bound(i)
                    bar_graph.CData(j,:)=[1 0 0];
                end

                %flag spikes
                bins=[l_bound(i):r_bound(i)];
                %edge_bins=times_struct.edges(bins);
                
                if cdt(i)<0
                    edge_bins=times_struct.edges([bins,bins(end)+1]);
                elseif cdt(i)>0
                    edge_bins=times_struct.edges([bins(1)-1,bins]);
                else
                    error('Conduction time not supported.')
                end
                
                ISI_bins=(and(times_struct.ISI>=edge_bins(1),times_struct.ISI<=edge_bins(end)));
                up_times{i}=times_struct.up(ISI_bins);
                down_times{i}=times_struct.down(ISI_bins);
            end
            conduction_time(i)=cdt(i);
            if cdt(i)>0
                direction(i)=1;
            else
                direction(i)=0;
            end
        else
            max_area(i)=areas(i);
            bar_graph.FaceColor='flat';
            for j=l_bound(i):r_bound(i)
                bar_graph.CData(j,:)=[1 0 0];
            end
            
            %flag spikes
            bins=[l_bound(i):r_bound(i)];
            %edge_bins=times_struct.edges(bins);
            
            if cdt(i)<0
                edge_bins=times_struct.edges([bins,bins(end)+1]);
            elseif cdt(i)>0
                edge_bins=times_struct.edges([bins(1)-1,bins]);
            else
                error('Conduction time not supported.')
            end
            
            ISI_bins=(and(times_struct.ISI>=edge_bins(1),times_struct.ISI<=edge_bins(end)));
            up_times{i}=times_struct.up(ISI_bins);
            down_times{i}=times_struct.down(ISI_bins);
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