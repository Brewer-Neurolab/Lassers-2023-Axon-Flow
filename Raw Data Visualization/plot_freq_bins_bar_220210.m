function [boi]=plot_freq_bins_bar_220210(time_idx,index,chan_2_plot,bin_width_sec,fs)

bins_edges=[time_idx(1):time_idx(bin_width_sec*fs):time_idx(end),time_idx(end)];
index=index/1e3;
vec_2_plot=[];
for pts=1:length(bins_edges)-1
    freq=sum(index>bins_edges(pts)&index<=bins_edges(pts+1))/bin_width_sec;
    vec_2_plot(pts)=freq;
end

figure('units','normalized','outerposition',[0 0 1 1])
hold on
centers=convert_edges_2_centers(bins_edges);
bar(centers,vec_2_plot)
title(strcat(chan_2_plot,'-',string(bin_width_sec),'s bins'))
xlabel("Seconds")
ylabel("Frequency (HZ)")

avg_spike_rate=length(index)/time_idx(end);
binned_avg_spike_rate=mean(vec_2_plot(vec_2_plot~=0));

if ~isnan(avg_spike_rate)
    yline(avg_spike_rate,'g')
end
if ~isnan(binned_avg_spike_rate)
    yline(binned_avg_spike_rate,'m')
end

%bins of intrest
range_avg=[binned_avg_spike_rate-(0.2*binned_avg_spike_rate),binned_avg_spike_rate+(0.2*binned_avg_spike_rate)];
boi=centers(vec_2_plot>=range_avg(1)&vec_2_plot<=range_avg(2));

hold off
end