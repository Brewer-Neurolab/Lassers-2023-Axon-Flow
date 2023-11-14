function write_new_tunnel_pseudo_times(pseudo_sub_dir, spike_idx, time_idx,is_ff, channels_to_write)
% you can restore original pseudo_times from
% create_pseudo_times_files_220216.m

spike_times=time_idx(spike_idx)*1e3;
clus_num=ones(length(spike_times),1);
cluster_class=[clus_num,spike_times'];

for i=1:length(channels_to_write)
    if is_ff
        save((pseudo_sub_dir+'\times_'+string(channels_to_write(i))+'_ff.mat'),"cluster_class")
    else
        save((pseudo_sub_dir+'\times_'+string(channels_to_write(i))+'_fb.mat'),"cluster_class")
    end
end

end