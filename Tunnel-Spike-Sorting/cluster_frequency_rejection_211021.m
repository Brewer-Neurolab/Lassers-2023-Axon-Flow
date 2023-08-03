function cluster_frequency_rejection_211021(folder_address,rec_leng)

cap_freq=10;
max_clus_size=cap_freq*rec_leng;

cd(folder_address)
times_dir=dir('times_*.mat');
times_name=string({times_dir.name});

for i=1:length(times_name)
    load(times_name(i))
    num_classes=max(cluster_class(:,1));
    for j=1:num_classes
        num_spikes=sum(cluster_class(:,1)==j);
        if num_spikes>max_clus_size
            cluster_class(cluster_class(:,1)==j)=0;
            save(times_name(i),'spikes','cluster_class','par','gui_status','temp','forced','inspk','ipermut')
        end
    end
end

end