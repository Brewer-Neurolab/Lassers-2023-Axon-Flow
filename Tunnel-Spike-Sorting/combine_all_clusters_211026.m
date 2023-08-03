function combine_all_clusters_211026(folder_address)
cd(folder_address)
times_dir=dir('times_*.mat');
times_name=string({times_dir.name});

for j=1:length(times_name)
    file=load(times_name(j));
    file.cluster_class(file.cluster_class(:,1)~=0,1)=1;
    save(times_name(j),'-struct','file')
end

end