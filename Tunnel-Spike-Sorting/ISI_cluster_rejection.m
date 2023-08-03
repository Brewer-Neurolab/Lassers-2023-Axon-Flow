%% Waveclus Cluster rejection based on ISI

fast_thresh=0.1;

times_files=dir;
times_files=times_files(contains(string({times_files.name}),'times'));
times_files=string({times_files.name});

%compute ISI
num_clus=length(unique(cluster_class(:,1)));

%calculate ISI
for i=1:num_clus
    pts=cluster_class(cluster_class(:,1)==i,2)./1e3;
    ISI{i}=diff(pts);
end

%find clusters with #isi<3ms more than 10% of spikes
clus2del=[];
for i=1:length(ISI)
    num_fast=sum(ISI{i}<0.003);
    perc_fast=num_fast/length(ISI{i});
    if_pass=perc_fast<fast_thresh;
    if ~if_pass && ~isnan(perc_fast)
        clus2del=[clus2del,i];
    end
end

idx_to_delete=logical(sum(cluster_class(:,1)==clus2del,2));
cluster_class(idx_to_delete,:)=[];
spikes(idx_to_delete,:)=[];
inspk(idx_to_delete,:)=[];
ipermut(idx_to_delete)=[];
forced(idx_to_delete)=[];

