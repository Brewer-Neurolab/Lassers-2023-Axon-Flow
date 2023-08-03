function combine_hi_lo_wave_clus_211020(lo_times,hi_times,lo_name)
%inputs should be directories represented as char or string
%when two thresholds are set for the same channel, this concatonates the
%times files produced by wave clus

lo=load(lo_times);
hi=load(hi_times);

% renumber clusters
largest_lo_cluster=max(lo.cluster_class(:,1));
%ignore 0 cluster, keep 0 as 0
hi.cluster_class(hi.cluster_class(:,1)~=0,1)=...
    hi.cluster_class(hi.cluster_class(:,1)~=0,1)+largest_lo_cluster;

%concatonate all the variables of each times file into one file 
spikes=[lo.spikes;hi.spikes];
cluster_class=[lo.cluster_class;hi.cluster_class];
par=[lo.par;hi.par];
gui_status=[lo.gui_status;hi.gui_status];
Temp=[lo.Temp,hi.Temp];
forced=[lo.forced,hi.forced];
inspk={lo.inspk;hi.inspk};
ipermut={lo.ipermut,hi.ipermut};

split_name=string(strsplit(lo_name,{'_','.'}));
name=strcat(split_name(1),'_',split_name(2),'.mat');

save(name,'spikes','cluster_class','par','gui_status','Temp','forced','inspk','ipermut')
end