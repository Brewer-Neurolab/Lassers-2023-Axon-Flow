function hi_lo_wave_clus_batch_211020(folder_address,electrodes,rec_leng)

%low wave clus

disp('Getting Low Spikes...')
for i=1:length(folder_address)
    disp(strcat('Getting Low Spikes ',string(i)))
    cd(folder_address{i})
    mkdir lo
    Get_spikes('Files.txt')
end

disp('Clustering Low Spikes...')
for i=1:length(folder_address)
    disp(strcat('Clustering Low Spikes ',string(i)))
    cd(folder_address{i})
    Do_clustering('Files_spikes.txt')
end

%rename 
for i=1:length(folder_address)
    cd(folder_address{i})
    spikes_dir=dir('*_spikes.mat');
    spikes_name=string({spikes_dir.name});
    for j=1:length(spikes_name)
        split_name=strsplit(spikes_name(j),{'.'});
        movefile(spikes_name(j),strcat(split_name(1),'_lo','.mat'))
        movefile(strcat(split_name(1),'_lo','.mat'),'.\lo')
    end
    times_dir=dir('times_*.mat');
    times_name=string({times_dir.name});
    for j=1:length(times_name)
        split_name=strsplit(times_name(j),{'.'});
        movefile(times_name(j),strcat(split_name(1),'_lo','.mat'))
        movefile(strcat(split_name(1),'_lo','.mat'),'.\lo')
    end
    figs=dir('fig2print*.png');
    figs_name=string({figs.name});
    for j=1:length(figs_name)
        split_name=strsplit(figs_name(j),{'.'});
        movefile(figs_name(j),strcat(split_name(1),'_lo','.png'))
        movefile(strcat(split_name(1),'_lo','.png'),'.\lo')
    end
end

%high wave clus

%set new parameters for high clustering
param.stdmin=50.1;
param.stdmax=500;

disp('Getting High Spikes...')
for i=1:length(folder_address)
    disp(strcat('Getting High Spikes ',string(i)))
    cd(folder_address{i})
    mkdir hi
    Get_spikes('Files.txt','par',param)
end

disp('Clustering High Spikes...')
for i=1:length(folder_address)
    disp(strcat('Clustering High Spikes ',string(i)))
    cd(folder_address{i})
    Do_clustering('Files_spikes.txt','par',param)
end

%rename 
for i=1:length(folder_address)
    cd(folder_address{i})
    spikes_dir=dir('*_spikes.mat');
    spikes_name=string({spikes_dir.name});
    for j=1:length(spikes_name)
        split_name=strsplit(spikes_name(j),{'.'});
        movefile(spikes_name(j),strcat(split_name(1),'_hi','.mat'))
        movefile(strcat(split_name(1),'_hi','.mat'),'.\hi')
    end
    times_dir=dir('times_*.mat');
    times_name=string({times_dir.name});
    for j=1:length(times_name)
        split_name=strsplit(times_name(j),{'.'});
        movefile(times_name(j),strcat(split_name(1),'_hi','.mat'))
        movefile(strcat(split_name(1),'_hi','.mat'),'.\hi')
    end
    figs=dir('fig2print*.png');
    figs_name=string({figs.name});
    for j=1:length(figs_name)
        split_name=strsplit(figs_name(j),{'.'});
        movefile(figs_name(j),strcat(split_name(1),'_hi','.png'))
        movefile(strcat(split_name(1),'_hi','.png'),'.\hi')
    end
end

% assemble matching high and low times files and merge cluster classes
for i=1:length(folder_address)
    cd(folder_address{i})
    lo_dir=dir('.\lo\times_*');
    lo_names=string({lo_dir.name});
    hi_dir=dir('.\hi\times_*');
    hi_names=string({hi_dir.name});
    for j=1:length(electrodes)
        lo_electrode=lo_names(contains(lo_names,strcat('_',electrodes{j},'_')));
        hi_electrode=hi_names(contains(hi_names,strcat('_',electrodes{j},'_')));
        if ~isempty(lo_electrode) && ~isempty(hi_electrode)
            load_lo_str=strcat('.\lo\',lo_electrode);
            load_hi_str=strcat('.\hi\',hi_electrode);
            combine_hi_lo_wave_clus_211020(load_lo_str,load_hi_str,lo_electrode)
        elseif ~isempty(lo_electrode) && isempty(hi_electrode)
            vars={'spikes','cluster_class','par','gui_status','Temp','forced','inspk','ipermut'};
            clear(vars{:})
            load_lo_str=strcat('.\lo\',lo_electrode);
            load(load_lo_str)
            split_name=string(strsplit(lo_electrode,{'_','.'}));
            name=strcat(split_name(1),'_',split_name(2),'.mat');
            save(name,'spikes','cluster_class','par','gui_status','Temp','forced','inspk','ipermut')
        elseif isempty(lo_electrode) && ~isempty(hi_electrode)
            vars={'spikes','cluster_class','par','gui_status','Temp','forced','inspk','ipermut'};
            clear(vars{:})
            load_hi_str=strcat('.\hi\',hi_electrode);
            load(load_hi_str)
            split_name=string(strsplit(hi_electrode,{'_','.'}));
            name=strcat(split_name(1),'_',split_name(2),'.mat');
            save(name,'spikes','cluster_class','par','gui_status','Temp','forced','inspk','ipermut')
        end
    end
end

% 0 clusters over 3000 spikes, too much for NMI
% for i=1:length(folder_address)
%     cluster_frequency_rejection_211021(folder_address{i},rec_leng)
% end

end