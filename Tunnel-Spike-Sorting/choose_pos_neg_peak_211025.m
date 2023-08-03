function choose_pos_neg_peak_211025(folder_address,electrodes)
%looks at tihe positive and negative peaks from detection and determines if
%there are spikes within 1.5 ms dead time, tallest spike is kept others are
%discarded

%% Positive spikes
disp('Getting Low Spikes (Pos)...')

%set new parameters for low clustering
param.stdmin=11;
param.stdmax=50;
for i=1:length(folder_address)
    disp(strcat('Getting Low Spikes ',string(i)))
    cd(folder_address{i})
    mkdir lo
    param.detection='pos';
    Get_spikes('Files.txt','par',param)
end

%rename 
for i=1:length(folder_address)
    cd(folder_address{i})
    spikes_dir=dir('*_spikes.mat');
    spikes_name=string({spikes_dir.name});
    for j=1:length(spikes_name)
        split_name=strsplit(spikes_name(j),{'.'});
        movefile(spikes_name(j),strcat(split_name(1),'_lo_pos','.mat'))
        movefile(strcat(split_name(1),'_lo_pos','.mat'),'.\lo')
    end
end

%high wave clus

%set new parameters for high clustering
param.stdmin=50.1;
param.stdmax=500;

disp('Getting High Spikes (Pos)...')
for i=1:length(folder_address)
    disp(strcat('Getting High Spikes ',string(i)))
    cd(folder_address{i})
    mkdir hi
    param.detection='pos';
    Get_spikes('Files.txt','par',param)
end

%rename 
for i=1:length(folder_address)
    cd(folder_address{i})
    spikes_dir=dir('*_spikes.mat');
    spikes_name=string({spikes_dir.name});
    for j=1:length(spikes_name)
        split_name=strsplit(spikes_name(j),{'.'});
        movefile(spikes_name(j),strcat(split_name(1),'_hi_pos','.mat'))
        movefile(strcat(split_name(1),'_hi_pos','.mat'),'.\hi')
    end
end

%% Negative
disp('Getting Low Spikes (Neg)...')

%set new parameters for low clustering
param.stdmin=11;
param.stdmax=50;
for i=1:length(folder_address)
    disp(strcat('Getting Low Spikes ',string(i)))
    cd(folder_address{i})
    mkdir lo
    param.detection='neg';
    Get_spikes('Files.txt','par',param)
end

%rename 
for i=1:length(folder_address)
    cd(folder_address{i})
    spikes_dir=dir('*_spikes.mat');
    spikes_name=string({spikes_dir.name});
    for j=1:length(spikes_name)
        split_name=strsplit(spikes_name(j),{'.'});
        movefile(spikes_name(j),strcat(split_name(1),'_lo_neg','.mat'))
        movefile(strcat(split_name(1),'_lo_neg','.mat'),'.\lo')
    end
end

%high wave clus

%set new parameters for high clustering
param.stdmin=50.1;
param.stdmax=500;

disp('Getting High Spikes (Neg)...')
for i=1:length(folder_address)
    disp(strcat('Getting High Spikes ',string(i)))
    cd(folder_address{i})
    mkdir hi
    param.detection='neg';
    Get_spikes('Files.txt','par',param)
end

%rename 
for i=1:length(folder_address)
    cd(folder_address{i})
    spikes_dir=dir('*_spikes.mat');
    spikes_name=string({spikes_dir.name});
    for j=1:length(spikes_name)
        split_name=strsplit(spikes_name(j),{'.'});
        movefile(spikes_name(j),strcat(split_name(1),'_hi_neg','.mat'))
        movefile(strcat(split_name(1),'_hi_neg','.mat'),'.\hi')
    end
end

%% Compare positive and negative spikes

for i=1:length(folder_address)
    cd(folder_address{i})
    for j=1:length(electrodes)
        %load and manipulate relevant data
%         pos_lo=dir(strcat('\lo\',electrodes{j},'_spikes_lo_pos.mat'));
%         pos_lo=string({pos_lo.name});
        pos_lo=load(strcat('.\lo\',electrodes{j},'_spikes_lo_pos.mat'));
        pos_lo_spikes=pos_lo.spikes;
        pos_lo_index=pos_lo.index;
        pos_lo_mid_idx=pos_lo.par.w_pre;
        pos_lo_mid_val=pos_lo_spikes(:,pos_lo_mid_idx);
%         pos_hi=dir(strcat('\hi\',electrodes{j},'*_spikes_hi_pos'));
%         pos_hi=string({pos_hi.name});
        pos_hi=load('.\hi\',electrodes{j},'_spikes_hi_pos');
        pos_hi_spikes=pos_hi.spikes;
        pos_hi_index=pos_hi.index;
        pos_hi_mid_idx=pos_hi.par.w_pre;
        pos_hi_mid_val=pos_hi_spikes(:,pos_hi_mid_idx);
%         neg_lo=dir(strcat('\lo\',electrodes{j},'*_spikes_lo_neg'));
%         neg_lo=string({neg_lo.name});
        neg_lo=load('.\lo\',electrodes{j},'_spikes_lo_neg');
        neg_lo_spikes=neg_lo.spikes;
        neg_lo_index=neg_lo.index;
        neg_lo_mid_idx=neg_lo.par.w_pre;
        neg_lo_mid_val=neg_lo_spikes(:,neg_lo_mid_idx);
%         neg_hi=dir(strcat('\hi\',electrodes{j},'*_spikes_hi_neg'));
%         neg_hi=string({neg_hi.name});
        neg_hi=load('.\hi\',electrodes{j},'_spikes_hi_neg');
        neg_hi_spikes=neg_hi.spikes;
        neg_hi_index=neg_hi.index;
        neg_hi_mid_idx=neg_hi.par.w_pre;
        neg_hi_mid_val=neg_hi_spikes(:,neg_hi_mid_idx);
        
        %find unique spike idx within tolerance
        spikes_idx=[pos_lo_index,pos_hi_index,neg_lo_index,neg_hi_index];
        spikes_amp=[pos_lo_mid_val,pos_hi_mid_val,neg_lo_mid_val,neg_hi_mid_val];
        spikes_spikes=[pos_lo_spikes;pos_hi_spikes;neg_lo_spikes;neg_hi_spikes];
        dead_time=1.5; %ms
        [~,~,spikes_tol]=uniquetol(spikes_idx,dead_time);
        to_delete=zeros(1,length(spikes_tol));
        for unique_times=1:max(spikes_tol)
            time_idx=spikes_tol==unique_times;
            [~,largest_amp_idx]=max(abs(spikes_amp(time_idx)));
            time_idx_ones=find(time_idx==1);
            time_idx(time_idx_ones(largest_amp_idx))=0;
            to_delete=or(to_delete,time_idx);
        end
        spikes_idx(to_delete)=0;
        [~,c_spikes]=size(pos_lo_spikes);
        spikes_spikes(to_delete,:)=zeros(1,c_spikes);
        pos_lo_index=spikes_idx(1:length(pos_lo_index));
        pos_hi_index=spikes_idx(length(pos_lo_index)+1:length(pos_hi_index));
        current_idx=length(pos_lo_index)+1+length(pos_hi_index);
        neg_lo_index=spikes_idx(current_idx:current_idx+length(neg_lo_index));
        current_idx=current_idx+length(neg_lo_index);
        neg_hi_index=spikes_idx(current_idx:current_idx+length(neg_hi_index));
        
        pos_lo_spikes=spikes_spikes(1:length(pos_lo_spikes(:,1)),:);
        pos_hi_spikes=spikes_spikes(length(pos_lo_spikes(:,1))+1:length(pos_hi_spikes(:,1)),:);
        current_idx=length(pos_lo_spikes(:,1))+1+length(pos_hi_spikes(:,1));
        neg_lo_spikes=spikes_spikes(current_idx:current_idx+length(neg_lo_spikes(:,1)),:);
        current_idx=current_idx+length(neg_lo_spikes(:,1));
        neg_hi_spikes=spikes_spikes(current_idx:current_idx+length(neg_hi_spikes(:,1)),:);
        
        %remove zero rows
        pos_lo_index(pos_lo_index==0)=[];
        pos_hi_index(pos_hi_index==0)=[];
        neg_lo_index(neg_lo_index==0)=[];
        neg_hi_index(neg_hi_index==0)=[];
        
        sum_spikes=sum(pos_lo_spikes,2);
        pos_lo_spikes(sum_spikes==0,:)=[];
        sum_spikes=sum(pos_hi_spikes,2);
        pos_hi_spikes(sum_spikes==0,:)=[];
        sum_spikes=sum(neg_lo_spikes,2);
        neg_lo_spikes(sum_spikes==0,:)=[];
        sum_spikes=sum(neg_hi_spikes,2);
        neg_hi_spikes(sum_spikes==0,:)=[];
        
        %assign new fields
        pos_lo.spikes=pos_lo_spikes;
        pos_lo.index=pos_lo_index;
        pos_hi.spikes=pos_hi_spikes;
        pos_hi.index=pos_hi_index;
        neg_lo.spikes=neg_lo_spikes;
        neg_lo.index=neg_lo_index;
        neg_hi.spikes=neg_hi_spikes;
        neg_hi.index=neg_hi_index;
        
        %save new fields
%         save(strcat('\lo\',electrodes{j},'*_spikes_lo_pos.mat'),'-struct',pos_lo)
%         save(strcat('\hi\',electrodes{j},'*_spikes_hi_pos.mat'),'-struct',pos_hi)
%         save(strcat('\lo\',electrodes{j},'*_spikes_lo_neg.mat'),'-struct',neg_lo)
%         save(strcat('\hi\',electrodes{j},'*_spikes_hi_pos.mat'),'-struct',neg_hi)
    end
end


end