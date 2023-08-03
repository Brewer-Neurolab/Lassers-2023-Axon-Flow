function assign_fid_211118(folder_address,parent_dir,is_prestim)

cd(folder_address)
current_dir=dir;
dir_names=string({current_dir.name});
folders=isfolder(dir_names);
folders([1,2])=0;
dir_names=dir_names(folders);

for i=1:length(dir_names)
    fid=i;
    filename=dir_names(i);
    if is_prestim
        save(strcat(parent_dir,'\FID',string(i),'_Spon.mat'),'fid','filename')
    elseif ~is_prestim
        save(strcat(parent_dir,'\FID',string(i),'_Stim.mat'),'fid','filename')
    else
        error('Define if this is prestimulation or poststimulation.')
    end
    
end

end