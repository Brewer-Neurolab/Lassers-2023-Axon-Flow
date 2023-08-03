function wave_clus_GUI_quick_open_211004(folder_address)
% This function opens the times files in the Wave Clus GUI for manual
% fixing of clusters. Upon closing the Wave Clus GUI this function will
% open a the next times files specified in folder_address. The input
% folder_address should be a string or char vector pertaining to a
% directory where times files are located.

cd(folder_address)
close all
times_folder=dir();
times_files=times_folder(contains(string({times_folder.name}),'times'));
times_files=string({times_files.name});

for i=1:length(times_files)
    wave_clus(char(times_files(i)))
    wave_gui=gcf;
    uiwait(wave_gui)
end

end