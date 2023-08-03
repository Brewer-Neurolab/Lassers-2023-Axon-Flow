%% Plot Raw Data and Spikes (FF and FB)
% Ensure you are in the correct directory that houses the .mat files for
% raw data and spikes
%visualizes data from spikes table data structure

clear
clc
close all
parent_dir='D:\Brewer lab data\HFS\test 4';
spikes_table=load(strcat(parent_dir,'\','allregion_unit_matched_spon.mat'));
spikes_table=spikes_table.allregion_unit_matched_spon;
%cd('D:\Brewer lab data\HFS\Filtered Spon MCD Files\03-Aug-2021_C')
cd('D:\Brewer lab data\HFS\test 4\Filtered Spon MCD Files\22-Oct-2021_B')
raw_dir=dir;
raw_dir=raw_dir([raw_dir.isdir]);
raw_dir=raw_dir(3:end);

%convert spikes table spike indecies from sample to time indecies

ff_colors=[255 255 0; 255 0 255; 0 255 255; 246 255 0]./255;
fb_colors=[0 0 0; 255 149 0; 128 255 0; 138 77 55]./255;

%%
for i=1:length(spikes_table)
    %cd('D:\Brewer lab data\HFS\Filtered Spon MCD Files\03-Aug-2021_C')
    cd('D:\Brewer lab data\HFS\Filtered Stim Data MCD Files\03-Aug-2021_A')
    cd(strcat('.\',raw_dir(i).name))
    for j=1:length(spikes_table{i}.Subregion)
        current_pair=spikes_table{i}.("Electrode Pairs")(j);
        channels=strsplit(current_pair,{' ','-'});
        for k=1:length(channels)
            data=load(strcat(channels(k),'.mat'));
            data=data.data;
            %data=data./1000000000000;
            data=highpass(data,300,25000);
            x=[1:length(data)];
            x=x./25000;
            figure('Position',[0 0 1920 700])
            hold on
            plot(x,data)
            xlim([0,x(end)])
            ylim([-500,500])
            ylabel("uV")
            xlabel("seconds")
            title(strcat(channels(k), " FID: ",string(i)))
            if k==1
                for up=1:length(spikes_table{i}.up_ff{j})
                    if ~isempty(find(spikes_table{i}.up_ff{j}{up}))
                        xline(find(spikes_table{i}.up_ff{j}{up})./25000,'Color',ff_colors(up,:))
                        x_dots=find(spikes_table{i}.up_ff{j}{up})./25000;
                        plot(x_dots,repmat(500,length(x_dots)),'or')
                    end
                end
                for up=1:length(spikes_table{i}.up_fb{j})
                    if ~isempty(find(spikes_table{i}.up_fb{j}{up}))
                        xline(find(spikes_table{i}.up_fb{j}{up})./25000,'Color',fb_colors(up,:))
                        x_dots=find(spikes_table{i}.up_fb{j}{up})./25000;
                        plot(x_dots,repmat(500,length(x_dots)),'og')
                    end
                end
            elseif k==2
                for up=1:length(spikes_table{i}.down_ff{j})
                    if ~isempty(find(spikes_table{i}.down_ff{j}{up}))
                        xline(find(spikes_table{i}.down_ff{j}{up})./25000,'Color',ff_colors(up,:))
                        x_dots=find(spikes_table{i}.down_ff{j}{up})./25000;
                        plot(x_dots,repmat(500,length(x_dots)),'or')
                    end
                end
                for up=1:length(spikes_table{i}.down_fb{j})
                    if ~isempty(find(spikes_table{i}.down_fb{j}{up}))
                        xline(find(spikes_table{i}.down_fb{j}{up})./25000,'Color',fb_colors(up,:))
                        x_dots=find(spikes_table{i}.down_fb{j}{up})./25000;
                        plot(x_dots,repmat(500,length(x_dots)),'og')
                    end
                end
            else
                error('Three or more tunnels not supported')
            end
            hold off
            saveas(gcf,strcat('D:\Brewer lab data\HFS\Filtered Spon MCD Files\Raw prestim pics','\',channels(k), " FID ",string(i),'.png'))
            saveas(gcf,strcat('D:\Brewer lab data\HFS\Filtered Spon MCD Files\Raw prestim pics','\',channels(k), " FID ",string(i),'.fig'))
            close all
        end
    end
end

%%
%Subplot of paired data
close all
for i=1:length(spikes_table)
    cd('D:\Brewer lab data\HFS\Filtered Spon MCD Files\03-Aug-2021_C')
    cd(strcat('.\',raw_dir(i).name))
    for j=1:length(spikes_table{i}.Subregion)
        current_pair=spikes_table{i}.("Electrode Pairs")(j);
        channels=strsplit(current_pair,{' ','-'});
        %subplot(2,1)
        ax1=subplot(2,1,1);
        ax2=subplot(2,1,2);
        axes=[ax1;ax2];
        for k=1:length(channels)
            subplot(axes(k))
            data=load(strcat(channels(k),'.mat'));
            data=data.data;
            %data=data./1000000000000;
            data=highpass(data,300,25000);
            x=[1:length(data)];
            %x=x./25000;
            %figure('Position',[0 0 1920 700])
            hold on
            plot(x,data)
            xlim([0,x(end)])
            ylim([-500,500])
            ylabel("uV")
            xlabel("samples")
            title(strcat(channels(k), " FID: ",string(i)))
            if k==1
                for up=1:length(spikes_table{i}.up_ff{j})
                    if ~isempty(find(spikes_table{i}.up_ff{j}{up}))
                        xline(find(spikes_table{i}.up_ff{j}{up}),'Color',ff_colors(up,:))
                        x_dots=find(spikes_table{i}.up_ff{j}{up});
                        plot(x_dots,repmat(500,length(x_dots)),'or')
                    end
                end
                for up=1:length(spikes_table{i}.up_fb{j})
                    if ~isempty(find(spikes_table{i}.up_fb{j}{up}))
                        xline(find(spikes_table{i}.up_fb{j}{up}),'Color',fb_colors(up,:))
                        x_dots=find(spikes_table{i}.up_fb{j}{up});
                        plot(x_dots,repmat(500,length(x_dots)),'og')
                    end
                end
            elseif k==2
                for up=1:length(spikes_table{i}.down_ff{j})
                    if ~isempty(find(spikes_table{i}.down_ff{j}{up}))
                        xline(find(spikes_table{i}.down_ff{j}{up}),'Color',ff_colors(up,:))
                        x_dots=find(spikes_table{i}.down_ff{j}{up});
                        plot(x_dots,repmat(500,length(x_dots)),'or')
                    end
                end
                for up=1:length(spikes_table{i}.down_fb{j})
                    if ~isempty(find(spikes_table{i}.down_fb{j}{up}))
                        xline(find(spikes_table{i}.down_fb{j}{up}),'Color',fb_colors(up,:))
                        x_dots=find(spikes_table{i}.down_fb{j}{up});
                        plot(x_dots,repmat(500,length(x_dots)),'og')
                    end
                end
            else
                error('Three or more tunnels not supported')
            end
            hold off
            %saveas(gcf,strcat('D:\Brewer lab data\HFS\Filtered Stim Data MCD Files\Raw Poststim Pics','\',channels(k), " FID ",string(i),'.png'))
            %close all
        end
        saveas(gcf,strcat('D:\Brewer lab data\HFS\Filtered Spon MCD Files\Raw prestim pics','\',channels(1),channels(2), " FID ",string(i),'.png'))
        close all
    end
end