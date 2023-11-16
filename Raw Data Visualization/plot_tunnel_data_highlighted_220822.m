%% Plot burst-highlighted raw data for axons in tunnels

clear
clc

% stim_dir="D:\Brewer lab data\HFS\No Stim\23-Nov-2021_B\"; % no stim
% stim_dir="D:\Brewer lab data\HFS\Theta Stim\10-May-2022_A\"; % 5hfs
% stim_dir="D:\Brewer lab data\HFS\HFS Stim\24-Nov-2021_A\"; %40 HFS

% testing
% stim_dir="D:\Brewer lab data\HFS\Theta Stim\31-Aug-2022_A\";
stim_dir="D:\Brewer lab data\HFS\Theta Stim\27-Sep-2022_A\";

stim_dir_struct=dir(stim_dir);
stim_files={stim_dir_struct.name};
is_dir=[stim_dir_struct.isdir];
stim_folders=stim_files(is_dir);
stim_folders=stim_folders(3:end)';

%testing
load("D:\Brewer lab data\HFS\Theta Stim\27-Sep-2022_A\allregion_unit_matched_stim.mat")

% load("D:\Brewer lab data\HFS\No Stim\23-Nov-2021_A\allregion_unit_matched_stim.mat") % nostim
% load("D:\Brewer lab data\HFS\Theta Stim\10-May-2022_A\allregion_unit_matched_stim.mat") % 5 HFS
% load('D:\Brewer lab data\HFS\HFS Stim\24-Nov-2021_A\allregion_unit_matched_stim.mat') %40 hfs
[spike_burst_dyn_table_stim_orig]=compute_spike_burst_dynamics_axons_220822(allregion_unit_matched_stim);



blues=["#0000fe","#4d4cff","#3f00ff","#9683ec"];
reds=["#ff0101","#ff0101","#db2d44","#fb607f"];

fs=25000;
t=[1/fs:1/fs:300];
%% Get Subregion

% subregion_OI=1;
subregions=["EC-DG","DG-CA3","CA3-CA1","CA1-EC"];
for subregion_OI=2%1:length(subregions)

    spike_burst_dyn_table_stim=spike_burst_dyn_table_stim_orig(spike_burst_dyn_table_stim_orig.regi==subregion_OI,:);



    for fi=1:9

        axons_fi=spike_burst_dyn_table_stim(spike_burst_dyn_table_stim.fi==fi,:);

        if isempty(axons_fi)
            continue
        end

        clear tunnel_names
        for i=1:height(axons_fi)
            tunnel_names(i,:)=strsplit(axons_fi.channel_name(i),{'-'});
        end
        tunnel_names=unique(tunnel_names(:,1));

        % get electrode pairs, which one matches
        allregion_subregion=allregion_unit_matched_stim{fi}(allregion_unit_matched_stim{fi}.regi==subregion_OI,:);
        elec_pairs=allregion_subregion.("Electrode Pairs");

        for i=1:length(elec_pairs)
            elec_pairs_split(i,:)=strsplit(elec_pairs(i),{'-'});
        end

        %     bursting_bounds=spike_burst_dyn_table_stim.BurstBounds(spike_burst_dyn_table_stim.fi==fi,:);

        for ti=1:length(tunnel_names)
            elecs_idx=any(strcmp(elec_pairs_split,tunnel_names(ti)),2);
            elecs=elec_pairs_split(elecs_idx,:);


            % up tunnel
            figure('units','normalized','outerposition',[0 0 1 1])
            load(stim_dir+stim_folders{fi}+"\"+elecs(1)+".mat")
            %         threshold=load(stim_dir+stim_folders{fi}+"\"+elecs(1)+".mat")
            data=bandpass(data,[300,3000],fs);
            plot(t,data,'k')
            hold on

            y_low=min(data);
            y_low=y_low(1)-10;
            y_hi=max(data);
            y_hi=y_hi(1)+10;
            ylim([y_low,y_hi])
            ylabel("uV")
            xlabel("seconds")

            for ai=1:length(allregion_subregion.up_ff{elecs_idx})

                spikes=allregion_subregion.up_ff{elecs_idx}{ai};
                plot(t(spikes),data(spikes),'.','Color',hex2rgb(char(reds(ai))))

                if any(contains(spike_burst_dyn_table_stim.channel_name(contains(spike_burst_dyn_table_stim.channel_name,tunnel_names(ti))...
                        & spike_burst_dyn_table_stim.fi==fi & spike_burst_dyn_table_stim.if_ff==1,:),tunnel_names(ti)))
                    burst_bounds=spike_burst_dyn_table_stim.BurstBounds(contains(spike_burst_dyn_table_stim.channel_name,tunnel_names(ti))...
                        & spike_burst_dyn_table_stim.fi==fi & spike_burst_dyn_table_stim.if_ff==1,:);tunnel_names(ti);
                    div_by_sr=@(x) x./fs;
                    burst_bounds=cellfun(div_by_sr,burst_bounds,'UniformOutput',false);
                    for bounds=1:length(burst_bounds)
                        for bursts=1:length(burst_bounds{bounds})
                            f=[1 2 3 4];
                            v=[burst_bounds{bounds}(bursts,1) y_low;burst_bounds{bounds}(bursts,2) y_low;burst_bounds{bounds}(bursts,2) y_hi;burst_bounds{bounds}(bursts,1) y_hi];
                            patch('Faces',f,'Vertices',v,'FaceColor',hex2rgb(char(reds(ai))),'FaceAlpha',0.3)
                        end
                    end
                end

            end

            for ai=1:length(allregion_subregion.up_fb{elecs_idx})

                spikes=allregion_subregion.up_fb{elecs_idx}{ai};
                plot(t(spikes),data(spikes),'.','Color',hex2rgb(char(blues(ai))))

                if any(contains(spike_burst_dyn_table_stim.channel_name(contains(spike_burst_dyn_table_stim.channel_name,tunnel_names(ti))...
                        & spike_burst_dyn_table_stim.fi==fi & spike_burst_dyn_table_stim.if_ff==0,:),tunnel_names(ti)))
                    burst_bounds=spike_burst_dyn_table_stim.BurstBounds(contains(spike_burst_dyn_table_stim.channel_name,tunnel_names(ti))...
                        & spike_burst_dyn_table_stim.fi==fi & spike_burst_dyn_table_stim.if_ff==0,:);tunnel_names(ti);
                    div_by_sr=@(x) x./fs;
                    burst_bounds=cellfun(div_by_sr,burst_bounds,'UniformOutput',false);
                    for bounds=1:length(burst_bounds)
                        for bursts=1:length(burst_bounds{bounds})
                            f=[1 2 3 4];
                            v=[burst_bounds{bounds}(bursts,1) y_low;burst_bounds{bounds}(bursts,2) y_low;burst_bounds{bounds}(bursts,2) y_hi;burst_bounds{bounds}(bursts,1) y_hi];
                            patch('Faces',f,'Vertices',v,'FaceColor',hex2rgb(char(blues(ai))),'FaceAlpha',0.3)
                        end
                    end
                end

            end

            hold off

            saveas(gcf,"D:\Brewer lab data\HFS\Tunnel Vis\HFS5 7SD\"+string(fi)+" "+string(elecs(1))+" "+subregions(subregion_OI)+".fig")

            % down tunnel
            figure('units','normalized','outerposition',[0 0 1 1])
            load(stim_dir+stim_folders{fi}+"\"+elecs(2)+".mat")
            data=bandpass(data,[300,3000],fs);
            plot(t,data,'k')
            hold on

            for ai=1:length(allregion_subregion.down_ff{elecs_idx})

                spikes=allregion_subregion.down_ff{elecs_idx}{ai};
                plot(t(spikes),data(spikes),'.','Color',hex2rgb(char(reds(ai))))

                if any(contains(spike_burst_dyn_table_stim.channel_name(contains(spike_burst_dyn_table_stim.channel_name,tunnel_names(ti))...
                        & spike_burst_dyn_table_stim.fi==fi & spike_burst_dyn_table_stim.if_ff==1,:),tunnel_names(ti)))
                    burst_bounds=spike_burst_dyn_table_stim.BurstBounds(contains(spike_burst_dyn_table_stim.channel_name,tunnel_names(ti))...
                        & spike_burst_dyn_table_stim.fi==fi & spike_burst_dyn_table_stim.if_ff==1,:);tunnel_names(ti);
                    div_by_sr=@(x) x./fs;
                    burst_bounds=cellfun(div_by_sr,burst_bounds,'UniformOutput',false);
                    for bounds=1:length(burst_bounds)
                        for bursts=1:length(burst_bounds{bounds})
                            f=[1 2 3 4];
                            v=[burst_bounds{bounds}(bursts,1) y_low;burst_bounds{bounds}(bursts,2) y_low;burst_bounds{bounds}(bursts,2) y_hi;burst_bounds{bounds}(bursts,1) y_hi];
                            patch('Faces',f,'Vertices',v,'FaceColor',hex2rgb(char(reds(ai))),'FaceAlpha',0.3)
                        end
                    end
                end

            end

            for ai=1:length(allregion_subregion.down_fb{elecs_idx})

                spikes=allregion_subregion.down_fb{elecs_idx}{ai};
                plot(t(spikes),data(spikes),'.','Color',hex2rgb(char(blues(ai))))

                if any(contains(spike_burst_dyn_table_stim.channel_name(contains(spike_burst_dyn_table_stim.channel_name,tunnel_names(ti))...
                        & spike_burst_dyn_table_stim.fi==fi & spike_burst_dyn_table_stim.if_ff==0,:),tunnel_names(ti)))
                    burst_bounds=spike_burst_dyn_table_stim.BurstBounds(contains(spike_burst_dyn_table_stim.channel_name,tunnel_names(ti))...
                        & spike_burst_dyn_table_stim.fi==fi & spike_burst_dyn_table_stim.if_ff==0,:);tunnel_names(ti);
                    div_by_sr=@(x) x./fs;
                    burst_bounds=cellfun(div_by_sr,burst_bounds,'UniformOutput',false);
                    for bounds=1:length(burst_bounds)
                        for bursts=1:length(burst_bounds{bounds})
                            f=[1 2 3 4];
                            v=[burst_bounds{bounds}(bursts,1) y_low;burst_bounds{bounds}(bursts,2) y_low;burst_bounds{bounds}(bursts,2) y_hi;burst_bounds{bounds}(bursts,1) y_hi];
                            patch('Faces',f,'Vertices',v,'FaceColor',hex2rgb(char(blues(ai))),'FaceAlpha',0.3)
                        end
                    end
                end

            end

            hold off

            saveas(gcf,"D:\Brewer lab data\HFS\Tunnel Vis\HFS5 7SD\"+string(fi)+" "+string(elecs(2))+" "+subregions(subregion_OI)+".fig")
        end

    end
end