function [spike_burst_dyn_table]=compute_spike_burst_dynamics_axons_220822(allregion_unit_matched)

% This folder contains the contains the mea structs 
%load(strcat(parent_dir,'\','allregion_unit_matched'))
% Parameters
fs = 25e3; 

% Initializaing
ISI_cell={}; IBI_cell={}; IBSR_cell={}; spnb_cell={}; bd_cell={}; chan_name_cell=[];
SR_vec=[]; regi_vec=[]; fi_vec=[]; ff_vec=[];
rowi=1;
for fi=1:length(allregion_unit_matched)
    fid_allregion=allregion_unit_matched{fi};
    for regi = 1:length(unique(allregion_unit_matched{fi}.Subregion)) 
        region_allregion=fid_allregion(fid_allregion.regi==regi,:);
        electrodes_list=[];
        for elec_pairs=1:length(region_allregion.("Electrode Pairs"))
            chan_names=strsplit(region_allregion.("Electrode Pairs")(elec_pairs),"-");
            electrodes_list=[electrodes_list;chan_names'];
        end
        for chani = 1:length(region_allregion.("Electrode Pairs"))
            chani_allregion=region_allregion(region_allregion.chani==chani,:);
            chan_names = strsplit(chani_allregion.("Electrode Pairs"),"-");
            %selecting data for computation
            chan2comp = chan_names{1};
            chan_idx = find_in_channel_list(electrodes_list, chan2comp);
            timestamps=[chani_allregion.up_ff{1:end},chani_allregion.up_fb{1:end}];
            is_ff_vec=[ones(1,length(chani_allregion.up_ff{1:end})),zeros(1,length(chani_allregion.up_fb{1:end}))];

            for ui =1:length(timestamps)
                % Compute spike-burst measures
                ISI=[];
                if ~isempty(timestamps{ui})
                [ISI, ibi_sec, IBSR, spnb, bd_msec, SR_val, burst_bounds] = compute_spike_burst_dynamics_220211...
                    (timestamps{ui}');
                end
                if isempty(ISI), continue; end
                ISI_cell{rowi} = ISI;
                IBI_cell{rowi} = ibi_sec;
                IBSR_cell{rowi} = IBSR;
                spnb_cell{rowi} = spnb;
                bd_cell{rowi} = bd_msec;
                SR_vec(rowi) = SR_val;
                fi_vec(rowi) = fi; regi_vec(rowi) = regi;
                burst_bounds_cell{rowi}=burst_bounds;
                chan_name_cell = [chan_name_cell, string(chan2comp)+"-"+ui];
                ff_vec(rowi) = is_ff_vec(ui);
                rowi=rowi+1;
                disp("chan: "+chan2comp+" processed")
            end
        end
    end
    
    disp(fi+" processed")
end
spike_burst_dyn_table = table(fi_vec',regi_vec',chan_name_cell', ff_vec', SR_vec', ISI_cell', IBI_cell', IBSR_cell', spnb_cell', bd_cell',burst_bounds_cell',....
    'VariableNames',{'fi','regi','channel_name', 'if_ff', 'SpikeRate','ISI','IBI','IntraBurstSpikeRate','SpikeperBurst','BurstDuration','BurstBounds'});

end