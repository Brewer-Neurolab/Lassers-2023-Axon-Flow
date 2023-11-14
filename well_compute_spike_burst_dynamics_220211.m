function [spike_burst_dyn_table]=well_compute_spike_burst_dynamics_220211(allregion_unit_matched)

%despite the name, this function computes spike dynamics for the wells

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
    for regi = 1:4
        region_allregion=fid_allregion(fid_allregion.regi==regi,:);
        electrodes_list=region_allregion.Electrode;
        %chani_allregion=region_allregion(region_allregion.chani==chani,:);
        chan_names = electrodes_list;
        %selecting data for computation
        chan2comp = chan_names;
        %chan_idx = find_in_channel_list(electrodes_list, chan2comp);
        timestamps=region_allregion.("Spike Train");
        %is_ff_vec=[ones(1,length(chani_allregion.up_ff{1:end})),zeros(1,length(chani_allregion.up_fb{1:end}))];

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
            chan_name_cell = [chan_name_cell, string(chan2comp(ui))];
            %ff_vec(rowi) = is_ff_vec(ui);
            rowi=rowi+1;
            disp("chan: "+chan2comp(ui)+" processed")
        end
    end
    
    disp(fi+" processed")
end
spike_burst_dyn_table = table(fi_vec',regi_vec',chan_name_cell', SR_vec', ISI_cell', IBI_cell', IBSR_cell', spnb_cell', bd_cell',burst_bounds_cell',....
    'VariableNames',{'fi','regi','channel_name', 'SpikeRate','ISI','IBI','IntraBurstSpikeRate','SpikeperBurst','BurstDuration','BurstBounds'});

end