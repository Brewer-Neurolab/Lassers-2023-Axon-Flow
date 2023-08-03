function [binCounts, NMI, xbin_ISI,b, keep_time ] = get_spike_delay_hist_210712(down_channel, up_channel, if_flip, if_hist, length_rec_sec)
%creates ISI histogram between ISI of cw_channel and ccw_channel

%Changed from 0.1 to 0.2 1/21/21
min_difference = 0.1; %ms ,ie, velocity = 2 m/s
max_difference = 1; %ms ,ie, velocity = 0.2 m/s
ISI = [];
unpaired = 0;
paired_spikes = 0;
ambig = 0;  
Histogram = 0;

%keep time structure stores ISI and up/down channel times
keep_time.ISI=[];
keep_time.up=[];
keep_time.down=[];
up_keep=[];
down_keep=[];
keep_time.edges=[];

%Not indicative of any real value, evaluate conduction time histograms
%yourself, 20 bins from this 40kHz measure is fine
% fs = 40e3;
fs = 25e3;

if if_flip
    temp = down_channel;
    down_channel = up_channel;
    up_channel = temp;
end
tot_spikes_up = length(up_channel);


%calculating interspike intervals  
for n = 1:tot_spikes_up  
    ISI_entire_record = up_channel(n,1) - down_channel(:,1);
    boundary_for_pair = find(ISI_entire_record > -max_difference & ISI_entire_record < max_difference);              %max limits on difference
        if length(boundary_for_pair) > 1                                                                              %more than 1 paired spike (Yash code new)
            ambig = ambig +1;
            ISI = [ISI NaN ];
            up_keep=[up_keep NaN ];
            down_keep=[down_keep NaN ];
        elseif isempty(boundary_for_pair)                                                                               %unpaired spike
            unpaired = unpaired +1;
            ISI = [ISI NaN ];
            up_keep=[up_keep NaN ];
            down_keep=[down_keep NaN ];
        else
            if ISI_entire_record(boundary_for_pair)> min_difference || ISI_entire_record(boundary_for_pair) < -min_difference  %min limit on difference 
                paired_spikes = paired_spikes + 1;
                ISI = [ISI ISI_entire_record(boundary_for_pair)];
                up_keep=[up_keep up_channel(n,1) ];
                down_keep=[down_keep down_channel(boundary_for_pair,1) ];
            else
                ambig = ambig+1;                                                                                        %ISI < min_difference
                ISI = [ISI NaN ];
                up_keep=[up_keep NaN ];
                down_keep=[down_keep NaN ];
            end
        end
end

binw = 2*1e3/fs;                                                           % Change bin-width here
                                                                           % For 40kHz, bin-width is 0.025 ms

if if_hist
    ISI = -ISI;                                                            %making positive lag means FF channels leading
    xbin_ISI = (-1+binw/2):binw:(1-binw/2);
%     Histogram = hist(ISI,xbin_ISI);
    %% Dividing by no. of sec
    xbin_edges = convertCenters2edges(xbin_ISI);
    Histogram = histogram(ISI, xbin_edges,'Visible','off');
    binCounts = Histogram.BinCounts;
    binCounts = binCounts./(length_rec_sec);
    b=bar(xbin_ISI,binCounts);
%     hist(ISI,xbin_ISI);
    set(gca,'fontsize',13)
    xlabel('<- Feedback            Conduction time (ms)         Feedforward->');
    ylabel('No. of Spike Pairs');
    xticks(xbin_ISI)
end

NMI = paired_spikes./max([length(down_channel), length(up_channel)]);

keep_time.ISI=ISI;
keep_time.up=up_keep;
keep_time.down=down_keep;
xbin_edges=xbin_edges(xbin_edges~=0);
keep_time.edges=xbin_edges;

end
