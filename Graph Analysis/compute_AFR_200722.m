function [AFR_vec, time_vec] = compute_AFR_200722 (peak_train, bin_s, fs)
% Computes Average firing rate from the peak_train vector given by
% spike_train within time-bin of width bin_s (sec)

binsample = bin_s*fs;            % bin [# of samples]
bin_num = floor(length(peak_train)/binsample);
len_s = length(peak_train)/fs;
time_vec = (0.5:1:bin_num-0.5).*(len_s/(bin_num));

AFR_vec = zeros(1, bin_num);  % 60*number of bins

for j = 1:bin_num  % Fill in the bins with the proper spike rate [spikes/sec]
    spikesxbin = length(find(peak_train(((j-1)*binsample+1):(j*binsample))));
    AFR_vec(j) = spikesxbin/bin_s; % [spikes/sec]
end

end