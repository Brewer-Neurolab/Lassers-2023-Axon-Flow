function [is_pos,max_r,save_fig]=seed_tunnels_with_direction_220321(source,target,fs)
% This function uses xcorr to determine the direction of spike lags to
% determine if a well is ff or fb from a tunnel. This is used if a tunnel
% cannot be determined from NMI but an electrode is still active. This is
% to be used solely with AFR computation and spike rate correlation.

time_vec=[2/1000:2/1000:300];
save_fig=0;

full_well=ismembertol(time_vec,source/1000,0.0005,'DataScale',1);
full_tunnel=ismembertol(time_vec,target/1000,0.0005,'DataScale',1);

max_lag_ms=20; %ms
% max_lag_fs=(max_lag_ms/1000)*fs;
max_lag_fs=max_lag_ms; % time vec is going by 1/1000 steps

[r,lags]=xcorr(full_well,full_tunnel,max_lag_fs);

max_corr_idx=r==max(r);
max_corr_lag=lags(max_corr_idx);
max_r=max(r);

figure
stem(lags,r)
xlabel('Lags (ms)')
ylabel('Events')

h=adtest(r);

if length(max_corr_lag)>1 || max_r(1)<5 || h~=1
    is_pos=[];
    max_r=[];
    return
end

if max_corr_lag>0
    is_pos=1;
    save_fig=1;
elseif max_corr_lag<0
    is_pos=0;
    save_fig=1;
else
    is_pos=[]; %undefined
end
end