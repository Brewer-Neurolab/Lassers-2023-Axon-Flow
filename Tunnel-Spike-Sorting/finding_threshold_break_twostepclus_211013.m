%% Peak distribution to help find largest peaks for break off
%manually load data and spikes files

min_peak_distance=1.5; %in ms
min_sd=11;

data=highpass(data,300,25000);
x=[1:length(data)];
x=x./25000;
plot(x,data,'Color','k')

[pospeaks,poslocs]=findpeaks(data,x,'MinPeakDistance',1.5/1000);
[negpeaks,neglocs]=findpeaks(-data,x,'MinPeakDistance',1.5/1000);

pos_above_thresh=find(abs(pospeaks)>=threshold);
pospeaks=pospeaks(pos_above_thresh);
poslocs=poslocs(pos_above_thresh);

neg_above_thresh=find(abs(negpeaks)>=threshold);
negpeaks=-negpeaks(neg_above_thresh);
neglocs=neglocs(neg_above_thresh);

plot(x,data,'Color','k')
hold on
plot(poslocs,pospeaks,'r.','MarkerSize',10)
plot(neglocs,negpeaks,'r.','MarkerSize',10)
hold off

%%
figure
histogram(pospeaks,30)
title('pos. peaks')
figure
histogram(-negpeaks,30)
title('neg. peaks')