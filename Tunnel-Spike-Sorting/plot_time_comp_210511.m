function plot_time_comp_210511(spon_areas_table_ff,stim_areas_table_ff,spon_areas_table_fb,stim_areas_table_fb,subregion)
[pre_ff,post_ff]=pre_post_time_comp_210511(spon_areas_table_ff,stim_areas_table_ff);
[pre_fb,post_fb]=pre_post_time_comp_210511(spon_areas_table_fb,stim_areas_table_fb);

figure
hold on
x_pre_post=categorical({'Prestimulation' 'Poststimulation'});
x_pre_post=reordercats(x_pre_post,{'Prestimulation' 'Poststimulation'});
for i=1:length(pre_ff)
    p1=plot(x_pre_post,[pre_ff(i),post_ff(i)],'-ob','MarkerFaceColor','b');
end
for i=1:length(pre_fb)
    p2=plot(x_pre_post,[pre_fb(i),post_fb(i)],'-or','MarkerFaceColor','r');
end
hold off
legend([p1(1),p2(1)],'Feedforward','Feedback')
title({subregion+' Pre and Post Stim Change in Conduction Spike Count Per Second','Feedforward n= '+string(length(pre_ff)),'Feedback n= '+string(length(post_fb))})
xlabel('States of Stimulation')
ylabel('Spikes/Second')

figure
hold on
for i=1:length(pre_ff)
    p1=plot(x_pre_post,[pre_ff(i),post_ff(i)],'-ob','MarkerFaceColor','b');
end
hold off
title({subregion+' Feedforward Pre and Post Stim Change in Conduction Spike Count Per Second','Feedforward n= '+string(length(pre_ff))})
xlabel('States of Stimulation')
ylabel('Spikes/Second')

figure
hold on
for i=1:length(pre_fb)
    p2=plot(x_pre_post,[pre_fb(i),post_fb(i)],'-or','MarkerFaceColor','r');
end
hold off
title({subregion+' Feedback Pre and Post Stim Change in Conduction Spike Count Per Second','Feedback n= '+string(length(post_fb))})
xlabel('States of Stimulation')
ylabel('Spikes/Second')
end