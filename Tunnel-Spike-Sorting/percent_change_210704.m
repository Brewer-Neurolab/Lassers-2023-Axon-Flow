function percent_change_210704(table1, table2, table3, table4, subregion)

% Creates pre stim post stim dot charts for upreg and down reg for either
% ff or fb axons and plots percent change
% Table 1/3 correspond to prestim and respectively correlate to table 2/4,
% ie table 1 is prestim ff table 2 is post stim ff, table 3 is prestim fb
% table 4 post stim fb

% table1_temp=table1;
% table2_temp=table2;
% table3_temp=table3;
% table4_temp=table4;

% for i=1:length(table1)
%     [rwtb,~]=size(table1{i});
%     for j=1:rwtb
%         stim_greater=table1{i}.("Spike Area"){j}<=table2{i}.("Spike Area"){j};
%         table1{i}.("Spike Area"){j}=table1{i}.("Spike Area"){j}(stim_greater);
%         table1{i}.("Conduction Time"){j}=table1{i}.("Conduction Time"){j}(stim_greater);
%         table1{i}.("Unit Pairs"){j}=table1{i}.("Unit Pairs"){j}(stim_greater);
%         table1{i}.Direction{j}=table1{i}.Direction{j}(stim_greater);
%         
%         table2{i}.("Spike Area"){j}=table2{i}.("Spike Area"){j}(stim_greater);
%         table2{i}.("Conduction Time"){j}=table2{i}.("Conduction Time"){j}(stim_greater);
%         table2{i}.("Unit Pairs"){j}=table2{i}.("Unit Pairs"){j}(stim_greater);
%         table2{i}.Direction{j}=table2{i}.Direction{j}(stim_greater);
%         
%         stim_greater=table3{i}.("Spike Area"){j}<=table4{i}.("Spike Area"){j};
%         table3{i}.("Spike Area"){j}=table3{i}.("Spike Area"){j}(stim_greater);
%         table3{i}.("Conduction Time"){j}=table3{i}.("Conduction Time"){j}(stim_greater);
%         table3{i}.("Unit Pairs"){j}=table3{i}.("Unit Pairs"){j}(stim_greater);
%         table3{i}.Direction{j}=table3{i}.Direction{j}(stim_greater);
%         
%         table4{i}.("Spike Area"){j}=table4{i}.("Spike Area"){j}(stim_greater);
%         table4{i}.("Conduction Time"){j}=table4{i}.("Conduction Time"){j}(stim_greater);
%         table4{i}.("Unit Pairs"){j}=table4{i}.("Unit Pairs"){j}(stim_greater);
%         table4{i}.Direction{j}=table4{i}.Direction{j}(stim_greater);
%         
%     end
%     
% end

[pre_ff, post_ff]=pre_post_time_comp_210511(table1, table2);
[pre_fb, post_fb]=pre_post_time_comp_210511(table3, table4);

p1=line([0,0],[0,0],'Color',[0 0 1]);
p2=line([0,0],[0,0],'Color',[1 0 0]);

% Perc Chg line all
figure
hold on
x_pre_post=categorical({'Prestimulation' 'Poststimulation'});
x_pre_post=reordercats(x_pre_post,{'Prestimulation' 'Poststimulation'});
for i=1:length(pre_ff)
    p1=plot(x_pre_post,[0, (post_ff(i)-pre_ff(i))/pre_ff(i)],'-ob','MarkerFaceColor','b');
end
for i=1:length(pre_fb)
    p2=plot(x_pre_post,[0,(post_fb(i)-pre_fb(i))/pre_fb(i)],'-or','MarkerFaceColor','r');
end
hold off
legend([p1(1),p2(1)],'Feedforward','Feedback')
title({subregion+' Fractional Change in Conduction Spike Hz','Feedforward n= '+string(length(pre_ff)),'Feedback n= '+string(length(post_fb))})
xlabel('States of Stimulation')
ylabel('Fractional Change (Hz)')

% perc change dist
figure
hold on
for i=1:length(pre_ff)
    p1=plot((post_ff(i)-pre_ff(i))/pre_ff(i),i,'-ob','MarkerFaceColor','b');
end
for i=1:length(pre_fb)
    p2=plot((post_fb(i)-pre_fb(i))/pre_fb(i),i,'-or','MarkerFaceColor','r');
end
hold off
legend([p1(1),p2(1)],'Feedforward','Feedback')
title({subregion+' Fractional Change in Conduction Spike Hz','Feedforward n= '+string(length(pre_ff)),'Feedback n= '+string(length(post_fb))})
ylabel('Axon n')
xlabel('Fractional Change (Hz)')

% perc chg line ff
figure
hold on
for i=1:length(pre_ff)
    p1=plot(x_pre_post,[0, (post_ff(i)-pre_ff(i))/pre_ff(i)],'-ob','MarkerFaceColor','b');
end
hold off
title({subregion+' Fractional Change in Conduction Spike Hz','Feedforward n= '+string(length(pre_ff))})
xlabel('States of Stimulation')
ylabel('Fractional Change (Hz)')

% perc change dist ff
figure
hold on
for i=1:length(pre_ff)
    p1=plot((post_ff(i)-pre_ff(i))/pre_ff(i),i,'-ob','MarkerFaceColor','b');
end
hold off
title({subregion+' Fractional Change in Conduction Spike Hz','Feedforward n= '+string(length(pre_ff))})
ylabel('Axon n')
xlabel('Fractional Change (Hz)')

% perc change line fb
figure
hold on
for i=1:length(pre_fb)
    p2=plot(x_pre_post,[0,(post_fb(i)-pre_fb(i))/pre_fb(i)],'-or','MarkerFaceColor','r');
end
hold off
title({subregion+' Fractional Change in Conduction Spike Hz','Feedback n= '+string(length(post_fb))})
xlabel('States of Stimulation')
ylabel('Fractional Change (Hz)')

% perc change dist fb
figure
hold on
for i=1:length(pre_fb)
    p2=plot((post_fb(i)-pre_fb(i))/pre_fb(i),i,'-or','MarkerFaceColor','r');
end
hold off
%legend([p1(1),p2(1)],'Feedforward','Feedback')
title({subregion+' Fractional Change in Conduction Spike Hz','Feedback n= '+string(length(post_fb))})
ylabel('Axon n')
xlabel('Fractional Change (Hz)')
end