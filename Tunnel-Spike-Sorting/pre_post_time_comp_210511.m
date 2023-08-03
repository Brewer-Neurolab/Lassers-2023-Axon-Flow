function [pre_stim_vals,post_stim_vals]=pre_post_time_comp_210511(table1,table2)
% Tables 1 and 2 should be pre stim and post stim in that order
% Assumes the tables are already zero padded
pre_stim_vals=[];
post_stim_vals=[];
for i=1:length(table1)
    pre_stim_vals=[pre_stim_vals,[table1{i}.("Spike Area"){1:end}]];
    post_stim_vals=[post_stim_vals,[table2{i}.("Spike Area"){1:end}]];
end

pre_stim_vals=pre_stim_vals./300;
post_stim_vals=post_stim_vals./300;
end