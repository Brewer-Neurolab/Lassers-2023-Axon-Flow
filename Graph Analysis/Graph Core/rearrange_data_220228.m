function rearrange_data_220228(datafolder)

%datafolder = "D:\Brewer lab data\HFS\Temporal Analysis\No Stim\AFR Output\";
% Rearrange data
source_fb = load(datafolder+'FB_source_AFR_table.mat');
source_ff = load(datafolder+'FF_source_AFR_table.mat');
% Creating if_FF column
source_ff.AFR_table.if_ff = ones(height(source_ff.AFR_table),1);
source_fb.AFR_table.if_ff = zeros(height(source_fb.AFR_table),1);

AFR_table = [ source_ff.AFR_table; source_fb.AFR_table];
save(datafolder+'concatenated_source_table.mat','AFR_table')

% Loading target table
target_ff = load(datafolder+'FF_target_AFR_table.mat');
target_fb = load(datafolder+'FB_target_AFR_table.mat');
% Creating if_FF column
target_ff.AFR_table.if_ff = ones(height(target_ff.AFR_table),1);
target_fb.AFR_table.if_ff = zeros(height(target_fb.AFR_table),1);

AFR_table = [target_ff.AFR_table; target_fb.AFR_table];
save(datafolder+'concatenated_target_table.mat','AFR_table')

end
