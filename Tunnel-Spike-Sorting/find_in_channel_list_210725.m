function [idx_vec] = find_in_channel_list_210725(channel_list,chan_name_vec)
%Returns the index of chan_name in channel_list
channel_list = replace(channel_list, ' ','');
if class(chan_name_vec) == "char"
    chan_name_vec = {char(chan_name_vec)};
end

for chani = 1:length(chan_name_vec)
    chan_name = replace(chan_name_vec{chani}, ' ','');
    idx = find(strcmpi(cellstr(channel_list), chan_name));
    if ~isempty(idx)
        idx_vec(chani) = idx;
    else
        idx_vec(chani) = 0;
    end
end

end