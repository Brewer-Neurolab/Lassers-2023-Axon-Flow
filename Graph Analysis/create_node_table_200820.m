function NodeTable = create_node_table_200820
% Creates node table for the graph object with attributes Name,No, X, Y

if exist('NodeTable.mat','file')
    load('NodeTable')
    return
else
    wd = pwd;
    cd("C:\Users\lasss\Documents\Research\Brewer Lab work\Code\AFR_p1s_connections\spikes-20220228T181146Z-001\spikes")
    subregList = ["EC","DG","CA3","CA1","tunnels"];
    name_vec = {}; coordinate_vec = {};
    for mea_n = subregList
        load(char("mea_"+mea_n))
        name_vec = [name_vec; mea.channel_names];
        coordinate_vec = [coordinate_vec; mea.coordinates];
    end
    coordinate_vec = cell2mat(coordinate_vec);
    NodeTable = table((1:length(name_vec))', (name_vec), coordinate_vec(:,1),...
        coordinate_vec(:,2),'VariableNames',{'No','Name','X','Y'});
    cd(wd);
    save('NodeTable','NodeTable')
end

end
