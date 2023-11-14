function plot_tunnels_in_graph_220309(unit_match_struct, endNodes, isSlope)

get_well_node_no = @(NodeTable, well_name) find(strcmp(NodeTable.Name, well_name));

tunnelEndNodes=string(endNodes.Edges.EndNodes);

if isSlope
    max_lwidth = 2;
else
    max_lwidth = 1;
end

%% Create NodeTable, and EdgeTable
NodeTable = create_node_table_200820;
diEndNodes = []; diCode =[];
unEndNodes = []; unCode =[];
edge_line_widths=[];
for regi=1:4
    for chani=1:5
        split_channel_names = strsplit(unit_match_struct{regi}{chani,2},'-');
        
        % Getting upstream and downstream channel names
        % Edited by SBL to only have up channel be the first listed in name
        % be up channel
        up_chan_name = split_channel_names{1};
        down_chan_name = split_channel_names{2};
        
%         max_FF_LW=10*max(endNodes.Edges.Weight((contains(tunnelEndNodes(:,2),up_chan_name))...
%                         | (contains(tunnelEndNodes(:,1),down_chan_name))& endNodes.Edges.Code==0.1))/max_lwidth;
%         max_FB_LW=10*max(endNodes.Edges.Weight((contains(tunnelEndNodes(:,2),down_chan_name))...
%                         | (contains(tunnelEndNodes(:,1),up_chan_name)) & endNodes.Edges.Code==0.6))/max_lwidth;

        max_FF_LW=10*max(endNodes.Edges.Weight(any(strcmp(tunnelEndNodes(:,2),up_chan_name),2)...
                        | any(strcmp(tunnelEndNodes,down_chan_name),2) & endNodes.Edges.Code==0.1))/max_lwidth;
        max_FB_LW=10*max(endNodes.Edges.Weight(any(strcmp(tunnelEndNodes,down_chan_name),2)...
                        | any(strcmp(tunnelEndNodes,up_chan_name),2) & endNodes.Edges.Code==0.6))/max_lwidth;
                    
        if isempty(unit_match_struct{regi}{chani,3})                        % Empty tunnels
            unEndNodes = [unEndNodes; [get_well_node_no(NodeTable, up_chan_name), ...
                    get_well_node_no(NodeTable, down_chan_name)]];
            unCode = [unCode; 1];
        else
            isPopulated=0;
            isFirst=0;
            for ui=1:length(unit_match_struct{regi}{chani,3})
                if unit_match_struct{regi}{chani,3}(ui) ...
                        && (any(strcmp(tunnelEndNodes(:,2),up_chan_name))...
                        || any(strcmp(tunnelEndNodes(:,1),down_chan_name)))                     % FF tunnels
                    diEndNodes = [diEndNodes; [get_well_node_no(NodeTable, up_chan_name), ...
                        get_well_node_no(NodeTable, down_chan_name)]];
                    diCode = [diCode; 0.1];
                    isPopulated=1;
                    if isFirst
                        unEndNodes(end,:)=[];
                        unCode(end)=[];
                        isFirst=0;
                    end
                    edge_line_widths=[edge_line_widths;max_FF_LW];
                elseif ~unit_match_struct{regi}{chani,3}(ui) ...
                        && (any(strcmp(tunnelEndNodes(:,2),down_chan_name))...
                        || any(strcmp(tunnelEndNodes(:,1),up_chan_name)))               % FB tunnels
                    diEndNodes = [diEndNodes; [get_well_node_no(NodeTable, down_chan_name), ...
                        get_well_node_no(NodeTable, up_chan_name)]];
                    diCode = [diCode; 0.6];
                    isPopulated=1;
                    if isFirst
                        unEndNodes(end,:)=[];
                        unCode(end)=[];
                        isFirst=0;
                    end
                    edge_line_widths=[edge_line_widths;max_FB_LW];
                elseif ~isPopulated
                    unEndNodes = [unEndNodes; [get_well_node_no(NodeTable, up_chan_name), ...
                    get_well_node_no(NodeTable, down_chan_name)]];
                    unCode = [unCode; 1];
                    isPopulated=1;
                    if ui==1
                        isFirst=1;
                    end
                end
            end
        end
    end
end

diEdgeTable = table(diEndNodes, diCode, 'VariableNames',{'EndNodes', 'Code'});
if isempty(unEndNodes) 
    unEdgeTable = table(unEndNodes, [], 'VariableNames',{'EndNodes', 'Code'});
else
    unEdgeTable = table(unEndNodes, unCode, 'VariableNames',{'EndNodes', 'Code'});
end

%% Plotting Filled tunnels (directed graph)

G = digraph(diEdgeTable, NodeTable);

[~,weight_idx]=sort(diEdgeTable.EndNodes(:,1));
edge_line_widths=edge_line_widths(weight_idx);

edge_obj = plot(G,'XData',G.Nodes.X,'YData',G.Nodes.Y,...
         'NodeColor','black');
edge_obj.EdgeCData = G.Edges.Code;
edge_obj.LineWidth = edge_line_widths;

arrow_obj = plot( G,'XData',G.Nodes.X,'YData',G.Nodes.Y,...
     'NodeColor','black');
arrow_obj.LineStyle = ':';
arrow_obj.EdgeColor = 'black';
arrow_obj.ArrowSize = 20;
arrow_obj.EdgeAlpha = 0.8;

%% Plotting empty tunnels (undirected graphs)
if ~isempty(unEdgeTable)
    G = graph(unEdgeTable, NodeTable);
    edge_obj = plot( G,'XData',G.Nodes.X,'YData',G.Nodes.Y,...
             'NodeColor','black');
    edge_obj.EdgeCData = G.Edges.Code;
    edge_obj.LineWidth = 2;
    set(gca,'XTick',[],'YTick',[])
end
hold off
end

