function nodes = updateNodesPositions(nodesPositions, t, nodes)

for i = 1:length(nodes)
    nodes(i).pos = nodesPositions{nodes(i).nodeIdx}(t, :);
end

end