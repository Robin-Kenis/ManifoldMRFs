function E = EdgesGridGraph(N1, N2)
    HEdges  = (1:N1-1).'+(0:N2-1)*N1;       % First elements in horizontal edges
    HEdges  = [HEdges(:), HEdges(:)+1];     % Matrix containing horizontal edges
    VEdges  = (1:N1:N1*N2-N1).'+(0:N1-1);   % First elements in vertical edges
    VEdges  = [VEdges(:), VEdges(:)+N1];    % Matrix containing vertical edges
    E       = [HEdges; VEdges];             % Matrix containing all edges
end