function [mu_i, mu_ij] = graph_to_measures(pg, AR, AdR, d)

    Nvert   = pg.NumNodes;
    nepg    = nodeEstimates(pg);
    nepg    = [nepg(:, 1:2), cos(nepg(:, 3)), sin(nepg(:, 3))];
    Nvar    = size(nepg, 2);
    if ~isempty(AR)
        R_i     = size(AR, 1) - size(AR, 2);
    else
        R_i     = nchoosek(2*ceil(d/2)+Nvar, Nvar);
    end
    mu_i    = zeros(R_i, Nvert);
    for i = 1:Nvert
        x0              = [1; nepg(i, :).'];
        v               = Vandermonde_vector(x0, d);
        mu_i(:, i)      = reduce_Vandermonde(v, AR);
    end

    Nedges  = pg.NumEdges;
    Edges   = edgeNodePairs(pg);
    if ~isempty(AdR)
        R_ij    = size(AdR, 1) - size(AdR, 2);
    else
        R_ij    = nchoosek(2*ceil(d/2)+2*Nvar, 2*Nvar);
    end
    mu_ij   = zeros(R_ij, Nedges);
    for ij = 1:Nedges
        i               = Edges(ij, 1);
        j               = Edges(ij, 2);
        x0              = [1; nepg(i, :).'; nepg(j, :).'];
        v               = Vandermonde_vector(x0, d);
        mu_ij(:, ij)    = reduce_Vandermonde(v, AdR);
    end

end

function F = reduce_Vandermonde(F, A)
    for i = size(A, 2):-1:1
        Acol    = A(:, i);
        AMI     = find(Acol,1,'last');
        F(AMI, :) = [];
    end
end