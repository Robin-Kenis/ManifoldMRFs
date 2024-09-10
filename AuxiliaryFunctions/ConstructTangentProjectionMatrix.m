function Proj = ConstructTangentProjectionMatrix(Nvar, Dvar, D, G)
    % Takes a matrix G encoding the projection onto the tangent space of
    % some manifold (as well as the dimensions of the embedding space and 
    % said tangent space, resp Nvar and Dvar, and the degree D) and returns
    % a matrix projecting onto the tangent space of the manifold, used in
    % encoding the Lipschitz constraint
    G       = G.';
    P       = nchoosek(D + Nvar, Nvar);
    N       = nchoosek(D + Nvar, Nvar);
    [GM, W] = GradientMatrix(D, Nvar);
    GM      = GM(:, 1:Nvar);
    W       = W(:, 1:Nvar);
    
    ind1    = sparse([], [], [], P, Dvar);
    ind2    = sparse([], [], [], P, Dvar);
    vals    = sparse([], [], [], P, Dvar);
    kidx    = 1;
    for h = 1:size(GM, 1)
        for i = 1:Nvar
            dind = GM(h, :);
            if dind(i) ~= 0
                deg = IndexToDegree(dind(i), Nvar);
                for j = 1:Dvar
                    if G(i, j) ~= 0
                        deg(abs(G(i, j)))   = deg(abs(G(i, j))) + 1;
                        inddeg              = DegreeToIndex(deg, 2*ceil(D/2));
                        ind1(kidx, j)       = inddeg;
                        ind2(kidx, j)       = h;
                        vals(kidx, j)       = sign(G(i, j))*W(h, i); 
                        deg(abs(G(i, j)))   = deg(abs(G(i, j))) - 1;
                        kidx = kidx + 1;
                    end
                end
            end
        end
    end
    Proj = cell(Dvar, 1);
    for j = 1:Dvar
        bools = ind1(:, j) ~= 0;
        Proj{j} = sparse(ind1(bools, j), ind2(bools, j), vals(bools, j), P, N);
    end
end