function dr2 = PrimFeas(X)
    dr2         = zeros(size(X, 3), 1);
    for dr = 1:size(X, 3)
        d           = eig(value(X(:, :, dr)));
        dmax2       = max(d)^2;
        dr2(dr)     = 1 - dmax2/norm(d, 'fro')^2;
    end
end