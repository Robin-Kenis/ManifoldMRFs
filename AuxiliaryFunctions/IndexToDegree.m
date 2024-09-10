function deg = IndexToDegree(idx, K)
    % Takes an index according to lexicographical ordering of K monomials
    % and computes the (multi)degree of that monomial
    deg     = zeros(1, K);
    D       = 0;
    while idx > nchoosek(K+D, D)  
        D = D + 1;
    end
    deg(K) = D;

    if D == 0
        return
    end
    idadd = nchoosek(K+D-1, D-1);
    for i = 1:K-1
        d = 0;
        while idx > idadd + nchoosek(K-i+d, d)
            d = d + 1;
        end
        deg(i) = D - d;
        D = d;
        if D ~= 0
            idadd = idadd + nchoosek(D+K-i-1, D-1);
        end
    end
    deg(K) = deg(K) - sum(deg(1:K-1));
end