function idx = DegreeToIndex(deg, D)
    % Takes a (multi)degree D and returns the index of the monomial of that
    % degree according to the lexicographical ordering
    idx     = 1;
    K       = length(deg);
    deg    = [D-sum(deg), deg];
    for i = 1:K
        D   = D - deg(i);
        if D ~= 0
            idx = idx + nchoosek(K+D-1, D-1);
        end
        K   = K - 1;
    end
end