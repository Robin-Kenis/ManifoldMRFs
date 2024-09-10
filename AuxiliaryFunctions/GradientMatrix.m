function [GM, W] = GradientMatrix(D, N)
    P       = nchoosek(D+N, N);
    GM      = zeros(P, N);
    W       = zeros(P, N);
    for i = 1:N
        for j   = 1:P
            deg = IndexToDegree(j, N);
            if deg(i) > 0
                W(j, i)         = deg(i);
                deg(i)          = deg(i) - 1;
                k               = DegreeToIndex(deg, D);
                GM(j, i)        = k;
            end
        end
    end
end