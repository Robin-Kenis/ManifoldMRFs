function V = Vandermonde_vector(g, d)
    
    n = length(g);
    K = zeros(d, n-1);
    V = g;
    for di = 2:d
        T = V*g.';
        K(di, 1) = sum(K(di-1, :)) + 1;
        for ni = 2:n-1
            K(di, ni) = K(di, 1) - sum(K(di-1, 1:ni-1));
        end
        for ni = 2:n
            V = [V; T(sum(K(di, 1:ni-1))+1:end, ni)];
        end
    end
    
end



% t = 0:0.01:2*pi;
% L1 = dual(F(10));
% p1 = cellfun(@(X) sum(full(X).*L1, 'all'), Amats_ij);
% z  = zeros(length(t), 1);
% for i = 1:size(t, 2)
% g = [1; sin(t(i)); cos(t(i)); 1];
% V = Vandermonde_vector(g, 16);
% z(i) = p1.'*V;
% end
% plot(t, z);

% t = 0:0.01:2*pi;
% L1 = dual(F(10));
% p1 = cellfun(@(X) sum(full(X).*L1, 'all'), Amats_ij);
% z  = zeros(length(t), 1);
% for i = 1:size(t, 2)
% g = [1; sin(t(i)); 0; cos(t(i)); -sin(t(i)); 0; cos(t(i))];
% V = Vandermonde_vector(g, 2);
% z(i) = p1.'*V;
% end
% plot(t, z);