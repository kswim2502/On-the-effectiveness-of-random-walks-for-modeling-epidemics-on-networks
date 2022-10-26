a = [];
n = 10;
p = 0.5;

for i = 1:100
    G = graph(true(n), 'omitselfloops');
    edgesToKeep = rand(numedges(G), 1) < p;
    G = graph(G.Edges(edgesToKeep, :));
    A = full(adjacency(G));
    [~,~,T] = MITabsorption(A,0.1);
    a = [a nnz(T)/length(T)^2];
end
