function Indicator = MITabsorptionGD(A,p,Repetition) %A is the adjacency and p is infection prob.
G = graph(A);
m = sum(A,'all')/2; %the number of edges
n = length(A); %the number of vertices

%%%%MIT to absorption via Geometric distribution
Indicator = zeros(n,1);
for ii=1:Repetition
    Samplings = geornd(p*ones(m,1))+1;
    G.Edges.Weight = Samplings;
    D = distances(G);
    Indicator = Indicator+max(D,[],2); %MIT matrix is symmetric
end
Indicator = 1/Repetition*Indicator;
end