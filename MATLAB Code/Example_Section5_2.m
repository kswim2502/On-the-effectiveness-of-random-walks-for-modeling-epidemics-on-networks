close all, clear, clc

D = dlmread('ia-infect-dublin.mtx');
G = graph(D(1:end, 1), D(1:end, 2));
[a b] = conncomp(G);
c = find(b == max(b));
x = find(a == c);
G = subgraph(G,x);
A = full(G.adjacency); A = A - diag(diag(A));

C{1} = A;

%% Comparisons for large graphs

for i = 1:length(C)
    [c, Time] = RankingNodes(C{i},0.1,300);

    [a1 b1] = sort(c(1,:));[a2 b2] = sort(c(2,:),'descend');[a3 b3] = sort(c(3,:));[a4 b4] = sort(c(4,:));[a5 b5] = sort(sum(C{i},2));
    [~, b1] = sort(b1);[~, b2] = sort(b2);[~, b3] = sort(b3);[~, b4] = sort(b4);[~, b5] = sort(b5);


    figure
    t = tiledlayout(2,3,'TileSpacing','Compact','Padding','none');

    nexttile
    plot(graph(C{i}),'NodeLabel',{},'NodeCData',b1,'MarkerSize',3,'Layout','force');
    colorbar
    title('Kemeny')

    nexttile
    plot(graph(C{i}),'NodeLabel',{},'NodeCData',b4,'MarkerSize',3,'Layout','force');
    colorbar;
    title('RWB')

    nexttile
    plot(graph(C{i}),'NodeLabel',{},'NodeCData',b5,'MarkerSize',3,'Layout','force');
    colorbar;
    title('Node Degree')

    nexttile
    plot(graph(C{i}),'NodeLabel',{},'NodeCData',b2,'MarkerSize',3,'Layout','force');
    colorbar;
    title('MICT')

    nexttile
    plot(graph(C{i}),'NodeLabel',{},'NodeCData',b3,'MarkerSize',3,'Layout','force');
    colorbar;
    title('RWC')

    colormap jet
end



