close all; clear; clc;

m=6; n=6;
A=LatticeGraph(m,n);
C{1} = A;

n1=7; n2=7;
A1=LatticeGraph(n1,n1); A2=LatticeGraph(n2,n2);
A=blkdiag(A1,A2);
A(1,n1^2+1)=1;A(n1^2+1,1)=1; A(2,n1^2+2)=1;A(n1^2+2,2)=1;
C{2} = A;

n1=3; n2=10;
A1=LatticeGraph(n1,n1); A2=LatticeGraph(n2,n2);
A=blkdiag(A1,A2);
A(1,n1^2+1)=1;A(n1^2+1,1)=1; A(2,n1^2+2)=1;A(n1^2+2,2)=1;
C{3} = A;


%%%% NOTE: To obtain the figures as in the article, we need to invert the
%%%% scale of MIT manually for each figure below
for i = 1:length(C)
    c = RankingNodes(C{i},0.1,300);
    
    figure
    t = tiledlayout(2,2,'TileSpacing','Compact','Padding','none');

    nexttile
    plot(graph(C{i}),'NodeLabel',{},'NodeCData',c(1,:),'MarkerSize',5,'Layout','force');
    colorbar
    title('Kemeny')

    nexttile
    plot(graph(C{i}),'NodeLabel',{},'NodeCData',c(4,:),'MarkerSize',5,'Layout','force');
    colorbar;
    title('RWB')

    nexttile
    plot(graph(C{i}),'NodeLabel',{},'NodeCData',c(2,:),'MarkerSize',5,'Layout','force');
    colorbar;
    title('MICT')

    nexttile
    plot(graph(C{i}),'NodeLabel',{},'NodeCData',c(3,:),'MarkerSize',5,'Layout','force');
    colorbar;
    title('RWC')

    colormap jet
end
