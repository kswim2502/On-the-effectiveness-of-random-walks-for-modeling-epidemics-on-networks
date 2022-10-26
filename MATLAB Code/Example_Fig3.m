close all; clear; clc;
%To obtain the figure in the article, we need to adjust the scale bar with
%max 100 in 'Matrix of MITs'
A=LatticeGraph(3,4);
plot(graph(A),"NodeLabel",{})

p = 0.1;
[~,M_sym]=MFPT(A,p);

M_mit = MITmatrix(A,p);

t = tiledlayout(1,2,'TileSpacing','Compact','Padding','none');
nexttile
heatmap(M_sym)
colormap jet
title('Matrix of MFPTs for T')

nexttile
heatmap(M_mit)
colormap jet
title('Matrix of MITs') %To obtain the figure in the article, we need to adjust the scale bar with max 100
