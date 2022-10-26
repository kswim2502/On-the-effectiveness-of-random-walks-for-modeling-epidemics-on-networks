clear, clc, close all

A1 = ones(3)-eye(3); A1(4,1) = 1; A1(1,4) = 1;
C{1} = A1;

load('paley9.mat');
C{2} = paley9;

load('petersen.mat');
C{3} = petersen;

A4 = zeros(12); A4(1,2:end)=ones(1,11); A4 = A4+A4';
C{4} = A4;

A5=LatticeGraph(3,4);
C{5} = A5;

Table = [];
for kk = 1:length(C)
    A = C{kk};
    p = 0.01; %infection probability
    m = sum(A,'all')/2; %the number of edges
    n = length(A); %the number of vertices
    Repetition = 300;
    MIT = MITmatrix(A,p);

    %%%%MCS
    tic
    MIT_Matrix_WO_GD=zeros(n); %MIT matrix without considering geometric dist.
    for MCS=1:1:Repetition
        Temp=[];
        for i=1:n
            t=1;
            inf{t}=zeros(n,1);
            inf{t}(i)=1; % i is infected at first
            while sum(inf{t})<n % Do the following untill everyone gets infected.
                inf{t+1}=inf{t};
                x = find(inf{t}==1); % infected individuals at time t
                for j=1:length(x)
                    y = find(A(x(j),:)==1); % neighbours of infected people
                    for k=1:length(y)
                        if inf{t}(y(k))==0 && rand(1)<=p
                            inf{t+1}(y(k))=1; Temp(i,y(k))=t;
                        end
                    end
                end
                t=t+1;
            end
        end
        MIT_Matrix_WO_GD=MIT_Matrix_WO_GD+Temp;
    end
    MIT_Matrix_WO_GD=1/Repetition*MIT_Matrix_WO_GD;
    t1 = toc;

    %%%%MCS via Geometric distribution
    tic
    MIT_Matrix_GD = zeros(n);
    for ii=1:Repetition
        G = graph(C{kk});
        Samplings = geornd(p*ones(m,1))+1;
        G.Edges.Weight = Samplings;
        Temp_Matrix = zeros(n);
        for i=1:n-1
            for j=i+1:n
                [~,w] = shortestpath(G,i,j);
                Temp_Matrix(i,j) = w;
            end
        end
        MIT_Matrix_GD = MIT_Matrix_GD+Temp_Matrix;
    end
    MIT_Matrix_GD = 1/Repetition*(MIT_Matrix_GD+MIT_Matrix_GD');
    t2 = toc;

    Mu1 = mean(abs(MIT-MIT_Matrix_WO_GD)./MIT, 'all', 'omitnan');
    Sig1 = var(abs(MIT-MIT_Matrix_WO_GD)./MIT, 0, 'all', 'omitnan');

    Mu2 = mean(abs(MIT-MIT_Matrix_GD)./MIT, 'all', 'omitnan');
    Sig2 = var(abs(MIT-MIT_Matrix_GD)./MIT, 0, 'all', 'omitnan');

    Table = [Table; Mu1 Sig1 Mu2 Sig2];
end