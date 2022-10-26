function [Ranking, Time]=RankingNodes(A,p,Repetition)
%%Kemeny's centrality based on node removal
Ind1=[];
tic
for i=1:length(A)
    B1=A;
    B1(i,:)=[]; B1(:,i)=[]; %the submatrix by deleting i-row and i-column
    G1 = graph(B1); [~,a] = conncomp(G1);
    if length(a) == 1
        Ind1=[Ind1 Kemeny(B1)];
    else
        Ind1=[Ind1 1e+10];
    end
end
t1 = toc;

%% MITs to absorption
tic
MITabsorption = MITabsorptionGD(A,p,Repetition);
t2 = toc;

%% RWC
tic
Acc = sum(A,2)'*ERM(A)-Kemeny(A); %x is a row vector whose entries are RWCs
RWC = 1./Acc;RWC=round(RWC,7);
t3 = toc;

%% RWB
tic
RWB = [];
RWB = newman_rwb2(A);RWB=round(RWB,7);
t4 = toc;
%% ranking
Ranking = [round(Ind1,7); round(MITabsorption',7); RWC; RWB];
Time = [t1 t2 t3 t4];

end