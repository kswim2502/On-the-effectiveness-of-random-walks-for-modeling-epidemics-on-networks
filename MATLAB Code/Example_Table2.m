clear, clc,
C = cell(0,1);
load('ia-workplace-contacts.mat');
C{1} = A;
load('ia-contacts_hypertext2009.mat');
C{2} = A;
load('ia-contact.mat');
C{3} = A;
load('contacts-prox-high-school-2013.mat');
C{4} = A;


Table = [];
for kk = 1:length(C)
    A = C{kk};
    %p = 0.2; %infection probability
    %p = 0.1; %infection probability
    p = 0.05; %infection probability
    m = sum(A,'all')/2; %the number of edges
    n = length(A); %the number of vertices
    Repetition = 100;
    %MIT = MITmatrix(A,p);

    %%%%MCS
    tic
    MIT_Matrix_WO_GD=zeros(n); %MIT matrix without considering geometric dist.
    for MCS=1:Repetition
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
    t1 = toc

    %%%%MCS via Geometric distribution
    tic
    MIT_Matrix_GD = zeros(n);
    for ii=1:Repetition
        G = graph(C{kk});
        Samplings = geornd(p*ones(m,1))+1;
        G.Edges.Weight = Samplings;
        D = distances(G);
        MIT_Matrix_GD = MIT_Matrix_GD+D;
    end
    MIT_Matrix_GD = 1/Repetition*MIT_Matrix_GD;
    t2 = toc

    Table = [Table; t1 t2];
end