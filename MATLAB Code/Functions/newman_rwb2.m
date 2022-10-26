function RWB = newman_rwb2(A)

%Transition matrix for random walk
n=size(A, 1);
D = diag(sum(A, 2));
vol=sum(diag(D));
T = inv(D)*A;

%Mean first passage matrix
[V, E] = eig(T.');
Eigs = nonzeros(E);
one = max(real(Eigs));
i = find(Eigs == one);
w = (1/sum(V(:, i)))*V(:, i)'; %stationary vector

W = diag(w);
Z = inv(eye(n) - T + ones(n, 1)*(w)); % fundamental matrix Z

M = (eye(n) - Z + ones(n)*diag(diag(Z)))*(inv(W)); %MFP matrix
M = M-diag(diag(M)); % convention: m_ii:=0

%Computing random walk betweenness
RWB = zeros(n, 1);

for i = 1:n
    J=find(A(i, :));
    for s=1:n
        for t = (s+1):n
            if (i~=s) && (i~=t)
                RWB(i)=RWB(i)+ (0.5/vol)*(sum(abs(M(s,J)-M(s,i)+M(t,i)-M(t,J)))); %M(J,s)-M(i,s)+M(i,t)-M(J,t)
            else
                RWB(i)=RWB(i)+1;
            end
        end
    end
end

RWB=RWB/(0.5*n*(n-1));
RWB=RWB';
end