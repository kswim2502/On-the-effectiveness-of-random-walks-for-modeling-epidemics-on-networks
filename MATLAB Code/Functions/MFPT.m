function [M_rw,M_sym]=MFPT(A,p) %M1 is the MFPT matrix for T_{RW}; M2 is the MFPT matrix for T_{sym}
n=size(A,1);
D=diag(sum(A));
T1=D^(-1)*A;
[w1,~]=eigs(T1',1); w1=abs(w1)/sum(abs(w1));
W1=diag(abs(w1));
Z = inv(eye(n) - T1 + ones(n, 1)*(w1)');
M_rw=(eye(n)-Z+ones(n)*diag(diag(Z)))*inv(W1);M_rw=M_rw-diag(diag(M_rw));

if nargin==2
    m=1/p;
else
    m=max(sum(A));
end
T2=1/m*A+diag(ones(n,1)-1/m*sum(A,2));
[w2,~]=eigs(T2',1);w2=abs(w2)/sum(abs(w2));
W2=diag(w2);
Z = inv(eye(n) - T2 + ones(n, 1)*(w2)');
M_sym=(eye(n)-Z+ones(n)*diag(diag(Z)))*inv(W2);M_sym=M_sym-diag(diag(M_sym));

%     m=max(sum(A));
%     T3=1/(m+a)*A+diag(ones(n,1)-1/(m+a)*sum(A,2));
%     [w3,~]=eigs(T3',1);w2=abs(w3)/sum(abs(w3));
%     Q3=eye(n)-T3;
%     W3=diag(abs(w3));Q3s=pinv(Q3);
%     M3=(eye(n)-Q3s+ones(n)*diag(diag(Q3s)))*inv(W3);
end