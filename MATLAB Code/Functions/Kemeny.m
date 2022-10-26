function k = Kemeny(A,type)

if (nargin <2) || isempty(type)
    type = 1;
end

if type == 1
    n=size(A,1);
    D=diag(sum(A));
    A=D^(-1)*A;
    a=eig(eye(n)-A);
    a=sort(a);
    k=sum(ones(n-1,1)./a(2:n));
else
    n=size(A,1);
    a=eig(eye(n)-A);
    a=sort(a);
    k=sum(ones(n-1,1)./a(2:n));
end
end