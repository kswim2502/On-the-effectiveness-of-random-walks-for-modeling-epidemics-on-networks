function R=ERM(A)
    n=length(A);
    L=diag(sum(A))-A;
    MPL=pinv(L);%Moore--Penrose inverse
    R=zeros(n);
    for i=1:n
       for j=1:n
          R(i,j)=MPL(i,i)+MPL(j,j)-MPL(i,j)-MPL(j,i);
       end
    end
end

