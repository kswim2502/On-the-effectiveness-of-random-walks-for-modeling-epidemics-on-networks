function A=LatticeGraph(m,n)
    if m==1
        A1=[1];
    else
        A1=diag(ones(m-1,1),1)+diag(ones(m-1,1),-1);
    end
    if n==1
        A2=[1];
    else
        A2=diag(ones(n-1,1),1)+diag(ones(n-1,1),-1);
    end
    
    A=kron(A1,eye(n))+kron(eye(m),A2);
end