function M = Algorithm1(A,p,repetition) %A is the adajcency matrix, p is infection probability
n=length(A);
M=zeros(n);%our desired mean first "infected" time matrix
for MCS=1:repetition 
    Temp=[];
    for i=1:n
        t=1;
        inf{t}=zeros(n,1);%inf{t} is a cell array each cell in which is a column vector indicating status of the individuals
        %0 means "non infected", and 1 means "infected".
        inf{t}(i)=1;%Initial infected person is $i$
        while sum(inf{t})<n % Do the following untill everyone gets infected
            inf{t+1}=inf{t};
            for j=1:n
                if inf{t}(j)==1
                    for k=1:n
                        if A(j,k)==1 & inf{t}(k)==0 & rand(1)<=p %We choose neibours of $j$ who are not infected at time t, and decide whether they are infected with probability p
                            inf{t+1}(k)=1; Temp(i,k)=t;
                        end
                    end
                end
            end
            t=t+1;
        end
    end
    M=M+Temp;
end
M=1/repetition*M;

end
