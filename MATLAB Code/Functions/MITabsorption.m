function [Indicator,Meantime_absorption,T]=MITabsorption(A,p) %Probability Infection Matrix
P=p*A;
n=length(A);

C = cell(1,0);
for k=1:n
    X = dec2bin(sum(nchoosek(2.^(0:n-1),k),2)) - '0';
    C{length(C)+1} = fliplr(X);
end
ver_binary = C;

% Note: we shall use binary strings as vertices of our transition matrix.
% We suppose they are given according to the number of ones
% 'ver_binary' is a cell array that contains $n$ cells; and
% 'ver_binary{i}' is the matrix each row of which is a $(0,1)$ vector with $i$ ones

T=[]; % this is our desired transition matrix
a0=0; % will be used for parameter to construct $T$
b0=2^n-1;% note: we exclude 0...0

% Here we construct $T$ row-block by row-block according to sizes of cells
% ii is row block index, and jj is column block index
% For each block in a row block, kk and ll is for rows and column, resp.
for ii=1:n
    T_blockrow=[];
    Temp1=ver_binary{ii};
    for jj=ii+1:n
        Temp2=ver_binary{jj};
        T_blockrow_0=zeros(size(Temp1,1),size(Temp2,1));
        for kk=1:size(Temp1,1)
            x_k=Temp1(kk,:);
            for ll=1:size(Temp2,1)
                x_l=Temp2(ll,:);
                if x_k*x_l'==nnz(x_k) & nnz(x_k)~=nnz(x_l) %Note if this condition is not satisfied, the corresponding entry is 0
                    Probability=1;
                    z0=find(x_l-x_k==1);%z0 is uninfected individuals in x_l in "previous" step
                    for qq=1:length(z0)
                        Temp_prob=x_k.*P(z0(qq),:);
                        Temp_prob=1-Temp_prob(Temp_prob~=0);
                        Probability=(Probability)*(1-prod(Temp_prob));
                        %Probability=(Probability)*(1-(1-p)^(x_k*A(z0(qq),:)'));%consider the complementary event
                    end
                    z1=find(x_l==0);%z1 is uninfected individuals in x_l in "previous" and "current" steps
                    for rr=1:length(z1)
                        Temp_prob=x_k.*P(z1(rr),:);
                        Temp_prob=1-Temp_prob(Temp_prob~=0);
                        Probability=(Probability)*prod(Temp_prob);
                        %Probability=(Probability)*(1-p)^(x_k*A(z1(rr),:)');
                    end
                    T_blockrow_0(kk,ll)=Probability;
                end
            end
        end
        T_blockrow=[T_blockrow T_blockrow_0];
    end
    T(a0+1:a0+size(Temp1,1),a0+size(Temp1,1)+1:b0)=T_blockrow;
    a0=a0+size(Temp1,1);
end

Q=T+diag(ones(b0,1)-sum(T,2));
T=Q;

T_temp=eye(b0)-Q;T_temp(:,end)=[];T_temp(end,:)=[];
Meantime_absorption=linsolve(T_temp,ones(length(T_temp),1));
Indicator=Meantime_absorption(1:n);
end

