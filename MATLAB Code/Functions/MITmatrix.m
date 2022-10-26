function M=MITmatrix(A,p)
[~,~,T] = MITabsorption(A,p); % T is the (2^n-1)x(2^n-1) transition matrix
n = log2(length(T)+1);

S = cell(1,0);
for k=1:n
    X = dec2bin(sum(nchoosek(2.^(0:n-1),k),2)) - '0';
    S{length(S)+1} = fliplr(X);
end
S = cell2mat(S');
M = [];
for i = 1:n
    Temp = T;
    x = find(S(:,i) == 1);
    Temp(x,:) = []; Temp(:,x) = [];
    v = sum(inv(eye(length(Temp))-Temp),2);
    if i == 1
        v = [0 ; v(1:n-1)];
    elseif i == n
        v = [v(1:n-1) ; 0];
    else
        v = [v(1:i-1); 0; v(i:n-1)];
    end
    M = [M v];
end

end