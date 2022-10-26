clear all
clc
A=ones(3)-eye(3);A(1,4)=1;A(4,1)=1;
plot(graph(A))

[M_rw,M_sym]=MFPT(A);

p=0.1;
repetition=300;
M_algorithm1 = Algorithm1(A,p,repetition);

M_Actual = MITmatrix(A,p);

[Indicator,Meantime_absorption,T]=MITabsorption(A,p);