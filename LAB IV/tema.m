%% 5.1
% close all;
% clc;
% A = [1 -1.5 0.7];
% B = [1 0.5];
% C = [1 -1 0.2];
% nk = 1;
% N = 1000;
% sigma = 1;
% lambda = 1;
% 
% [Mid,Did,Dva] = ISLAB_5A(A,B,C,nk,N,sigma,lambda);

%% 5.2
close all;
clc;
A = [1 -1.5 0.7];
B = [1 0.5];
C = [1 -1 0.2];
nk = 1;
N = 100000;
sigma = 1;
lambda = 1;

[Mid_MVI,Did_MVI,Dva_MVI] = ISLAB_5B(A,B,C,nk,N,sigma,lambda);