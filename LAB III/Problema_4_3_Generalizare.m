%% ARX[na,nb]

at = [ 0.2, 0.3,  0.5, 0.6];  % include 1 pe prima poziție
bt = [0.1, 0.3, -0.2, 0.2, 0.3]; % 0 la început pentru delay nk=1
na = length(at);  % numărul de coef AR (fără primul 1)
nb = length(bt);    % numărul de coef X
nk = 2;
[a_I,b_I,lambda_I,magi_I,phii_I,mag_I,phi_I,f_I] = ISLAB_4I(na,nb,nk,at,bt,[],[],[]); 
%% 
at1 = [0.2, 0.3, -0.1, 0.5, 0.6];  % include 1 pe prima poziție
bt1 = [0.4, 0.7, -0.2, 0.2, 0.3]; % 0 la început pentru delay nk=1
na1 = length(at1);  % numărul de coef AR (fără primul 1)
nb1 = length(bt1);    % numărul de coef X
nk1 = 0;
[a_I1,b_I1,lambda_I1,magi_I1,phii_I1,mag_I1,phi_I1,f_I1] = ISLAB_4J(na1,nb1,nk1,at1,bt1,[],[],[]); 