% GAIC_R2    Module that evaluates the optimum structural 
%            orders through Akaike criterion, generalized 
%            by Rissanen (GAIC_R), when identifying a 
%            model with 2 structural indexes 
%            (such as: ARX or ARMA). 
%
% Inputs:	Lambda # estimated noise variance matrix 
%               N      # size of measured data set
%
% Outputs:      na     # optimum row order
%           	nb     # optimum column order 
%           	GAICR  # values of GAIC_R criterion
%
% Explanation:  The matrix of estimated noise variances 
%               includes Na rows and Nb columns. The GAIC_R 
%               criterion is evaluated for each element (i,j) 
%               of matrix. To select the optimum structural 
%               orders na and nb, the smallest of these values 
%               is considered. Note that the structural indices 
%               result by decrementing na and nb with 1. 
%
% Author:   Dan Stefanoiu (*)
% Revised:  Dan Stefanoiu (*)
%           Lavinius Ioan Gliga (*)
%
% Created: January  12, 2016
% Revised: November 24, 2020
%          August    9, 2018
%
% Copyright: (*) "Politehnica" Unversity of Bucharest, ROMANIA
%                Department of Automatic Control & Computer Science
%

function [na,nb,GAICR] = GAIC_R2(Lambda,N) 

%
% BEGIN
% 
% Messages 
% ~~~~~~~~
	FN  = '<GAIC_R2>: ' ; 
	E1  = [FN 'Missing, empty or null N. Empty outputs. Exit.'] ; 
	E2  = [FN 'Missing or empty Lambda. Empty outputs. Exit.'] ;
	W1  = [FN 'Inconsistent Lambda. Coarse structure returned.'] ; 
%
% Faults preventing
% ~~~~~~~~~~~~~~~~~
na = [] ;
nb = [] ; 
GAICR = [] ; 
if (nargin < 1)
   war_err(E2) ; 
   return ;
end
if (isempty(Lambda))
   war_err(E2) ; 
   return ;
end 
Lambda = abs(Lambda) ; 
Lambda(~Lambda) = eps ; 
if (nargin < 2)
   war_err(E1) ; 
   return ;
end 
if (isempty(N))
   war_err(E1) ; 
   return ;
end 
N = abs(round(N(1))) ; 
if (~N)
   war_err(E1) ; 
   return ;
end 
N = log(N)/N ; 
% 
% Evaluating the GAIC_R criterion
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
[na,nb] = size(Lambda) ; 
GAICR = log(Lambda) + ... 
        N*((1:na)'*ones(1,nb) + ones(na,1)*(1:nb)-2) ; 
if ((na<2) && (nb<2)) 
   war_err(W1) ; 
elseif (na<2) 
   [~,nb] = min(GAICR) ; 
elseif (nb<2) 
   [~,na] = min(GAICR) ; 
else
   [Lambda,na] = min(GAICR) ; 
   [~,nb] = min(Lambda) ; 
   na = na(nb) ;
end 
%
% END
%
