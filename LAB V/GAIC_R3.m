% GAIC_R3    Module that evaluates the optimum structural 
%            orders through Akaike criterion, generalized 
%            by Rissanen (GAIC_R), when identifying a 
%            model with 3 structural indexes 
%            (such as: ARMAX or incomplete Box-Jenkins). 
%
% Inputs:	Lambda # estimated noise variance matrix 
%               N      # size of measured data set
%
% Outputs:      na     # optimum row order
%           	nb     # optimum column order 
%               nc     # optimum layer order
%           	GAICR  # values of GAIC_R criterion
%
% Explanation:  The 3D array of estimated noise variances 
%               includes Na rows, Nb columns and Nc layers. 
%               The GAIC_R criterion is evaluated for each 
%               element (i,j,k) of matrix. To select the optimum 
%               structural orders na, nb and nc, the smallest 
%               of these values is considered. Note that the 
%               structural indices result by decrementing 
%               na, nb and nc with 1. 
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

function [na,nb,nc,GAICR] = GAIC_R3(Lambda,N) 

%
% BEGIN
% 
% Messages 
% ~~~~~~~~
	FN  = '<GAIC_R3>: ' ; 
	E1  = [FN 'Missing, empty or null N. Empty outputs. Exit.'] ; 
	E2  = [FN 'Missing or empty Lambda. Empty outputs. Exit.'] ;
%
% Faults preventing
% ~~~~~~~~~~~~~~~~~
na = [] ;
nb = [] ; 
nc = [] ;
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
% 
% Evaluating the GAIC_R criterion
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
[~, ~,Nc] = size(Lambda) ; 
na = zeros(1, Nc);
nb = zeros(1, Nc);
Gmin = zeros(1, Nc);
for nc=1:Nc
   [Na,Nb,G] = GAIC_R2(Lambda(:,:,nc)*exp((nc-1)*log(N)/N),N) ; 
   na(nc) = Na ; 
   nb(nc) = Nb ; 
   GAICR = cat(3,GAICR,G) ; 
   Gmin(nc) = G(Na,Nb) ; 
end 
[~,nc] = min(Gmin) ; 
na = na(nc) ; 
nb = nb(nc) ; 
%
% END
%
