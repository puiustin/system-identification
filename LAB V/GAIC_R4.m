% GAIC_R4    Module that evaluates the optimum structural 
%            orders through Akaike criterion, generalized 
%            by Rissanen (GAIC_R), when identifying a 
%            model with 4 structural indexes 
%            (such as: Box-Jenkins). 
%
% Inputs:	Lambda # estimated noise variance matrix 
%               N      # size of measured data set
%
% Outputs:      nb     # optimum B polynomial order 
%               nc     # optimum C polynomial order
%               nd     # optimum D polynomial order 
%               nf     # optimum F polynomial order
%           	GAICR  # values of GAIC_R criterion
%
% Explanation:  The 4D array of estimated noise variances 
%               includes Nb rows, Nc columns, Nd layers and Nf structures
%               of size Nb,Nc,Nd. 
%               The GAIC_R criterion is evaluated for each 
%               element (i,j,k) of matrix. To select the optimum 
%               structural orders na, nb and nc, the smallest 
%               of these values is considered. Note that the 
%               structural indices result by decrementing 
%               nb, nc, nd and nf with 1. 
%
% Author:   Dan Stefanoiu (*)
% Revised:  Dan Stefanoiu (*)
%
% Created:      January  12, 2016
% Last upgrade: November 24, 2020
%
% Copyright: (*) "Politehnica" Unversity of Bucharest, ROMANIA
%                Department of Automatic Control & Computer Science
%

function [nb,nc,nd,nf,GAICR] = GAIC_R4(Lambda,N) 

%
% BEGIN
% 
% Messages 
% ~~~~~~~~
	FN  = '<GAIC_R3>: ' ; 
	E1  = [FN 'Missing, empty or null N. Empty outputs. Exit.'] ; 
	E2  = [FN 'Missing or empty Lambda. Empty outputs. Exit.'] ; 
        E3  = [FN 'No valid model detected. Empty outputs. Exit.'] ;
	W1  = [FN 'Inconsistent Lambda. Coarse structure returned.'] ; 
%
% Faults preventing
% ~~~~~~~~~~~~~~~~~
nb = [] ; 
nc = [] ;
nd = [] ;
nf = [];
GAICR = [] ;
if (nargin < 1)
   war_err(E2) ; 
   return ;
end ;
if (isempty(Lambda))
   war_err(E2) ; 
   return ;
end ; 
Lambda = abs(Lambda) ; 
Lambda(~Lambda) = eps ;  
if (nargin < 2)
   war_err(E1) ; 
   return ;
end ;
if (isempty(N))
   war_err(E1) ; 
   return ;
end ; 
N = abs(round(N(1))) ; 
if (~N)
   war_err(E1) ; 
   return ;
end ; 
% 
% Evaluating the GAIC_R criterion
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
[Nb,Nc,Nd,Nf] = size(Lambda) ; 
GAICR = [] ; 
Gmin = [] ; 
for nf=1:Nf
   [Nb,Nc,Nd,G] = GAIC_R3(Lambda(:,:,:,nf)*exp((nf-1)*log(N)/N),N) ;  
   nb = [nb Nb] ; 
   nc = [nc Nc] ; 
   nd = [nd Nd] ; 
   nf = [nf Nf] ; 
   GAICR = cat(4,GAICR,G) ; 
   Gmin = [Gmin G(Nb,Nc,Nf)] ; 
end ; 
[Gmin,nf] = min(Gmin) ;  
nb = nb(nf) ; 
nc = nc(nf) ;
nd = nd(nf) ; 
%
% END
%
