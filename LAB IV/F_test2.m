function [na,nb] = F_test2(Lambda,N) 
%
% F_TEST2    Module that evaluates the optimum structural 
%            orders through the F-test, when identifying 
%            a model with 2 structural indexes 
%            (such as ARX or ARMA). 
%
% Inputs:	Lambda # estimated noise variance matrix 
%               N      # size of measured data set
%
% Outputs:      na     # optimum row order
%           	nb     # optimum column order
%
% Explanation:  The matrix of estimated noise variances 
%               includes Na rows and Nb columns. The F-test 
%               is evaluated for each element (i,j) of matrix 
%               by considering adjacent elements (i+1,j) and 
%               (i,j+1). From the resulted values, only the 
%               ones inferior to 4/N are selected. To meet the 
%               Parsimony Principle, the biggest of these values 
%               is considered when returning the optimum 
%               structural orders na and nb. Note that the 
%               structural indices result by decrementing 
%               na and nb with 1. 
%
% Author:   Dan Stefanoiu (*)
% Revised:  Dan Stefanoiu (*)
%
% Last upgrade: March 17, 2004
%
% Copyright: (*) "Politehnica" Unversity of Bucharest, ROMANIA
%                Department of Automatic Control & Computer Science
%

%
% BEGIN
% 
% Messages 
% ~~~~~~~~
	FN  = '<F_TEST2>: ' ; 
	E1  = [FN 'Missing, empty or null N. Empty outputs. Exit.'] ; 
	E2  = [FN 'Missing or empty Lambda. Empty outputs. Exit.'] ; 
    E3  = [FN 'No valid model detected. Empty outputs. Exit.'] ;
	W1  = [FN 'Inconsistent Lambda. Coarse structure returned.'] ; 
%
% Faults preventing
% ~~~~~~~~~~~~~~~~~
na = [] ;
nb = [] ;
if (nargin < 2)
   war_err(E1) ; 
   return ;
end 
if (isempty(N))
   war_err(E1) ; 
   return ;
end  
N = abs(fix(N(1))) ; 
if (~N)
   war_err(E1) ; 
   return ;
end  
N = 4/N ; 
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
% 
% Evaluating the F-test
% ~~~~~~~~~~~~~~~~~~~~~
[na,nb] = size(Lambda) ; 
if ((na<2) && (nb<2)) 
   war_err(W1) ; 
elseif ((na<2) || (nb<2))
   if (na<2)
      Lambda = (Lambda(1:(nb-1))-Lambda(2:nb))./Lambda(2:nb) ; 
   else
      Lambda = (Lambda(1:(na-1))-Lambda(2:na))./Lambda(2:na) ;
   end  
   Lambda = Lambda - N ; 
   Lambda(Lambda>=0) = -Inf ; 
   if (na<2)
      [Lambda,nb] = max(Lambda) ; 
   else 
      [Lambda,na] = max(Lambda) ; 
   end 
   if (isinf(Lambda))
      na = [] ; 
      nb = [] ;
      war_err(E3) ; 
   end 
else 
   Na = na ; 
   Nb = nb ; 
   Fr = (Lambda(1:(na-1),:)-Lambda(2:na,:))./Lambda(2:na,:) - N ;
   Fc = (Lambda(:,1:(nb-1))-Lambda(:,2:nb))./Lambda(:,2:nb) - N ; 
   Fr(Fr>=0) = -Inf ; 
   Fc(Fc>=0) = -Inf ; 
   [Fr,i] = max(Fr) ; 
   if (Na<3)
      j = i ; 
      i = 1 ; 
   else
      [Fr,j] = max(Fr) ; 
      i = i(j) ; 
   end  
   [Fc,na] = max(Fc) ; 
   if (Nb<3)
      nb = 1 ; 
   else
      [Fc,nb] = max(Fc) ; 
      na = na(nb) ; 
   end 
   if (isinf(Fr) && isinf(Fc))
      na = [] ; 
      nb = [] ; 
      war_err(E3) ;
   elseif (Fr>Fc)
      na = i ; 
      nb = j ;
   end 
end 
%
% END
%

