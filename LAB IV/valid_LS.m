function vi = valid_LS(Model,Data) 
%
% VALID_LS   Module that evaluates the validation index of 
%            a model, by using the specified data 
%	     (appropriate for Least Squares Method). 
%
% Inputs:	Model # IDMODEL object representing the model 
%                       to validate
%        	Data  # IDDATA object representing the validation 
%                       data set
%
% Outputs:  	vi    # validation index: 
%                        0 = invalid model
%                        1 = weak validity 
%                        2 = regular validity
%                        3 = extended validity 
%
% Explanation:  The whitening test for LS Method is employed 
%               to provide the validity index. See theory for 
%               more details. 
%
% Author:   Dan Stefanoiu (*)
% Revised:  Dan Stefanoiu (*)
%
% Last upgrade: March 13, 2004
%
% Copyright: (*) "Politehnica" Unversity of Bucharest, ROMANIA
%                Department of Automatic Control & Computer Science
%

%
% BEGIN
% 
% Messages 
% ~~~~~~~~
	FN  = '<VALID_LS>: ' ; 
	E1  = [FN 'Missing or empty Model. Empty output. Exit.'] ; 
	E2  = [FN 'Missing or empty Data. Empty output. Exit.'] ; 
% 
% Faults preventing
% ~~~~~~~~~~~~~~~~~
vi = [] ;
if (nargin < 2)
   war_err(E2) ; 
   return ; 
end 
if (isempty(Data))
   war_err(E2) ; 
   return ; 
end  
if (nargin < 1)
   war_err(E1) ; 
   return ; 
end 
if (isempty(Model))
   war_err(E1) ; 
   return ; 
end  
% 
% Evaluating the prediction errors
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
e = pe(Model,Data) ; 
% 
% Evaluating the auto-correlation sequence
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
[e,k] = xcorr(e.y,'coeff') ; 
e = e(k>=0) ; 
k = length(e) ;
e = e*sqrt(k) ; 
% 
% Testing the validity
% ~~~~~~~~~~~~~~~~~~~~
vi = [(abs(e)<=2.17) ... 
      (abs(e)<=1.96) ... 
      (abs(e)<=1.808)] ; 
vi = sum(vi)/k ; 
vi = sum(vi >= [0.97 0.95 0.93]) ; 
%
% END
%

