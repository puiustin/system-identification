%
%	File SCALING.M
%
%	Function: SCALING
%
%	Synopsis: a = scaling(X,d) ; 
%
%	Returns the 2-length vector a that can be used to 
%	re-scale the graphical variation of X, such that 
%	the maximum/minimum values of axes are higher/lower 
%       by d*(max(X)-min(X)) comparing to max(X)/min(X). 
%	By default, d=0.1 (which means 10%). The argument 
%	d can be negative, allowing the user to perform 
%	zoom on a specified zone in the graphic. 
%
%	Usually, X is a vector. If X is a matrix, its columns 
%	are understood as different variations of the same 
%	parameter (like in case of MATLAB function PLOT) 
%	and max/min are computed accordingly. 
%
%	The output is empty if something is wrong. 
%
%	Uses:	 WAR_ERR
%
%	Author:  Dan STEFANOIU
%   Revised: Lavinius Ioan GLIGA
%	Created: April 14, 2004
%	Revised: July  31, 2007
%            August 9, 2018
%

function  a = scaling(X,d) 

%
% BEGIN
%
% Messages 
% ~~~~~~~~
	FN  = '<SCALING>: ' ; 
	E1  = [FN 'Missing or empty input. Empty output. Exit.'] ; 
%
% Faults preventing 
% ~~~~~~~~~~~~~~~~~
	a = [] ; 
	if (nargin < 2) 
	   d = 0.1 ; 
	end 
	if (isempty(d))
	   d = 0.1 ; 
	end  
	d = d(1) ; 
	if (nargin < 1)
	   war_err(E1) ; 
	   return ; 
	end 
	if (isempty(X))
	   war_err(E1) ; 
	   return ; 
	end 
%
% Scaling 
% ~~~~~~~
	if (norm(X)<eps)
	   a = [-1 1] ; 
	else
	   a = min(min(X)) ; 
	   X = max(max(X)) ; 
	   d = d*(X-a) ; 
	   a = [a-d X+d] ; 
	end 
%
% END 
%
