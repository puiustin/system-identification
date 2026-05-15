%
%	File VECTORIZE.M
%
%	Function: VECTORIZE
%
%	Synopsis: v = vectorize(A) ; 
%
%	Vectorizes the matrix or vector A, by reshaping into a line vector v.
% 
%	In case of matrices, the columns are transposed and enumerated 
%	successively within the line vector. 
%
%	The output is empty if something is wrong. 
%
%	Uses:	 WAR_ERR
%
%	Author:  Dan STEFANOIU
%   2nd Revisor: Lavinius Ioan GLIGA
%	Created: March      1, 2003
%	Revised: September 22, 2007
%            February  26, 2018
%

function  v = vectorize(A)

%
% BEGIN
%
% Messages 
% ~~~~~~~~
	FN  = '<VECTORIZ>: ' ; 
	E1  = [FN 'Missing input. Empty output. Exit.'] ; 
%
% Faults preventing 
% ~~~~~~~~~~~~~~~~~
	v = [] ; 
	if (nargin < 1)
	   war_err(E1) ; 
	   return ; 
	end  
%
% Vectorizing 
% ~~~~~~~~~~~
	v = A ; 
	[m,n] = size(v) ; 
	if ((m>1) && (n>1)) 
	   v = reshape(v,1,m*n) ; 
	elseif (m>1)
	   v = v.' ; 
	end  
%
% END 
%

