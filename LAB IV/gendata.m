function [D,V,P] = gendata(A,B,C,nk,N,sigma,lambda) 
%
% GENDATA    Module that generates data from an 
%            ARX model with colored noise 
%            (or ARMAX model). 
%
% Inputs:	A      # true coefficients of AR part 
%                        ([1 -1.5 0.7], by default)
%        	B      # true coefficients of X part 
%                        ([1 0.5], by default)
%        	C      # true coefficients of MA part (noise filter)
%                        ([1 -1 0.2], by default)
%               nk     # intrinsic delay of process (1, by default)
%               N      # simulation period (250, by default)
%               sigma  # standard deviation of PRB input 
%                        (1, by default); 
%                        if null, ARMA or AR processes are considered
%               lambda # standard deviation of white noise 
%                        (1, by default)
%                        if null, noise free ARX or AR processes 
%                        are considered
%
% Outputs:  	D      # IDDATA object representing the 
%                        I/O generated data 
%               V      # IDDATA object representing the 
%                        I/O noise generated data 
%                        (white noise as input, colored noise 
%                        as output)
%               P      # IDMODEL object representing the 
%                        process that provided the data 
%
% Explanation:  An ARMAX model is stimulated with a PRB input 
%               with 2 values (+1 and -1). The generated input 
%               and the observed output are returned. Note that 
%               various processes can be used, by setting sigma 
%               and/or lambda to null. 
%
% Author:   Helena Haglund (*)
%           Bjorn Wittenmark (*)
% Revised:  Dan Stefanoiu (**)
%
% Last upgrade: (*)  January  5, 1997
%               (**) March   17, 2004
%
% Copyright: (*)  Lund Institute of Technology, SWEDEN
%                 Department of Automatic Control
%            (**) "Politehnica" Unversity of Bucharest, ROMANIA
%                 Department of Automatic Control & Computer Science
%

%
% BEGIN
% 
% Faults preventing
% ~~~~~~~~~~~~~~~~~
if (nargin < 7)
   lambda = 1 ;
end 
if (isempty(lambda))
   lambda = 1 ;
end  
lambda = abs(lambda(1)) ; 
if (nargin < 6)
   sigma = 1 ;
end 
if (isempty(sigma))
   sigma = 1 ;
end  
sigma = abs(sigma(1)) ; 
if (nargin < 5)
   N = 250 ;
end 
if (isempty(N))
   N = 250 ;
end  
N = abs(fix(N(1))) ; 
if (~N)
   N = 250 ;
end  
if (nargin < 4)
   nk = 1 ;
end 
if (isempty(nk))
   nk = 1 ;
end 
nk = abs(fix(nk(1))) ; 
if (~nk)
   nk = 1 ;
end  
if (nargin < 3)
   C = [1 -1 0.2] ;
end 
if (isempty(C))
   C = [1 -1 0.2] ;
end  
if (nargin < 2)
   B = [1 0.5] ;
end 
if (isempty(B))
   B = [1 0.5] ;
end  
if (nargin < 1)
   A = [1 -1.5 0.7] ;
end
if (isempty(A))
   A = [1 -1.5 0.7] ; 
end
A = roots(A) ; 
A(abs(A)>=1) = 1./A(abs(A)>=1) ; 	% Correct the stability. 
A = poly(A) ; 
% 
% Generating the Gaussian white noise
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
e = lambda*randn(N,1) ; 
% 
% Generating the Gaussian colored noise
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
V = filter(C,1,e) ; 
V = iddata(V,e) ; 	% White noise:   V.u
			% Colored noise: V.y
% 
% Generating the PRB input
% ~~~~~~~~~~~~~~~~~~~~~~~~
u = sigma*sign(randn(N,1)) ; 
% 
% Constructing the data provider model
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
P = idpoly(A,[zeros(1,nk) B],C,1,1,lambda*lambda) ; 
% 
% Generating the data
% ~~~~~~~~~~~~~~~~~~~
D = iddata(sim(P,[u e]),u) ;
%
% END
%

