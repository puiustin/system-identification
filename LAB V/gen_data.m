function [D,V,P] = gen_data(DP,N,sigma,lambda,bin) 
%
% GENDATA    Module that generates data from a 
%            process model. 
%
% Inputs:       DP     # IDMODEL object representing the 
%                        process that provides the data 
%               N      # simulation period (250, by default)
%               sigma  # standard deviation of PRB input 
%                        (1, by default) 
%                        if null, input inhibited
%               lambda # standard deviation of white noise 
%                        (1, by default)
%                        if null, noise free processes 
%                        are considered
%               bin    # flag indicating the type of input:
%                          0 -> Gaussian PRB
%                          1 -> flip-flop Gaussian PRB (default)
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
% Explanation:  The specified model is stimulated with a PRB signal. 
%               The generated input and the observed output are 
%               returned. Note that various processes can be used, 
%               by setting sigma and/or lambda to null. 
%
% Author:   Dan Stefanoiu (*)
%           Lavinius Ioan Gliga
%
% Created: April 9, 2004
% Revised: August 9, 2018
%
% Copyright: (*) "Politehnica" Unversity of Bucharest, ROMANIA
%                Department of Automatic Control & Computer Science
%

%
% BEGIN
% 
% Faults preventing
% ~~~~~~~~~~~~~~~~~
if (nargin < 5)
   bin = 1 ;
end
if (isempty(bin))
   bin = 1 ;
end 
bin = abs(sign(bin(1))) ; 
if (nargin < 4)
   lambda = 1 ;
end
if (isempty(lambda))
   lambda = 1 ;
end 
lambda = abs(lambda(1)) ; 
if (nargin < 3)
   sigma = 1 ;
end
if (isempty(sigma))
   sigma = 1 ;
end 
sigma = abs(sigma(1)) ; 
if (nargin < 2)
   N = 250 ;
end 
if (isempty(N))
   N = 250 ;
end 
N = abs(fix(N(1))) ; 
if (~N)
   N = 250 ;
end 
if (nargin < 1) 
   DP = idpoly([1 -1.5 0.7],[0 1 0.5],[1 -1 0.2]) ;
end
if (isempty(DP))
   DP = idpoly([1 -1.5 0.7],[0 1 0.5],[1 -1 0.2]) ;
end 
if (~isa(DP,'IDMODEL'))
   DP = idpoly([1 -1.5 0.7],[0 1 0.5],[1 -1 0.2]) ;
end 
P = roots(DP.a) ; 
P(abs(P)>=1) = 1./P(abs(P)>=1) ; 	% Correct the stability. 
DP.a = poly(P) ; 
P = DP ; 
% 
% Generating the Gaussian white noise
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
e = lambda*randn(N,1) ; 
% 
% Generating the Gaussian colored noise
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
V = filter(P.c,1,e) ; 
V = iddata(V,e) ; 	% White noise:   V.u
			% Colored noise: V.y
% 
% Generating the PRB input
% ~~~~~~~~~~~~~~~~~~~~~~~~
u = randn(N,1) ; 
u = sigma*(sign(u).^bin).*(u.^(1-bin)) ; 
% 
% Generating the data
% ~~~~~~~~~~~~~~~~~~~
D = iddata(sim(P,[u e]),u) ;
%
% END
%
