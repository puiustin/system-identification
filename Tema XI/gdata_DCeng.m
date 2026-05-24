function [D,V,P] = gdata_DCeng(cv,K0,T0,Tmax,Ts,U,lambda) 
%
% GDATA_DCENG   Module that generates data from the physical model
%               of a DC engine with constant or variable parameters. 
%
% Inputs:       cv     # flag indicating the type of process: 
%			    cv=0 -> constant parameters (default)
%			    cv=1 -> variable parameters
%               K0     # constant gain (4, by default)
%               T0     # constant time constant (0.5 s, by default)
%               Tmax   # duration of simulation period 
%                        (80 s, by default)
%               Ts     # sampling period 
%                        (0.1 s, by default)
%               U      # Amplitude of square wave input 
%                        (0.5, by default)
%               lambda # standard deviation of white noise 
%                        (1, by default)
%                        if null, a noise free process 
%                        is considered
%
% Outputs:  	D      # IDDATA object representing the 
%                        I/O generated data 
%               V      # IDDATA object representing the 
%                        I/O noise generated data 
%                        (white noise as input, colored noise 
%                        as output)
%               P      # LTI-TF object representing the 
%                        process that provided the data 
%			 (see P.num, P.den for transfer 
%                        function polynomials in case the 
%                        parameters K and T are constant; 
%                        otherwise, P.num{1}=K(t), while 
%                        P.den{1}=T(t))
%
% Explanation:  The continuous model of a DC engine is expressed 
%               by a second order transfer function with a null 
%               pole, gain K and time constant T: 
%                                        K
%                             H(s) = --------- .
%                                     s(1+sT)
%               The model with constant or variable parameters is 
%               employed to generate the identification data. 
%
% Author:   Dan Stefanoiu (*)
% Revised:  Lavinius Ioan Gliga (*)
%
% Created: April 29, 2004 
% Revised: January 30, 2012
%          August 9, 2018
%
% Copyright: (*) "Politehnica" University of Bucharest, ROMANIA
%                Department of Automatic Control & Computer Science

%
% BEGIN
%
% Constants
% ~~~~~~~~~
nu = 1/7 ;			% Fraction representing the ratio between 
                                % the period of input square wave and 
                                % the simulation period. 
a = [3 -3 1] ;                  % Coefficients of parabola that 
                                % defines the variable gain K: 
                                % K(t) = K0*a*
                                %        [t^2 ; t*Tmax ; Tmax^2]/(Tmax^2) .
b = [0.5 10*pi] ;               % Relative amplitude (b(1)) and 
                                % pulsation of harmonic wave (b(2))
                                % that defines the variable 
                                % time constant T: 
                                % T(t) = T0*(1+b(1)*sin(b(2)*t/Tmax)).
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
   U = 0.5 ;
end
if (isempty(U))
   U = 0.5 ;
end 
U = U(1) ; 
if (abs(U)<eps)
   U = 0.5 ;
end 
if (nargin < 5)
   Ts = 0.1 ;
end 
if (isempty(Ts))
   Ts = 0.1 ;
end 
Ts = abs(Ts(1)) ; 
if (Ts<eps)
   Ts = 0.1 ;
end 
if (nargin < 4)
   Tmax = 80 ;
end
if (isempty(Ts))
   Tmax = 80 ;
end 
Tmax = abs(Tmax(1)) ; 
if (Tmax<eps)
   Tmax = 80 ;
end  
Tmax = max(250*Ts,Tmax) ; 
if (nargin < 3) 
   T0 = 0.5 ;
end 
if (isempty(T0))
   T0 = 0.5 ;
end 
T0 = T0(1) ; 
if (abs(T0)<eps)
   T0 = 0.5 ; 
end 
if (nargin < 2) 
   K0 = 4 ;
end
if (isempty(K0))
   K0 = 4 ;
end 
K0 = K0(1) ; 
if (abs(K0)<eps)
   K0 = 4 ; 
end 
if (nargin < 1) 
   cv = 0 ;
end 
if (isempty(cv))
   cv = 0 ;
end 
cv = cv(1) ; 
% 
% Building the continuous model with constant parameters 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
P = tf(K0,[T0 1 0],0) ;
% 
% Generating the Gaussian white noise
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
N = round(Tmax/Ts) ; 		% Set the number of samples. 
Tmax = N*Ts ; 			% Correct Tmax. 
N = N+1 ;                       % Correct number of samples. 
e = lambda*randn(N,1) ; 
% 
% Generating the square wave input
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
D = ones(round(N*nu/2),1) ; 	% Set half period wave. 
U = U*[D ; -D] ; 		% Set input on the whole period. 
u = kron(ones(fix(1/nu),1),U) ; % Set almost all input. 
u = [u ; U] ; 			% Complete the input. 
u = u(1:N) ;                    % Initial null values. 
% 
% Generating the data
% ~~~~~~~~~~~~~~~~~~~
t = 0:Ts:Tmax ; 		% Set time axis. 
if (cv)				% Case: variable parameters. 
   n = t/Tmax ;                 % Set the normalized time. 
                                % Variable gain. 
   K = K0*(a(1)*(n.^2)+a(2)*n+a(3)) ; 
   T = T0*(1+b(1)*sin(b(2)*n)) ;% Variable time constant. 
   u = [0 ; 0 ; u] ;            % Initial null input. 
   D = zeros(N + 2, 1);		% Initial null output. 
   for n=1:N			% Generate the noisy output. 
                                % Continuous time transfer function. 
      V = tf(K(n),[T(n) 1 0],0) ; 
      V = c2d(V,Ts) ;           % Discretizing the transfer function. 
                                % Output ideal data. 
      V = V.num{1}*[0 ; u(n+1) ; u(n)]- ... 
          V.den{1}*[0 ; D(n+1) ; D(n)] ; 
      D(n + 2) = V ; 
   end  
   u(1:2) = [] ;                % Remove initial input.  
   D(1:2) = [] ;                % Remove initial output.
   P.num{1} = K ;               % Record K variation. 
   P.den{1} = T ;               % Record T variation. 
else				% Case: constant parameters. 
   D = lsim(P,u,t) ; 	        % Generate the ideal output. 
end 
D = D+e ;                       % Output noisy data. 
D = iddata(D,u,Ts) ; 	        % Pack input-output data. 
V = iddata(e,e,Ts) ;		% Record the noise. 
 		                % White noise: V.u or V.y.
%
% END
%
