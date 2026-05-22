function ISLAB_2A(C,A,N,tau_max,nr)
%
% ISLAB_2A   Module that computes and displays true and estimated 
%            covariance functions of an ARMA[1,1] process:
%            Ay=Ce. 
%
% Inputs:   C         # MA polynomial (vector [1 c])
%           A         # AR polynomial (vector [1 a])
%           N         # number of data (100, by default)
%           tau_max   # maximum time shift (50, by default)
%           nr        # number of realizations (1, by default)
%
% Outputs:	---------
%
% Explanation:	The number of data N used in the estimation 
%               and the maximum time shift tau_max can also 
%               be used to study the effect of the individual 
%               realizations by repeating the estimation 
%               procedure nr times.
%
% Author:   Helena Haglund (*)
% Revised:  Bjorn Wittenmark (*)
%           Dan Stefanoiu (**)
%           Lavinius Ioan Gliga (**)
%
% Last upgrade: (*)  January   3, 1997
%               (**) March     8, 2004
%                    February 10, 2012
%                    February 26, 2018
%
% Copyright: (*)  Lund Institute of Technology, SWEDEN
%                 Department of Automatic Control
%            (**) "Politehnica" Unversity of Bucharest, ROMANIA
%                 Department of Automatic Control & Computer Science
%

% setarea valorilor implicite pentru parametri
lam=1;
if nargin<5
   nr=1;
end
if isempty(nr)
   nr=1;
end
if nargin<4
   tau_max=50;
end
if isempty(tau_max)
   tau_max=50;
end
if nargin<3
   N=100;
end
if isempty(N)
   N=100;
end
if nargin<2
   C=1;
end
if isempty(C)
   C=1;
end
if nargin<1
   A=1;
end
if isempty(A)
   A=1;
end
% 
% prevenirea erorilor prin verificarea lungimii polinoamelor
% ~~~~~~~~~~~~~~~~~
if length(A)<2
   a=0; 
else
   a=A(2);
end
if length(C)<2
   c=0;
else
   c=C(2);
end
% 
% calculul raspunsului estimat
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tau=[1:tau_max]';
tau_vec=[0;tau];
e=randn(N,nr); % generare zgomot alb
sys= filt(C,A,1); % definire sistem liniar discret
%r_e = zeros( length(e) - tau_max + 1, nr);
for n=1:nr
  yn=lsim(sys,e(:,n)); % simulare raspuns sistem la zgomot
  [tmp,~]=xcov(yn,tau_max,'unbiased'); % estimare functie de autocovarianta
  r_e(:,n)=tmp(tau_max+1:end);
end
% 
% calculul raspunsului teoretic (valori adevarate)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~
r_t=lam/(1-a^2)*[(1+c^2-2*a*c); ...
                 (c-a)*(1-a*c)*(-a).^(tau-1)];
% 
% afisarea grafica a rezultatelor
% ~~~~~~~~
NN=min(N,50);
subplot(211)
  plot(tau_vec,r_t,'r-',tau_vec,r_e,'b--'); % comparatie teoretic vs estimat
  hold on;
  legend('True','Estimated');
  title(['Covariance functions']);
  xlabel('k');
subplot(212)
  plot(yn(1:NN)); % afisare o realizare a procesului
  grid;
  title(['Realization (',num2str(NN),' samples)']);
  xlabel('Discrete time');
%
% END
%
