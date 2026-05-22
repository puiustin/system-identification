function [a,b,lambda,magi,phii,mag,phi,f] = ISLAB_4J(na,nb,nk, at,bt,K,N,nr) 
%ARMAX[na,nb] + uf (intrare nefiltrata)
% ISLAB_3D   Module that estimates the parameters of an 
%           OE[2,2] model, with the help of 
%            Least Squares (LS) Method. 
%
% Inputs:	at  # true coefficients of AR part 
%                     ([-0.4 -0.32], by default)
%               bt  # true coefficents of X part 
%                     ([0.5 0.03], by default)
%               K   # number of frequency nodes (50, by default)
%               N   # simulation period (100, by default)
%               nr  # number of realizations (100, by default)
%
% Outputs:	a      # LS estimates of AR coefficients 
%                        (nr-by-2 matrix)
%               b      # LS estimates of X coefficents 
%                        (nr-by-2 matrix)
%               lambda # LS estimates of noise standard deviation 
%                        (nr-length vector)
%               magi   # magnitude of ideal frequency response 
%                        (noise free) (K-length vector)
%               phii   # phase of ideal frequency response 
%                        (noise free) (K-length vector)
%               mag    # magnitude of all nr estimations of 
%                        frequency response (K-by-nr matrix) 
%               phi    # phase of all nr estimations of 
%                        frequency response (K-by-nr matrix) 
%               f      # frequency axis (semi-logarithmic) 
%                        between 10^(-2) and pi (K-length vector) 
%
% Explanation:	The parameters of an ARX[2,2] model are 
%               estimated by means of LS Method. This simulator 
%               shows the estimation errors produced on 
%               Bode diagram of model, as well as the estimation 
%               of noise variance. 
%
% Author:   Dan Stefanoiu (*)
% Revised:  Dan Stefanoiu (*)
%           Lavinius Ioan Gliga (**)
%
% Last upgrade: March 8, 2004 / January 21, 2012
%               April 3, 2018
%
% Copyright: (*) "Politehnica" Unversity of Bucharest, ROMANIA
%                Department of Automatic Control & Computer Science
%

%
% BEGIN
% 

global FIG ;			% Figure number handler 
FIG = 1;    

% Constants
% ~~~~~~~~~
lambdat = 1 ;			% True variance of white noise. 
		 
Ts = 1 ; 			% Sampling period. 
% 
% Faults preventing
% ~~~~~~~~~~~~~~~~~
if (nargin < 8)
   nr = 100 ;
end
if (isempty(nr))
   nr = 100 ;
end 
nr = abs(fix(nr(1))) ; 
if (~nr)
   nr = 100 ;
end 
if (nargin < 7)
   N = 100 ;
end
if (isempty(N))
   N = 100 ;
end  
N = abs(fix(N(1))) ; 
if (~N)
   N = 100 ;
end 
if (nargin < 6)
   K = 50 ;
end 
if (isempty(K))
   K = 50 ;
end 
K = abs(fix(K(1))) ; 
if (~K)
   K = 50 ;
end  
if (nargin < 5)
    nk = 1;
end
if(isempty(nk))
    nk = 1;
end

nk = abs(fix(nk(1)));
if(~nk)
    nk = 1;
end
if (nargin < 4)
   nb = 1;
end 
if (isempty(nb))
   nb = 1;
end  
nb = abs(fix(nb(1))) ; 
if (~nb)
   nb = 1 ;
end  

if (nargin < 3)
   na = 1;
end 
if (isempty(na))
   na = 1;
end  
K = abs(fix(na(1))) ; 
if (~na)
   na = 1 ;
end



if (nargin < 2)
   bt = 0.5;
end 
if (isempty(bt))
   bt = 0.5;
end  
bt = [zeros(1,nk) bt] ;

if (nargin < 1)
   at = -0.4;
end
if (isempty(at))
   at = -0.4  ; 
end
at = [1 at] ; 
a = roots(at) ;
for i = 1 : length(a)
if (abs(a(i))>=1)
   war_err(['<ISLAB_4D>: Model unstable. ' ...
            'Reciprocal model considered instead.']) ; 
   a(i) = 1/a(i) ; 
end 
end

at = poly(a) ; 
% 
% Generating the ideal frequency response
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
f = logspace(-2,pi,K) ; 	% The frequency axis. 
[magi,phii] = dbode(bt,at,Ts,f) ; 
%p = -0.8 ;                                % phii measured in degrees. 
% 
% Generating the Gaussian white noise
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
e = lambdat*randn(N,nr) ;  
% 
% Generating the PRB input
% ~~~~~~~~~~~~~~~~~~~~~~~~
u = sign(randn(N,nr)) ; 
% 
% Filtering and normalizing the input
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%u = filter(1,[1 p],u) ; 
u = u./sqrt(ones(N,1)*sum(u.*u)/N) ; 
% 
% Generating the nr realizations
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
y = filter(bt,at,u)  + e; 
% 
% 
% Estimating the ARX parameters 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
a = zeros(nr, na) ; 
b = zeros(nr, nb) ;
lambda = zeros(nr, 1) ;
mag = zeros(K, nr) ; 
phi = zeros(K, nr) ; 

%am modificat aproximarea zgomotului oentru a putea observa o diferenta
%intre cele doua valori (cea estimata si cea reala)
e = 1.01 *e;
for p = 1:nr
    %pentru modelul oe...vom utiliza y-e in loc de y pentru ca zgomotul
    %afecteaza direct iesirea
    data = iddata(y(:,p) - e(:,p),u(:,p),Ts);

   model = arx(data,[na nb nk]);
   r = [model.A(2:end) model.B((nk+1):end)]';
  
   a_vec = model.A(2:end);
    if length(a_vec) < na
    a_vec = [a_vec, zeros(1, na - length(a_vec))];
    end
    a(p,:) = a_vec;
    b_vec = model.B((nk+1):end);
    if length(b_vec) < nb
    b_vec = [b_vec, zeros(1, nb - length(b_vec))];
    end
    b(p,:) = b_vec;
     %  matricea de regresori pentru eroare 
   e_est = y(:, p);   % prima coloana = y(t)

   %  termenii AR: y(t-1), y(t-2), ..., y(t-na)
   for i = 1:na
       e_est = [e_est, [zeros(i, 1); y(1:(N - i), p) - e(1:(N-i), p)]];
   end

   % termenii X (input): -u(t-nk), -u(t-nk-1), ...
   for i = 1:nb
       shift = i + nk - 1;   % intarziere totală (inclusiv nk)
       e_est = [e_est, -[zeros(shift, 1); u(1:(N - shift), p)]];
   end

   % Vectorul parametrilor: [1; a; b]
   r_vec = [1; a(p, :)'; b(p, :)'];

   % Estimarea erorii 
   e_est = e_est * r_vec;

   % Estimarea deviatiei 
   lambda(p) = norm(e_est) / sqrt(N - na - nb);

% 
% Estimating the frequency response 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   [R,r] = dbode([0 b(p,:)],[1 a(p,:)],Ts,f) ; 
   mag(:, p) = R ; 
   phi(:, p) = r ; 
end 
% 
% Evaluating the average of parameters
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
am = mean(a) ; 
bm = mean(b) ; 
% 
% Evaluating the average of estimation errors
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
magstd = (magi*ones(1,nr)-mag)' ; 
phistd = (phii*ones(1,nr)-phi)' ;
magm = mean(magstd)' ; 
phim = mean(phistd)' ;
% 
% Evaluating the standard deviation
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
magstd = std(magstd,1)' ; 
phistd = std(phistd,1)' ;
% 
% Plotting
% ~~~~~~~~
figure(FIG),clf ; 
   fig_look(FIG,2) ; 
   subplot(211) ; 
      semilogx(f,magm,'-b', ... 
               f,magm-magstd,':r',f,magm+magstd,':r') ; 
      a = axis ; 
      axis([a(1:3) a(4)+0.1*(a(4)-a(3))]) ; 
title(['Estimating an OE[' num2str(na) ',' num2str(nb), '] model ' ... 
          'by the Least Squares Method.']) ; 
      xlabel('Normalized frequency [rad/s] (log)') ; 
      ylabel('FR magnitude') ; 
      set(FIG,'DefaultTextHorizontalAlignment','left') ;
      legend('estimation error', ...
             'standard deviation tube') ; 
   subplot(212) ;
      set(FIG,'DefaultTextHorizontalAlignment','center') ;
      semilogx(f,phim,'-b', ... 
               f,phim-phistd,':r',f,phim+phistd,':r') ; 
      xlabel('Normalized frequency [rad/s] (log)') ; 
      ylabel('FR phase [deg]') ; 
      set(FIG,'DefaultTextHorizontalAlignment','left') ;
      legend('estimation error', ... 
             'standard deviation tube') ; 
FIG = FIG+1 ;
figure(FIG),clf ; 
   fig_look(FIG,2) ; 
   set(FIG,'DefaultTextHorizontalAlignment','center') ;
   plot(1:nr,lambda.^2,'-m',1:nr,ones(nr,1),'-b') ; 
   a = axis ;
   e = a(4)-a(3) ; 
   axis([0 nr+1 a(3)-0.1*3 a(4)+0.25*e]) ; 
   title(['Estimating an OE[' num2str(na) ',' num2str(nb), '] model ' ... 
          'by the Least Squares Method.']) ; 
   xlabel('Realization index') ; 
   ylabel('Noise variance') ; 
   set(FIG,'DefaultTextHorizontalAlignment','left') ; 
   legend('estimated','true') ; 
   text(nr/10,a(4)+0.15*e,... 
        ['        True parameters: ' ... 
          sprintf('%9.4f',[at(2:end) bt((nk+1):end)])]) ; 
   text(nr/10,a(4)+0.05*e,... 
        ['Estimated parameters: ' ... 
          sprintf('%9.4f',[am bm])]) ; 
FIG = FIG+1 ;
%
% END
%