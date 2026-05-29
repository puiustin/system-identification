function [D,V,P] = IOdata_DCeng(ip,cv,K0,T0,Tp0,Tmax,Ts,U,lambda)
%
% IODATA_DCENG   Genereaza date I/O pentru un motor de curent continuu
%                cu constanta de timp parazita Tp.
%
% Model continuu:
%
%                  K
% H(s) = ------------------------
%        s*(1+T*s)*(1+Tp*s)
%
% Scop:
%   Functia genereaza date de identificare pentru estimarea parametrilor:
%   K, T si Tp.
%
% Intrari:
%   ip     # tipul semnalului de intrare:
%            0 -> pulsuri dreptunghiulare clasice
%            1 -> pulsuri dreptunghiulare cu durata pseudo-aleatoare
%            2 -> SPAB cu valori -U si +U
%   cv     # tipul parametrilor:
%            0 -> parametri constanti
%            1 -> parametri variabili
%   K0     # amplificarea de baza
%   T0     # constanta principala de timp
%   Tp0    # constanta de timp parazita
%   Tmax   # durata simularii
%   Ts     # perioada de esantionare
%   U      # amplitudinea intrarii
%   lambda # deviatia standard a zgomotului alb
%
% Iesiri:
%   D      # obiect iddata cu datele intrare-iesire
%   V      # obiect iddata cu zgomotul alb generat
%   P      # structura cu modelul real si parametrii reali K, T, Tp
%

%
% BEGIN
%

%
% Constants
% ~~~~~~~~~
nu = 1/7 ;          % raport intre perioada semnalului dreptunghiular si Tmax
a = [3 -3 1] ;      % coeficienti pentru variatia lui K
b = [0.5 10*pi] ;   % amplitudine relativa si pulsatie pentru variatia lui T

%
% Faults preventing
% ~~~~~~~~~~~~~~~~~
if (nargin < 9)
   lambda = 1 ;
end
if (isempty(lambda))
   lambda = 1 ;
end
lambda = abs(lambda(1)) ;

if (nargin < 8)
   U = 0.5 ;
end
if (isempty(U))
   U = 0.5 ;
end
U = U(1) ;
if (abs(U)<eps)
   U = 0.5 ;
end

if (nargin < 7)
   Ts = 0.1 ;
end
if (isempty(Ts))
   Ts = 0.1 ;
end
Ts = abs(Ts(1)) ;
if (Ts<eps)
   Ts = 0.1 ;
end

if (nargin < 6)
   Tmax = 80 ;
end
if (isempty(Tmax))
   Tmax = 80 ;
end
Tmax = abs(Tmax(1)) ;
if (Tmax<eps)
   Tmax = 80 ;
end
Tmax = max(250*Ts,Tmax) ;

if (nargin < 4)
   T0 = 0.5 ;
end
if (isempty(T0))
   T0 = 0.5 ;
end
T0 = abs(T0(1)) ;
if (T0<eps)
   T0 = 0.5 ;
end

if (nargin < 5)
   Tp0 = T0/10 ;
end
if (isempty(Tp0))
   Tp0 = T0/10 ;
end
Tp0 = abs(Tp0(1)) ;
if (Tp0<eps)
   Tp0 = T0/10 ;
end

if (nargin < 3)
   K0 = 4 ;
end
if (isempty(K0))
   K0 = 4 ;
end
K0 = K0(1) ;
if (abs(K0)<eps)
   K0 = 4 ;
end

if (nargin < 2)
   cv = 0 ;
end
if (isempty(cv))
   cv = 0 ;
end
cv = round(abs(cv(1))) ;

if (nargin < 1)
   ip = 0 ;
end
if (isempty(ip))
   ip = 0 ;
end
ip = round(abs(ip(1))) ;

%
% Generarea zgomotului alb
% ~~~~~~~~~~~~~~~~~~~~~~~~
N = round(Tmax/Ts) ;
Tmax = N*Ts ;
N = N+1 ;
e = lambda*randn(N,1) ;

%
% Generarea semnalului de intrare
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
switch ip

   case 0
      % Pulsuri dreptunghiulare clasice.
      Dp = ones(round(N*nu/2),1) ;
      Up = U*[Dp ; -Dp] ;
      u = kron(ones(fix(1/nu),1),Up) ;
      u = [u ; Up] ;
      u = u(1:N) ;

   case 1
      % Pulsuri dreptunghiulare cu durata pseudo-aleatoare.
      u = zeros(N,1) ;
      sgn = 1 ;
      k = 1 ;

      while (k <= N)
         L = randi([round(N*nu/4) round(N*nu)]) ;
         idx = k:min(N,k+L-1) ;
         u(idx) = sgn*U ;
         sgn = -sgn ;
         k = k+L ;
      end

   otherwise
      % SPAB: semnal pseudo-aleator binar cu valori -U si +U.
      u = zeros(N,1) ;
      reg = ones(1,7) ;

      for k = 1:N
         newbit = xor(reg(7),reg(6)) ;
         outbit = reg(end) ;
         reg = [newbit reg(1:end-1)] ;

         if (outbit)
            u(k) = U ;
         else
            u(k) = -U ;
         end
      end

end

%
% Generarea datelor de iesire
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~
t = 0:Ts:Tmax ;

if (cv)

   %
   % Cazul cu parametri variabili.
   % K si T variaza in timp, iar Tp este mentinut constant.
   %

   n = t/Tmax ;

   K  = K0*(a(1)*(n.^2)+a(2)*n+a(3)) ;
   T  = T0*(1+b(1)*sin(b(2)*n)) ;
   Tp = Tp0*ones(size(t)) ;

   u2 = [0 ; 0 ; 0 ; u] ;
   y2 = zeros(N+3,1) ;

   for k = 1:N

      % Se construieste modelul continuu la fiecare esantion,
      % apoi se discretizeaza.
      H = tf(K(k),conv([T(k) 1 0],[Tp(k) 1]),0) ;
      H = c2d(H,Ts) ;

      num = H.num{1}(:).' ;
      den = H.den{1}(:).' ;

      if (length(num) < 4)
         num = [zeros(1,4-length(num)) num] ;
      end

      if (length(den) < 4)
         den = [zeros(1,4-length(den)) den] ;
      end

      % Ecuatia recurenta a modelului discret de ordin 3.
      y2(k+3) = num(1)*u2(k+3) + num(2)*u2(k+2) + ...
                num(3)*u2(k+1) + num(4)*u2(k) - ...
                den(2)*y2(k+2) - den(3)*y2(k+1) - ...
                den(4)*y2(k) ;

   end

   y = y2(4:end) ;

   P = struct ;
   P.K = K(:) ;
   P.T = T(:) ;
   P.Tp = Tp(:) ;
   P.num{1} = K(:)' ;
   P.den{1} = T(:)' ;
   P.par{1} = Tp(:)' ;

else

   %
   % Cazul cu parametri constanti.
   % Se genereaza iesirea folosind modelul continuu cu K0, T0 si Tp0.
   %

   H = tf(K0,conv([T0 1 0],[Tp0 1]),0) ;
   y = lsim(H,u,t) ;

   P = struct ;
   P.H = H ;
   P.K = K0*ones(N,1) ;
   P.T = T0*ones(N,1) ;
   P.Tp = Tp0*ones(N,1) ;
   P.num{1} = K0*ones(1,N) ;
   P.den{1} = T0*ones(1,N) ;
   P.par{1} = Tp0*ones(1,N) ;

end

%
% Adaugarea zgomotului la iesire
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
y = y(:) + e ;

%
% Impachetarea datelor in obiecte iddata
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
D = iddata(y,u,Ts) ;
V = iddata(e,e,Ts) ;

%
% END
%