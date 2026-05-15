function [a,l2]=spefac(r)
%
% SPEFAC	Spectrum factorization. 
%
% Inputs:	r   # auto-covariance sequence (vector)
%
% Outputs:	a   # coefficients of spectral factor (AR model - vector)
%         	l2  # variance of AR model
%
% Explanation:	Given the sequence r(0) ... r(n), 
%             	compute the sequence of coefficients a(1) ... a(n) 
%             	and the variance l2 such that: 
%             	l2*[z^n    + a1*z^(n-1) + ... + a(n)]*
%             	   [z^(-n) + a1*z^(1-n) + ... + a(n)] = 
%                = r(0)+r(1)*[z+z^(-1)]+...+r(n)*[z^(n)+z^(-n)]. 
%
% Note: The sequence r is assumed to be positively definite.
%
% Author:   Helena Haglund (*)
% Revised:  Bjorn Wittenmark (*)
%           Dan Stefanoiu (**)
%           Lavinius Ioan Gliga (**)
%
% Last upgrade: (*)  January 3, 1997
%               (**) March   8, 2004
%               (**) February 26, 2018
%
% Copyright: (*)  Lund Institute of Technology, SWEDEN
%                 Department of Automatic Control
%            (**) "Politehnica" Unversity of Bucharest, ROMANIA
%                 Department of Automatic Control & Computer Science
%

%
% BEGIN
%
if r(1)<max(abs(r))
    error('### r is not positive definite') ; 
end
[n,n2]=size(r);
if n<n2
    n=n2;
    r=r';
end
a=r'/sqrt(r(1));
da=1;
k=0;
while da>1e-14
    k=k+1;
    aa = zeros(n, n);
    for i=1:n
        aa(i,:)=[a(i:n) zeros(1,i-1)]+[zeros(1,i-1) a(1:n+1-i)];
    end
    x=2*(aa\r);
    a1=(a+x')/2;
    da=norm(a1-a);
    a=a1;
    if k==50
        error('### Convergence too slow.');
    end
end
l2=a(1)^2; 
a = a/a(1);
% 
% END
%