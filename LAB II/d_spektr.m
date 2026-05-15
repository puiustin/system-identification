function [w,fi]=d_spektr(A,B,sigma2)
%
% D_SPEKTRUM	Evaluates the spectral density of linear systems 
%               (filters) output when stimulated by the white noise. 
%
% Inputs:	A        # transfer function denominator (polynomial)
%        	B        # transfer function numerator (polynomial)
%        	sigma2   # white noise variance
%
% Outputs:	w        # omega (pulsation axis)
%         	fi       # spectral density of filter output
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
lgw1=-2;
w=logspace(lgw1,pi)';
Hp=freqz(B,A,w);
wm=-w;
Hm=freqz(B,A,wm);
fi=sigma2*(Hp.*Hm)/(2*pi);
% 
% END
% 