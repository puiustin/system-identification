function [K, T, Tp] = Hz2Hs_DCeng(M, Ts)
% HZ2HS_DCENG   Converts discrete-time model parameters to physical parameters
%               for a DC engine with a parasitic time constant Tp.
%
%               Continuous model: H(s) = K / [s(1+Ts)(1+Tps)]
%               Discrete model: Hd(z) = (b1*z^-1 + b2*z^-2 + b3*z^-3) /
%                                       (1 + a1*z^-1 + a2*z^-2 + a3*z^-3)
%
% Inputs:       M  # IDMODEL object (ARX or OE)
%               Ts # Sampling period
%
% Outputs:      K  # Estimated gain
%               T  # Estimated main time constant
%               Tp # Estimated parasitic time constant

% Extract coefficients
if isprop(M, 'f') && ~isempty(M.f) % OE model: denominator is 1 + f1*z^-1 + f2*z^-2 + f3*z^-3
    a = M.f ;
else % ARX model: denominator is 1 + a1*z^-1 + a2*z^-2 + a3*z^-3
    a = M.a ;
end
b = M.b ;

% Pad coefficients with zeros in case they are shorter than expected
if length(a) < 4, a = [a, zeros(1, 4 - length(a))]; end
if length(b) < 4, b = [b, zeros(1, 4 - length(b))]; end

a1 = a(2); a2 = a(3); a3 = a(4);
b1 = b(2); b2 = b(3); b3 = b(4);

% Physical Gain K
% Using the property that for H(s) with an integrator, 
% K = lim_{z->1} [(z-1)/Ts * Hd(z)]
% Hd(z) = B(z) / [(1-z^-1) * (1 + (a1+1)z^-1 - a3z^-2)]
K = (b1 + b2 + b3) / (Ts * (2 + a1 - a3));

% Time constants T and Tp
% The non-unity poles of Hd(z) are roots of: u^2 + (1+a1)u - a3 = 0
% These roots correspond to e^(-Ts/T) and e^(-Ts/Tp) in ideal discretization,
% but Problem 12.3 specifies Euler-Pade: z = 1 - Ts/T.
delta = (1 + a1)^2 + 4*a3;
if delta < 0
    warning('Complex poles detected. Physical parameters may be inaccurate.');
    delta = 0;
end

u1 = (-(1 + a1) + sqrt(delta)) / 2;
u2 = (-(1 + a1) - sqrt(delta)) / 2;

% We assume T > Tp, so x = 1 - Ts/T should be larger (closer to 1)
x = max(u1, u2);
y = min(u1, u2);

T = Ts / (1 - x);
Tp = Ts / (1 - y);

end
