function [D, V, P] = IOdata_DCeng(type, K0, T0, Tp0, Tmax, Ts, U, lambda)
% IODATA_DCENG  Generates identification data for a DC engine with 
%               a parasitic time constant Tp and various input signals.
%
% Inputs:       type   # type of input signal:
%                        1 -> Square wave (default)
%                        2 -> Pseudo-random pulses (random durations)
%                        3 -> SPAB (PRBS - random +/-U at each Ts)
%               K0     # Gain
%               T0     # Main time constant
%               Tp0    # Parasitic time constant
%               Tmax   # Simulation duration
%               Ts     # Sampling period
%               U      # Input amplitude
%               lambda # White noise standard deviation

% Defaults
if nargin < 8, lambda = 1; end
if nargin < 7, U = 0.5; end
if nargin < 6, Ts = 0.1; end
if nargin < 5, Tmax = 80; end
if nargin < 4, Tp0 = T0/10; end
if nargin < 3, T0 = 0.5; end
if nargin < 2, K0 = 4; end
if nargin < 1, type = 1; end

% Constants
N = round(Tmax/Ts);
Tmax = N*Ts;
N = N + 1;
t = (0:Ts:Tmax)';

% Continuous model: H(s) = K / [s(1+Ts)(1+Tps)]
% Denominator: s(1 + (T+Tp)s + T*Tp*s^2) = s + (T+Tp)s^2 + T*Tp*s^3
P = tf(K0, [T0*Tp0, T0+Tp0, 1, 0]);

% Generate Input Signal
switch type
    case 1 % Standard Square Wave
        nu = 1/7;
        Dhalf = ones(round(N*nu/2), 1);
        Ucycle = U * [Dhalf; -Dhalf];
        u = kron(ones(fix(1/nu)+1, 1), Ucycle);
        u = u(1:N);
        
    case 2 % Pseudo-random pulses
        u = zeros(N, 1);
        curr_val = U;
        k = 1;
        while k <= N
            len = randi([round(2/Ts), round(10/Ts)]); % random duration between 2s and 10s
            end_idx = min(k + len - 1, N);
            u(k:end_idx) = curr_val;
            curr_val = -curr_val;
            k = end_idx + 1;
        end
        
    case 3 % SPAB (Signal Pseudo-Aleator Binar)
        u = U * sign(randn(N, 1));
        
    otherwise
        error('Unknown input type');
end

% Generate Output
y_ideal = lsim(P, u, t);
e = lambda * randn(N, 1);
y_noisy = y_ideal + e;

% Pack Data
D = iddata(y_noisy, u, Ts);
V = iddata(e, [], Ts);
% Store true parameters in P
P_info.K = K0;
P_info.T = T0;
P_info.Tp = Tp0;
P = P_info;

end
