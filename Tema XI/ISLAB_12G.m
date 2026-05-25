function ISLAB_12G(K0, T0, Tp0, Tmax, Ts, U, lambda_vals)
% ISLAB_12G   Mini-simulator for Problem 12.3: Identification of 3 physical 
%             parameters (K, T, Tp) for motor diagnosis.
%
% Inputs:     K0, T0, Tp0, Tmax, Ts, U (standard parameters)
%             lambda_vals # vector of noise intensities to test (SNR analysis)

if nargin < 7, lambda_vals = [0.1 0.5 1.0]; end
if nargin < 6, U = 0.5; end
if nargin < 5, Ts = 0.1; end
if nargin < 4, Tmax = 100; end
if nargin < 3, Tp0 = 0.05; end
if nargin < 2, T0 = 0.5; end
if nargin < 1, K0 = 4; end

global FIG; if isempty(FIG), FIG = 1; end

fprintf('\n========================================================\n');
fprintf('ISLAB_12G: Identificarea a 3 parametri (K, T, Tp) - DIAGNOZA\n');
fprintf('========================================================\n');

for type = 1:3
    switch type
        case 1, type_str = 'Pulsuri Dreptunghiulare';
        case 2, type_str = 'Pulsuri Aleatoare';
        case 3, type_str = 'SPAB (PRBS)';
    end
    
    fprintf('\n>>> Testare semnal intrare: %s\n', type_str);
    
    for lambda = lambda_vals
        % Generate data
        [D, V, P_true] = IOdata_DCeng(type, K0, T0, Tp0, Tmax, Ts, U, lambda);
        
        % Identify OE model (3 poles, 3 zeros, 1 delay)
        % H(s) order 3 -> Hd(z) order 3
        M = oe(D, [3 3 1]);
        
        % Extract physical parameters
        [K_est, T_est, Tp_est] = Hz2Hs_DCeng(M, Ts);
        
        % SNR Calculation
        SNR = var(D.y - V.y) / var(V.y);
        
        fprintf('   lambda=%.2f (SNR=%.2f) | K=%.4f (err=%.1f%%) | T=%.4f (err=%.1f%%) | Tp=%.4f (err=%.1f%%)\n', ...
            lambda, SNR, K_est, abs(K_est-K0)/K0*100, T_est, abs(T_est-T0)/T0*100, Tp_est, abs(Tp_est-Tp0)/Tp0*100);
            
        % Diagnosis logic (Scenario for Requirement d)
        % 5 Wear levels based on Tp (incipient to severe)
        % Wear level 1: Tp < 1.5 * Tp0
        % Wear level 2: Tp < 2.5 * Tp0
        % Wear level 3: Tp < 4.0 * Tp0
        % Wear level 4: Tp < 6.0 * Tp0
        % Wear level 5: Tp >= 6.0 * Tp0 (Severe - Stop Motor!)
        
        ratio = Tp_est / (T0/10); % Relative to "healthy" Tp0 = T0/10
        if ratio < 1.5, level = 1; msg = 'Normal (Incipient)';
        elseif ratio < 2.5, level = 2; msg = 'Uzură Mică';
        elseif ratio < 4.0, level = 3; msg = 'Uzură Moderată';
        elseif ratio < 6.0, level = 4; msg = 'Uzură Avansată';
        else, level = 5; msg = 'CRITIC - REPARAȚIE NECESARĂ!';
        end
        
        % Display Plot only for last lambda in each type
        if lambda == lambda_vals(end)
            figure(FIG); clf;
            subplot(2,1,1);
            plot(D.u, '-b'); hold on; plot(D.y, '-r');
            title(sprintf('I/O Data (%s) - lambda=%.2f', type_str, lambda));
            legend('u', 'y_{meas}');
            
            subplot(2,1,2);
            y_sim = sim(M, D);
            plot(D.y, '-r'); hold on; plot(y_sim.y, '--b');
            title(sprintf('Diagnosis: Level %d - %s (Tp_{est}=%.4f)', level, msg, Tp_est));
            legend('Measured', 'Simulated');
            FIG = FIG + 1;
            drawnow;
        end
    end
end

fprintf('\nConcluzie Diagnoză:\n');
fprintf('Motorul trebuie oprit când Tp depășește pragul critic (ex: Tp > %0.3f s).\n', 6 * (T0/10));
fprintf('Identificarea este mai precisă pentru SNR ridicat (lambda mic) și semnale persistente (SPAB).\n');

end
