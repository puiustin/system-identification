function Mid = armax_e(D, si)
%ARMAX_E Identificarea modelului ARMAX folosind MCMMPE (MCMMP Extinsa)
%   Mid = armax_e(D, si)
%
%   Intrari:
%       D  - Obiect iddata (date intrare-iesire)
%       si - Vectorul indicilor structurali: [na nb nc nk]
%
%   Iesiri:
%       Mid - Modelul identificat (obiect idpoly)
%
%   Algoritm conform Tema 6:
%   Etapa 1: Estimarea zgomotului folosind un model ARX de ordin mare.
%   Etapa 2: Estimarea parametrilor finali folosind regresori extinsi.

    % 1. Extragerea datelor si a indicilor structurali [cite: 380]
    y = D.y;
    u = D.u;
    N = length(y);
    Ts = D.Ts;
    
    na = si(1);
    nb = si(2);
    nc = si(3);
    nk = si(4);

    %% ETAPA 1: Estimarea zgomotului (ARX de ordin mare) 
    % n_alpha, n_beta >> max(na, nb, nc) (SURSA DE EROARE)
    
    n_big = 3 * max([na, nb, nc]); 
    if n_big < 10, n_big = 10; end % Limita de siguranta
    
    % Construim modelul ARX aproximant: ARX[n_big, n_big, nk]
    % Folosim functia 'arx' din Matlab care aplica MCMMP
    M_arx = arx(D, [n_big n_big nk]); %SURSA EROARE
    
   
    % e_hat[n] = y[n] - y_hat_arx[n] 
    e_data = resid(M_arx, D); %SURSA EROARE
    e_hat = e_data.y;

    %% ETAPA 2: Estimarea finala (Parametrii ARMAX)
    % theta = [a1...ana b1...bnb c1...cnc]
    
    % Indexul de start pentru a avea istoric valid pentru toate variabilele
    start_idx = max([na, nb+nk, nc]) + 1;
    
    Phi = [];
    Y_vect = [];
    
    % phi[n] = [-y[n-1]... | u[n-nk]... | e_hat[n-1]...] 
    for n = start_idx:N
        % Regresorii iesirii (partea A)
        reg_y = -y(n-1 : -1 : n-na)';
        
        % Regresorii intrarii (partea B) 
        if nb > 0
            reg_u = u(n-nk : -1 : n-nk-nb+1)';
        else
            reg_u = [];
        end
        
        % Regresorii zgomotului (partea C) - folosind eroarea din Etapa 1
        %(SURSA DE EROARE) 
        if nc > 0
            reg_e = e_hat(n-1 : -1 : n-nc)';
        else
            reg_e = [];
        end
        
        % Asamblarea liniei curente
        Phi = [Phi; [reg_y, reg_u, reg_e]];
        
        % Vectorul tinta
        Y_vect = [Y_vect; y(n)];
    end
    
    % Rezolvarea sistemului prin CMMP (Theta = Phi \ Y) 
    Theta = Phi \ Y_vect;
    
    %% Construirea modelului final
    % Extragem coeficientii din vectorul Theta
    idx = 0;
    
    % Polinomul A: 1 + a1*q^-1 + ...
    A_poly = [1, Theta(idx+1 : idx+na)'];
    idx = idx + na;
    
    % Polinomul B
    if nb > 0
        B_poly = Theta(idx+1 : idx+nb)';
    else
        B_poly = []; % Explicit gol daca nb=0
    end
    idx = idx + nb;
    
    % Polinomul C: 1 + c1*q^-1 + ...
    C_poly = [1, Theta(idx+1 : idx+nc)'];
    
    % Crearea obiectului IDMODEL (idpoly)
    % Corectie pentru eroare: Nu trimitem argumentele 1, 1 (D si F) explicit
    % si setam InputDelay doar daca avem intrare (B nu e gol).
    
    if isempty(B_poly)
        % Cazul ARMA (fara intrare) - nu setam InputDelay
        Mid = idpoly(A_poly, [], C_poly, ...
                     'NoiseVariance', var(e_hat), ...
                     'Ts', Ts);
    else
        % Cazul ARMAX (cu intrare)
        Mid = idpoly(A_poly, B_poly, C_poly, ...
                     'NoiseVariance', var(e_hat), ...
                     'Ts', Ts, ...
                     'InputDelay', nk);
    end

    % Nota: idpoly seteaza automat D=1 si F=1 daca sunt omise pentru ARMAX.
end