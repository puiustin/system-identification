function Mid = bj_e(D, si)
%BJ_E Identificarea modelului Box-Jenkins folosind MCMMPE (via ARMAX)
%   Mid = bj_e(D, si)
%
%   Intrari:
%       D  - Date iddata
%       si - Structura BJ: [nb nc nd nf nk]
%
%   Algoritm:
%       1. Converteste problema BJ intr-una ARMAX extinsa.
%       2. Identifica ARMAX folosind 'armax_e'.
%       3. Extragem polinoamele B, C, D, F prin separarea radacinilor comune.

   
    % si = [nb nc nd nf nk]
    nb = si(1);
    nc = si(2);
    nd = si(3);
    nf = si(4);
    nk = si(5);

  
    % A_tot = F * D  => na_armax = nf + nd
    % B_tot = B * D  => nb_armax = nb + nd
    % C_tot = C * F  => nc_armax = nc + nf
    
    na_x = nf + nd;
    nb_x = nb + nd;
    nc_x = nc + nf;
    
    % Construim structura pentru ARMAX: [na_x nb_x nc_x nk]
    si_armax = [na_x, nb_x, nc_x, nk];
    
    % 3. Identificarea modelului intermediar ARMAX
    % Apeleaza functia armax_e
    M_armax = armax_e(D, si_armax);
    
    % Extragem polinoamele identificate (vectorii de coeficienti)
    At = M_armax.a; %  F * D
    Bt = M_armax.b; %  B * D
    Ct = M_armax.c; %  C * F
    
    % 4. Separarea Polinoamelor (Algoritmul radacinilor comune)
    
    % Gasim radacinile polinoamelor totale
    rA = roots(At);
    rB = roots(Bt);
    
    % --- Identificarea lui D(q) ---
    % D este partea comuna dintre At (F*D) si Bt (B*D).
    % Cautam cele 'nd' radacini din rA care sunt cele mai apropiate de rB.
    
    if nd > 0
        distante = zeros(length(rA), 1);
        % Pentru fiecare radacina din A, gasim distanta pana la cea mai
        % apropiata radacina din B.
        for i = 1:length(rA)
            d_min = min(abs(rA(i) - rB));
           distante(i) = d_min;
        end
        
        % Sortam radacinile lui A in functie de cat de aproape sunt de B
        [~, idx_sort] = sort(distante, 'ascend');
        
        % Cele mai apropiate 'nd' radacini sunt considerate radacinile lui D
        idx_D = idx_sort(1:nd);
        roots_D = rA(idx_D);
        
        % Restul radacinilor din A sunt ale lui F
        idx_F = idx_sort(nd+1:end);
        roots_F = rA(idx_F);
        
        % Reconstruim polinoamele din radacini (si fortam coef real)
        PolyD = real(poly(roots_D));
        PolyF = real(poly(roots_F));
    else
        % Cazul degenerat fara D (nd=0)
        PolyD = 1;
        PolyF = At;
    end
    
    % --- Identificarea lui B(q) ---
    % Avem Bt = B * D. Putem face impartirea polinoamelor (deconvolutie).
    % stim D estimat, aflam B.
    
    [PolyB, R_b] = deconv(Bt, PolyD);
    
    % --- Identificarea lui C(q) ---
    % Avem Ct = C * F. Stim F  aflam C.
    [PolyC, R_c] = deconv(Ct, PolyF);
    
   
    
    Mid = idpoly(1, PolyB, PolyC, PolyD, PolyF, ...
                 'NoiseVariance', M_armax.NoiseVariance, ...
                 'Ts', D.Ts, ...
                 'InputDelay', nk);

end