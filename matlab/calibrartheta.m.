%% === Nested Logit 2x2: solución CERRADA de thetas desde (1)-(5) ===
% - Dentro de cada nido (Peak/Off) NO se divide por mu (sigues tu eq. (2))
% - En el nivel superior SÍ entra mu (tu eq. (4))
% - Identificación: theta^{Peak,Car}_i = 0 para todo i
% - Entradas: tiempos g_{qm} (=T_*_hr), costos, betas por grupo i=1..5, mu,
%             shares observados p_ap, p_bp, p_aop, p_bop (NO-condicionales; suman P_travel_i),
%             y P_travel (si es común por grupoo). Si no das P_travel por grupo,
%             se usa 1 - sum(p_*) de cada grupo.

clear; clc

% ----- Parámetros (tu ejemplo) -----

% Sensibilidades al tiempo de viaje (h), por grupo y modo-periodo:
bet_T_ap  = [-0.6833  -1.1193  -1.6449  -2.3918  -3.8974];  % auto peak
bet_T_aop = [-0.3232 -0.5294 -0.7780 -1.1313 -1.8434];
bet_T_bp  = [-1.2096  -1.9815  -2.9119  -4.2342  -6.8994];  
bet_T_bop = [-0.7941 -1.3008 -1.9116 -2.7797 -4.5294];  % bus off-peak

% Sensibilidad a costo monetario ($), negativa por convención:
bet_cost  = [-1.70, -1.40, -1.20, -1.00, -0.72];
mu        = 0.25;

T_ap_hr  = [0.4932 0.4571 0.4367 0.4345 0.4869];
T_aop_hr = [0.4805 0.4687 0.4172 0.4508 0.4908];
T_bp_hr  = [1.0187 1.0179 1.1108 1.6632 1.8274];
T_bop_hr = [1.0930 1.0966 1.2033 1.6149 1.8492];

C_bp  = (6.75/39) * ones(1,5);
C_bop = C_bp;

L    = 6.399; c0c  = 0.3566;
Cauto = c0c * L;      % $/viaje
C_ap  = Cauto * ones(1,5);
C_aop = C_ap;

% Shares observados NO condicionales (picos + off; no incluyen NT explícito)
p_ap  = [0.1161 0.1365 0.2047 0.2296 0.1620];
p_aop = [0.0932 0.1291 0.1608 0.1663 0.0890];
p_bp  = [0.3292 0.3120 0.2664 0.2796 0.3650];
p_bop = [0.3615 0.3225 0.2681 0.2244 0.2840];

% Si tienes un P_travel común (aquí 0.90), úsalo para P_NT:
Ptravel = 0.90;   % si difiere por grupo, pásalo como vector 1x5
% ----- Pre-chequeos y construcción de P_peak / P_off / P_nt (SIN nargin) -----
p_sum       = p_ap + p_bp + p_aop + p_bop;         % share total "viaja" por grupo
P_peak_obs  = p_ap  + p_bp;
P_off_obs   = p_aop + p_bop;

% Si diste Ptravel (escalar o 1xG), úsalo; si no, infiere P_nt de los shares
if exist('Ptravel','var') && ~isempty(Ptravel)
    if isscalar(Ptravel)
        P_nt_obs = 1 - Ptravel*ones(size(P_peak_obs)); % mismo valor para todos los grupos
    else
        P_nt_obs = 1 - Ptravel(:).';                   % vector fila 1xG
    end
else
    P_nt_obs = max(0, 1 - p_sum);                     % por grupo, a partir de shares incond.
end

% sanity checks
assert(all(P_peak_obs>0 & P_off_obs>0), 'Pico y Off deben ser >0 por grupo.');
assert(all(P_nt_obs>=0),               'Algún P_nt salió negativo: revisa Ptravel/shares.');

% ----- Resolver thetas por grupo (cerrado) -----
G = 1:5;
th_ap  = zeros(1,5);          % identificación: 0
th_bp  = zeros(1,5);
th_aop = zeros(1,5);
th_bop = zeros(1,5);
th_nt  = zeros(1,5);

% ----- Pre-chequeos y construcción de P_peak / P_off / P_nt (SIN nargin) -----
p_sum       = p_ap + p_bp + p_aop + p_bop;      % share total "viaja" por grupo
P_peak_obs  = p_ap  + p_bp;
P_off_obs   = p_aop + p_bop;

% Si diste Ptravel (escalar o 1xG), úsalo; si no, infiere P_nt de los shares
if exist('Ptravel','var') && ~isempty(Ptravel)     % <— scripts: usa exist('var','var')
    if isscalar(Ptravel)
        P_nt_obs = 1 - Ptravel*ones(size(P_peak_obs));   % mismo valor en todos los grupos
    else
        P_nt_obs = 1 - Ptravel(:).';                     % vector fila 1xG
    end
else
    P_nt_obs = max(0, 1 - p_sum);                        % inferido por grupo
end

% sanity checks
assert(all(P_peak_obs>0 & P_off_obs>0), 'Pico y Off deben ser >0 por grupo.');
assert(all(P_nt_obs>=0),               'Algún P_nt salió negativo: revisa Ptravel/shares.');

% Aux para S^{qm}_i = lambda_i*cost + beta_i^{qm}*g
S_ap  = bet_cost.*C_ap  + bet_T_ap .*T_ap_hr;
S_bp  = bet_cost.*C_bp  + bet_T_bp .*T_bp_hr;
S_aop = bet_cost.*C_aop + bet_T_aop.*T_aop_hr;
S_bop = bet_cost.*C_bop + bet_T_bop.*T_bop_hr;

for i = G
    % 1) Diferencia PEAK (dentro-nido) de (2):
    r13 = p_ap(i)/p_bp(i);     % = P(Car|Peak)/P(Bus|Peak)
    Dp  = log(r13) - bet_cost(i)*(C_ap(i)-C_bp(i)) - bet_T_ap(i)*(T_ap_hr(i)-T_bp_hr(i));
    th_ap(i) = 0;                                      % identificación
    th_bp(i) = -Dp;                                    % cerrado

    % 2) A_peak conocido (3):
    A_peak = lse([ th_ap(i)+S_ap(i),  th_bp(i)+S_bp(i) ]);

    % 3) A_off a partir de (4) con ratio P_peak/P_off:
    A_off  = A_peak - (1/mu)*log( P_peak_obs(i) / P_off_obs(i) );

    % 4) Diferencia OFF (2) y nivel absoluto (3) → cerrado
    r24 = p_aop(i)/p_bop(i);
    Doff = log(r24) - bet_cost(i)*(C_aop(i)-C_bop(i)) - bet_T_aop(i)*(T_aop_hr(i)-T_bop_hr(i));

    Loff = lse([ S_aop(i),  S_bop(i) - Doff ]);       % log(e^{S_a}+e^{S_b-D})
    th_aop(i) = A_off - Loff;                         % cerrado
    th_bop(i) = th_aop(i) - Doff;                     % cerrado

    % 5) NT por (4): usando ratio P_nt/P_off (o P_nt/P_peak)
    th_nt(i) = A_off + (1/mu)*log( P_nt_obs(i) / P_off_obs(i) );
end

% ----- Chequeo: reconstruir shares con (2)-(5) -----
% (a) Condicionales dentro de nido:
U_ap  = th_ap  + S_ap;  U_bp  = th_bp  + S_bp;
U_aop = th_aop + S_aop; U_bop = th_bop + S_bop;

denP  = exp(U_ap)  + exp(U_bp);
denOP = exp(U_aop) + exp(U_bop);
PcarP = exp(U_ap)./denP;   PbusP = exp(U_bp)./denP;
PcarO = exp(U_aop)./denOP; PbusO = exp(U_bop)./denOP;

% (b) Niveles superiores:
A_peak = log(denP);
A_off  = log(denOP);
numP   = exp(mu*A_peak);
numO   = exp(mu*A_off);
numNT  = exp(mu*th_nt);
denN   = numP + numO + numNT;
PnestP = numP ./ denN;
PnestO = numO ./ denN;
PnestN = numNT./ denN;

% (c) Shares incondicionales modelados:
p_ap_m  = PnestP.*PcarP;
p_bp_m  = PnestP.*PbusP;
p_aop_m = PnestO.*PcarO;
p_bop_m = PnestO.*PbusO;

% ----- Reporte
fprintf('theta_ap  = [% .4f % .4f % .4f % .4f % .4f]\n', th_ap);
fprintf('theta_bp  = [% .4f % .4f % .4f % .4f % .4f]\n', th_bp);
fprintf('theta_aop = [% .4f % .4f % .4f % .4f % .4f]\n', th_aop);
fprintf('theta_bop = [% .4f % .4f % .4f % .4f % .4f]\n', th_bop);
fprintf('theta_nt  = [% .4f % .4f % .4f % .4f % .4f]\n', th_nt);

gap = [p_ap; p_bp; p_aop; p_bop] - [p_ap_m; p_bp_m; p_aop_m; p_bop_m];
fprintf('max |gap shares| por grupo: '); fprintf('%0.3g ', max(abs(gap),[],1)); fprintf('\n');

% ===== util =====
function y = lse(v)        % log-sum-exp estable
    vm = max(v);
    y  = vm + log(sum(exp(v - vm)));
end
