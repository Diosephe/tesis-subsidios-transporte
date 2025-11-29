% ====================================================
% Tesis
% Yoseph Barrera
% ====================================================
clear; close all; clc

%% Parámetros generales
Y_day   = 9997;                       % Demanda de viajes diaria (total)b
share   = [0.1661 0.1832 0.3047 0.2205 0.1255];  % Proporción por grupo (5 grupos)
Y_grp   = Y_day * share;              % Viajes/día por grupo

L       = 7.548;      % km por viaje
a       = 1.5;        % Ocupación auto [pax/auto]
Q       = 1400;       % Capacidad vial (veq/h)
h_peak  = 8;          % Horas punta
h_off   = 16;         % Horas valle
v_walk  = 4.5;        % km/h (caminar)
c0c     = 0.3566;     % $/km (operación auto)

%% Tiempos de transporte (paradas y abordaje)
tp      = 15/3600;    % h/bus (tiempo fijo por parada)
tsb     = 2.5/3600;   % h/pax (abordaje)
abpr    = 2;          % Parámetro BPR (multiplicador)
bbpr    = 4;          % Parámetro BPR (exponente)

%% Costos de operación de buses
% Costos lineales en el tamaño del bus k:
% Gb(k) = Gb0 + Gb1 * k  → costeo ligado a distancia recorrida (bus·km)
% Gv(k) = Gv0 + Gv1 * k  → costeo ligado a horas de servicio (bus·h)

% Componente por bus·km
Gb0     = 7.1638;      % Intercepto ($/bus·km, ajustado)
Gb1     = 0.0297875;   % Pendiente respecto a k

% Componente por bus·hora
Gv0     = 0.1317;      % Intercepto ($/bus·h, ajustado)
Gv1     = 0.0033;      % Pendiente respecto a k

% Factor de calibración de costos (escala a contabilidad real)
cost_amp = 4.5;

%% Preferencias (Nested Logit)
% Constantes específicas de alternativa (por grupo)
th_ap  = [  0.0000  0.0000  0.0000  0.0000  0.0000 ];
th_aop = [-0.3197 -0.2965 -0.8476 -1.6236 -2.3881];
th_bp  = [-2.1837 -1.4979 -1.1583  1.0270  4.5184];
th_bop = [-2.3512 -2.0009 -2.2554 -2.1158 -0.2420];
th_nt  = [-9.2127 -9.0412 -9.3796 -9.7765 -10.1794]; % no-travel

% Sensibilidades al tiempo de viaje (h), por grupo y modo-periodo:
bet_T_ap  = [-0.6833  -1.1193  -1.6449  -2.3918  -3.8974];  % auto punta
bet_T_aop = [-0.3232  -0.5294  -0.7780  -1.1313  -1.8434];  % auto valle
bet_T_bp  = [-1.2096  -1.9815  -2.9119  -4.2342  -6.8994];  % bus punta
bet_T_bop = [-0.7941  -1.3008  -1.9116  -2.7797  -4.5294];  % bus valle

% Sensibilidad al costo monetario ($)
bet_cost  = [-1.70, -1.40, -1.20, -1.00, -0.72];

% Parámetro de escala del nivel superior del nested logit
mu = 0.25;

%% Finanzas / implementación de políticas
dl_cost = 0;                      % $/km (pistas exclusivas, si aplica)
eta_rev = 0.45;                   % Proporción de ingresos de peajes al privado
mcpf    = 1.8;                   % Costo marginal de fondos públicos

%% Pesos de bienestar
w = ones(1,5);

%% === REF65 — Optimización (caso base con subsidio 65%) ==================

% Subsidio a la operación de buses
sPerc = 65;
s     = sPerc/100;

% Vector de decisión x = [f fop p Pb Pbop YA YB YAop YBop k]
x0 = [ 33, 12, 3.0, 0.6, 0.6, 500, 600, 100, 100, 1 ].';
lb = [  1,  1, 1.0, 0.0, 0.0, 100, 100,  50,  50,   0 ].';
ub = [ 70, 40, 4.0, 2.0, 2.0, 1000, 1000, 650, 650, 160 ].';

% Igualdad lineal: Pb == Pbop
Aeq = zeros(1, numel(x0));
Aeq(4) = 1;
Aeq(5) = -1;
beq = 0;

% Función objetivo: max SW ≡ min(-SW) (incluye dl_cost)
obj = @(x) obj_ref( ...
    x, ...
    L,a,Q,abpr,bbpr,tp,tsb,v_walk, ...
    th_ap,th_bp,th_aop,th_bop,th_nt, ...
    bet_T_ap,bet_T_bp,bet_T_aop,bet_T_bop, ...
    bet_cost, mu, c0c, ...
    Y_grp, h_peak, h_off, ...
    Gb0,Gb1,Gv0,Gv1, cost_amp, eta_rev, mcpf, dl_cost, w);

% Restricciones no lineales: capacidad, flujos y subsidio s
nonlcon = @(x) nonlcon_subshare( ...
    x, s, ...
    L,a,Q,abpr,bbpr,tp,tsb,v_walk, ...
    th_ap,th_bp,th_aop,th_bop,th_nt, ...
    bet_T_ap,bet_T_bp,bet_T_aop,bet_T_bop, bet_cost, mu, c0c, ...
    Y_grp, h_peak, h_off, ...
    Gb0,Gb1,Gv0,Gv1, cost_amp, eta_rev);

% Escala típica (para fmincon)
typX = max([abs(x0), abs(ub), abs(lb), 0.5*(ub-lb), ones(size(x0))],[],2);

% Opciones de fmincon
opts = optimoptions('fmincon', ...
    'Algorithm','sqp', ...
    'Display','iter', ...
    'ScaleProblem',true, ...
    'FiniteDifferenceType','central', ...
    'FiniteDifferenceStepSize',1e-4, ...
    'HessianApproximation','lbfgs', ...
    'TypicalX',[ 33, 12, 3.0, 0.6, 0.6, 500, 600, 200, 200, 45 ], ...
    'MaxFunctionEvaluations', 50000, ...
    'MaxIterations',          5000, ...
    'ConstraintTolerance',    1e-8, ...
    'OptimalityTolerance',    1e-8, ...
    'StepTolerance',          1e-10 );

% Ejecución
[x_opt65, fval65, exitflag65, output65] = fmincon(obj, x0, [],[], Aeq,beq, lb,ub, nonlcon, opts);

%% Reconstrucción de variables óptimas (REF65)
f   = x_opt65(1);  fop = x_opt65(2);  p = x_opt65(3);
Pb  = x_opt65(4);  Pbop= x_opt65(5);
YA  = x_opt65(6);  YB  = x_opt65(7);
YAop= x_opt65(8);  YBop= x_opt65(9);
k   = x_opt65(10);
k1  = (YB*L)/max(f,1e-12);   % Tamaño operativo de bus

% Probabilidades y tiempos
outR = demand_nl(f,p,0,Pb, YA,YB,  fop,p,0,Pbop, YAop,YBop, ...
                 k, L,a,Q,abpr,bbpr,tp,tsb,v_walk, ...
                 th_ap,th_bp,th_aop,th_bop,th_nt, ...
                 bet_T_ap,bet_T_bp,bet_T_aop,bet_T_bop, bet_cost, mu, c0c);

Pshare  = (YA + YB)*h_peak / Y_day;
OPshare = (YAop + YBop)*h_off / Y_day;
NTshare = 1 - Pshare - OPshare;
CPshare = YA/(YA + YB);       BPshare  = YB/(YA + YB);
COPshare= YAop/(YAop+YBop);   BOPshare = YBop/(YAop+YBop);

VelAp  = 1 ./ ta_km(f,  p, 0, Pb,   YA,  YB,  k, L, a, Q, abpr, bbpr, tp, tsb);
VelBp  = 1 ./ tb_km(f,  p, 0, Pb,   YA,  YB,  k, L, a, Q, abpr, bbpr, tp, tsb);
VelAop = 1 ./ ta_km(fop,p, 0, Pbop, YAop,YBop, k, L, a, Q, abpr, bbpr, tp, tsb);
VelBop = 1 ./ tb_km(fop,p, 0, Pbop, YAop,YBop, k, L, a, Q, abpr, bbpr, tp, tsb);

fprintf('\n=== REF65 (s=%d%%) ===\n', sPerc);
fprintf('f=%.4f fop=%.4f p=%.4f Pb=%.4f Pbop=%.4f YA=%.3f YB=%.3f YAop=%.3f YBop=%.3f k=%.3f\n', ...
        f,fop,p,Pb,Pbop,YA,YB,YAop,YBop,k);
fprintf('k1=%.3f | VelAp=%.3f VelBp=%.3f VelAop=%.3f VelBop=%.3f\n', k1, VelAp,VelBp,VelAop,VelBop);
fprintf('Shares: Peak=%.3f Off=%.3f NT=%.3f | Cond: CarP=%.3f BusP=%.3f CarOP=%.3f BusOP=%.3f\n', ...
        Pshare,OPshare,NTshare, CPshare,BPshare, COPshare,BOPshare);

%% Reporte de bienestar (CS y SW) para REF65
SW_opt = -fval65;  % fmincon minimiza -SW

out_opt = demand_nl( ...
    f, p, 0, Pb,    YA,   YB, ...
    fop, p, 0, Pbop,YAop, YBop, ...
    k, L, a, Q, abpr, bbpr, tp, tsb, v_walk, ...
    th_ap, th_bp, th_aop, th_bop, th_nt, ...
    bet_T_ap, bet_T_bp, bet_T_aop, bet_T_bop, bet_cost, mu, c0c);

CS_opt = consumer_surplus(out_opt, Y_grp, bet_cost, mu, w, th_nt);

[RevB_opt, RevC_opt] = revenues(Pb, Pbop, 0, 0, L, a, h_peak, h_off, ...
                                YB, YBop, YA, YAop, eta_rev);

CostB_opt = bus_oper_cost( ...
    f, p, 0, Pb,    YA, YB, ...
    fop, p, 0, Pbop,YAop, YBop, ...
    k, L, a, Q, abpr, bbpr, tp, tsb, ...
    h_peak, h_off, Gb0, Gb1, Gv0, Gv1, cost_amp);

SW_check = social_welfare(CS_opt, RevB_opt, RevC_opt, CostB_opt, mcpf, dl_cost);

% Chequeo de restricciones
[c_opt, ceq_opt] = nonlcon_subshare(x_opt65, s, ...
    L,a,Q,abpr,bbpr,tp,tsb,v_walk, ...
    th_ap,th_bp,th_aop,th_bop,th_nt, ...
    bet_T_ap,bet_T_bp,bet_T_aop,bet_T_bop, bet_cost, mu, c0c, ...
    Y_grp, h_peak, h_off, ...
    Gb0,Gb1,Gv0,Gv1, cost_amp, eta_rev);

viol_ineq = max([0; c_opt]);         % Máxima violación en c(x)≤0
viol_eq   = norm(ceq_opt, Inf);      % Norma infinito de ceq

fprintf('CS=%.2f | SW(from fval)=%.2f | SW(recalc)=%.2f | RevB=%.2f | RevC=%.2f | CostB=%.2f | DL=%.2f\n', ...
        CS_opt, SW_opt, SW_check, RevB_opt, RevC_opt, CostB_opt, dl_cost);
fprintf('Restricciones: max(c)<=0 → %.2e | ||ceq||_inf=%.2e | exitflag=%d\n', ...
        viol_ineq, viol_eq, exitflag65);

% Baseline de bienestar para comparar con otros escenarios
SWref  = SW_opt;
CSref  = CS_opt;
REFsol = struct('f',f,'fop',fop,'p',p,'Pb',Pb,'Pbop',Pbop, ...
                'YA',YA,'YB',YB,'YAop',YAop,'YBop',YBop, ...
                'k',k,'k1',k1,'sPerc',sPerc);

save('baseline_REF.mat','SWref','CSref','REFsol');
save('baseline_SUB65.mat','SWref','CSref','REFsol');

%% Elasticidades REF (oferta fija)
eps_guard = @(x) max(x,1e-12);

% Shares por grupo desde outR
S_P  = outR.P.nest.P;       % S_{P,i}
S_O  = outR.P.nest.O;       % S_{O,i}
s_bP = outR.P.cond.P.bus;   % s_{b|P,i}
s_bO = outR.P.cond.O.bus;   % s_{b|O,i}

% d pi_{q,i} / d P_b^{(q)}
dPi_dPb_P = S_P .* s_bP .* bet_cost .* ( (1 - s_bP) + mu .* (1 - S_P) .* s_bP );   % 1x5
dPi_dPb_O = S_O .* s_bO .* bet_cost .* ( (1 - s_bO) + mu .* (1 - S_O) .* s_bO );   % 1x5

% d YB^{(q)} / d P_b^{(q)}  (viajes/h)
dYB_dPb   = sum(Y_grp .* dPi_dPb_P) / h_peak;
dYBop_dPb = sum(Y_grp .* dPi_dPb_O) / h_off;

% Elasticidad propia agregada por franja
eps_peak = (Pb   / eps_guard(YB  )) * dYB_dPb;
eps_off  = (Pbop / eps_guard(YBop)) * dYBop_dPb;

% Elasticidad propia total diaria (si Pb = Pbop)
Ybus_day   = YB*h_peak + YBop*h_off;
dYday_dPb  = sum(Y_grp .* (dPi_dPb_P + dPi_dPb_O));
eps_total  = (Pb / eps_guard(Ybus_day)) * dYday_dPb;

% Elasticidades propias por grupo (punta/valle)
eps_peak_grp = Pb   .* bet_cost .* ( (1 - s_bP) + mu .* (1 - S_P) .* s_bP );  % 1x5
eps_off_grp  = Pbop .* bet_cost .* ( (1 - s_bO) + mu .* (1 - S_O) .* s_bO );  % 1x5

% Elasticidad diaria por grupo (Pb=Pbop)
brP = (1 - s_bP) + mu .* (1 - S_P) .* s_bP;
brO = (1 - s_bO) + mu .* (1 - S_O) .* s_bO;
den = eps_guard(S_P .* s_bP + S_O .* s_bO);
eps_daily_grp = Pb .* bet_cost .* ( (S_P .* s_bP .* brP + S_O .* s_bO .* brO) ./ den );  % 1x5

fprintf('\n== Elasticidades REF (oferta fija) ==\n');
fprintf('Agregadas:  eps_peak = %.3f | eps_off = %.3f | eps_total = %.3f\n', ...
        eps_peak, eps_off, eps_total);

grp = (1:numel(Y_grp)).';
T_eps_grp = table(grp, eps_peak_grp.', eps_off_grp.', eps_daily_grp.', ...
    'VariableNames', {'Grupo','Eps_Peak','Eps_OffPeak','Eps_Daily'});
disp(T_eps_grp);

% Chequeo numérico local (sin reoptimizar oferta)
bump = 0.01;  % +1 %
% Punta
out1P = demand_nl(f,p,0,(1+bump)*Pb, YA,YB,  fop,p,0,Pbop, YAop,YBop, ...
                  k, L,a,Q,abpr,bbpr,tp,tsb,v_walk, ...
                  th_ap,th_bp,th_aop,th_bop,th_nt, ...
                  bet_T_ap,bet_T_bp,bet_T_aop,bet_T_bop, bet_cost, mu, c0c);
YB0P_h = sum(Y_grp .* outR.P.nest.P .* outR.P.cond.P.bus) / h_peak;
YB1P_h = sum(Y_grp .* out1P.P.nest.P .* out1P.P.cond.P.bus) / h_peak;
eps_peak_num = ((YB1P_h - YB0P_h)/eps_guard(YB0P_h)) / bump;

% Valle
out1O = demand_nl(f,p,0,Pb, YA,YB,  fop,p,0,(1+bump)*Pbop, YAop,YBop, ...
                  k, L,a,Q,abpr,bbpr,tp,tsb,v_walk, ...
                  th_ap,th_bp,th_aop,th_bop,th_nt, ...
                  bet_T_ap,bet_T_bp,bet_T_aop,bet_T_bop, bet_cost, mu, c0c);
YB0O_h = sum(Y_grp .* outR.P.nest.O .* outR.P.cond.O.bus) / h_off;
YB1O_h = sum(Y_grp .* out1O.P.nest.O .* out1O.P.cond.O.bus) / h_off;
eps_off_num = ((YB1O_h - YB0O_h)/eps_guard(YB0O_h)) / bump;

fprintf('Chequeos numéricos (+1%%): peak = %.3f | off = %.3f\n', eps_peak_num, eps_off_num);

%% Elasticidades cruzadas auto→bus (semi-elasticidades; Pa=0 en REF)
WITH_CROSS = true;
if WITH_CROSS
    s_aP = 1 - s_bP;
    s_aO = 1 - s_bO;
    dPi_dPa_P = S_P .* s_bP .* bet_cost .* (L/a) .* ( -s_aP + mu .* (1 - S_P) .* s_aP );
    dPi_dPa_O = S_O .* s_bO .* bet_cost .* (L/a) .* ( -s_aO + mu .* (1 - S_O) .* s_aO );

    dYB_dPa   = sum(Y_grp .* dPi_dPa_P) / h_peak;
    dYBop_dPa = sum(Y_grp .* dPi_dPa_O) / h_off;

    semi_peak_car = (1 / eps_guard(YB  )) * dYB_dPa;
    semi_off_car  = (1 / eps_guard(YBop)) * dYBop_dPa;

    semi_peak_car_grp = bet_cost .* (L/a) .* ( -s_aP + mu .* (1 - S_P) .* s_aP );
    semi_off_car_grp  = bet_cost .* (L/a) .* ( -s_aO + mu .* (1 - S_O) .* s_aO );

    fprintf('Semi-elast. cruzadas: dln YB^P/dPa = %.4f | dln YB^O/dPa = %.4f\n', ...
            semi_peak_car, semi_off_car);

    T_cross_grp = table(grp, semi_peak_car_grp.', semi_off_car_grp.', ...
        'VariableNames', {'Grupo','Semi_dlnYB_P_dPa','Semi_dlnYB_O_dPa'});
    disp(T_cross_grp);
end

%% Funciones locales para el barrido de subsidio (SUBX)
function val = obj_noX(x, L,a,Q,abpr,bbpr,tp,tsb,v_walk, ...
    th_ap,th_bp,th_aop,th_bop,th_nt, ...
    bet_T_ap,bet_T_bp,bet_T_aop,bet_T_bop, bet_cost, mu, c0c, ...
    Y_grp, h_peak, h_off, ...
    Gb0,Gb1,Gv0,Gv1, cost_amp, eta_rev, mcpf, dl_cost, w)

    f=x(1); fop=x(2); p=x(3); Pb=x(4); Pbop=x(5);
    YA=x(6); YB=x(7); YAop=x(8); YBop=x(9); k=x(10);

    out = demand_nl(f,p,0,Pb, YA,YB,  fop,p,0,Pbop, YAop,YBop, ...
                    k, L,a,Q,abpr,bbpr,tp,tsb,v_walk, ...
                    th_ap,th_bp,th_aop,th_bop,th_nt, ...
                    bet_T_ap,bet_T_bp,bet_T_aop,bet_T_bop, bet_cost, mu, c0c);
    [RevB, RevC] = revenues(Pb, Pbop, 0, 0, L,a,h_peak,h_off, YB,YBop, YA,YAop, eta_rev);
    CostB = bus_oper_cost( ...
                f,p,0,Pb, YA,YB,  fop,p,0,Pbop, YAop,YBop, ...
                k, L,a,Q,abpr,bbpr,tp,tsb, h_peak,h_off, Gb0,Gb1,Gv0,Gv1, cost_amp);
    CS  = consumer_surplus(out, Y_grp, bet_cost, mu, w, th_nt);
    SW  = social_welfare(CS, RevB, RevC, CostB, mcpf, dl_cost);
    val = -SW;
end

function [c,ceq] = nonlcon_subshare(x, s, ...
    L,a,Q,abpr,bbpr,tp,tsb,v_walk, ...
    th_ap,th_bp,th_aop,th_bop,th_nt, ...
    bet_T_ap,bet_T_bp,bet_T_aop,bet_T_bop, bet_cost, mu, c0c, ...
    Y_grp, h_peak, h_off, ...
    Gb0,Gb1,Gv0,Gv1, cost_amp, eta_rev)

    f=x(1); fop=x(2); p=x(3); Pb=x(4); Pbop=x(5);
    YA=x(6); YB=x(7); YAop=x(8); YBop=x(9); k=x(10);

    % Capacidad por bus (ambos periodos)
    c1 = (YB  * L)/max(f,   1e-12) - k;
    c2 = (YBop* L)/max(fop, 1e-12) - k;

    % Capacidad de paradero
    capP = QP(f,   p, 0, Pb,   YA,   YB);
    capO = QP(fop, p, 0, Pbop, YAop, YBop);
    c3 = f   - capP;
    c4 = fop - capO;
    c  = [c1;c2;c3;c4];

    % Flujos consistentes con el nested logit
    out = demand_nl(f,p,0,Pb, YA,YB,  fop,p,0,Pbop, YAop,YBop, ...
                    k, L,a,Q,abpr,bbpr,tp,tsb,v_walk, ...
                    th_ap,th_bp,th_aop,th_bop,th_nt, ...
                    bet_T_ap,bet_T_bp,bet_T_aop,bet_T_bop, bet_cost, mu, c0c);
    e1 = YA   * h_peak - sum(Y_grp .* out.P.nest.P .* out.P.cond.P.car);
    e2 = YB   * h_peak - sum(Y_grp .* out.P.nest.P .* out.P.cond.P.bus);
    e3 = YAop * h_off  - sum(Y_grp .* out.P.nest.O .* out.P.cond.O.car);
    e4 = YBop * h_off  - sum(Y_grp .* out.P.nest.O .* out.P.cond.O.bus);

    % Balance con subsidio: RevB = (1-s)*CostB
    [RevB, ~] = revenues(Pb, Pbop, 0, 0, L,a,h_peak,h_off, YB,YBop, YA,YAop, eta_rev);
    CostB = bus_oper_cost( ...
                f,p,0,Pb, YA,YB,  fop,p,0,Pbop, YAop,YBop, ...
                k, L,a,Q,abpr,bbpr,tp,tsb, h_peak,h_off, Gb0,Gb1,Gv0,Gv1, cost_amp);
    e5 = RevB - (1 - s) * CostB;

    ceq = [e1; e2; e3; e4; e5];
end

%% SUBSIDY SWEEP (SUBX): barrido de subsidio s = 0:5:100
% RevB = (1-s)*CostB. Para s=100% se fija Pb=Pbop=0.

submax = 100; step = 5;
grilla   = 0:step:submax;                 % [0 5 10 ... 100]
limit  = numel(grilla);

% Preasignación
f1 = nan(limit,1); f1op = nan(limit,1); p1 = nan(limit,1);
Pb1 = nan(limit,1); Pb1op = nan(limit,1);
YA1 = nan(limit,1); YB1  = nan(limit,1);
YA1op = nan(limit,1); YB1op = nan(limit,1);
bussize = nan(limit,1);
Pshare = nan(limit,1); OPshare = nan(limit,1); NTshare = nan(limit,1);
CPshare = nan(limit,1); BPshare = nan(limit,1);
COPshare = nan(limit,1); BOPshare = nan(limit,1);
VelAp = nan(limit,1); VelBp = nan(limit,1); VelAop = nan(limit,1); VelBop = nan(limit,1);
SW1 = nan(limit,1); CS1 = nan(limit,1);

% Semilla inicial; luego se usa la solución anterior
seed = struct('f',33,'fop',12,'p',3.0,'Pb',0.6,'Pbop',0.6, ...
              'YA',500,'YB',600,'YAop',200,'YBop',200,'k',120);

for t = 1:limit
    sPerc = grilla(t);
    s     = sPerc/100;

    % x0: primera iteración usa seed; luego last solution
    if t==1
        x0 = [seed.f; seed.fop; seed.p; seed.Pb; seed.Pbop; ...
              seed.YA; seed.YB; seed.YAop; seed.YBop; seed.k];
    else
        x0 = [f1(t-1); f1op(t-1); p1(t-1); Pb1(t-1); Pb1op(t-1); ...
              YA1(t-1); YB1(t-1); YA1op(t-1); YB1op(t-1); bussize(t-1)];
    end

    % Cotas
    lb = [  1,  1, 1.0, 0.0, 0.0, 50, 50,  20,  20,   1 ].';
    ub = [ 70, 50, 4.0, 10.0, 10.0, 1000, 1000, 750, 750, 160 ].';

    % s=100% → Pb=Pbop=0
    if s==1
        lb(4)=0; ub(4)=0;
        lb(5)=0; ub(5)=0;
    end

    % Igualdad lineal: Pb == Pbop
    Aeq = zeros(1,10); Aeq(4)=1; Aeq(5)=-1; beq = 0;

    % Objetivo
    obj = @(x) obj_noX(x, L,a,Q,abpr,bbpr,tp,tsb,v_walk, ...
                         th_ap,th_bp,th_aop,th_bop,th_nt, ...
                         bet_T_ap,bet_T_bp,bet_T_aop,bet_T_bop, bet_cost, mu, c0c, ...
                         Y_grp, h_peak, h_off, ...
                         Gb0,Gb1,Gv0,Gv1, cost_amp, eta_rev, mcpf, dl_cost, w);

    % Restricciones no lineales con subsidio s
    nonlcon = @(x) nonlcon_subshare(x, s, ...
                         L,a,Q,abpr,bbpr,tp,tsb,v_walk, ...
                         th_ap,th_bp,th_aop,th_bop,th_nt, ...
                         bet_T_ap,bet_T_bp,bet_T_aop,bet_T_bop, bet_cost, mu, c0c, ...
                         Y_grp, h_peak, h_off, ...
                         Gb0,Gb1,Gv0,Gv1, cost_amp, eta_rev);

    opts = optimoptions('fmincon', ...
        'Algorithm','sqp', ...
        'Display','iter', ...
        'ScaleProblem',true, ...
        'FiniteDifferenceType','central', ...
        'FiniteDifferenceStepSize',1e-4, ...
        'HessianApproximation','lbfgs', ...
        'TypicalX',[ 33, 12, 3.0, 0.6, 0.6, 500, 600, 200, 200, 45 ], ...
        'MaxFunctionEvaluations', 50000, ...
        'MaxIterations',          50000, ...
        'ConstraintTolerance',    1e-8, ...
        'OptimalityTolerance',    1e-8, ...
        'StepTolerance',          1e-10 );

    [x, fval, exitflag, out] = fmincon(obj, x0, [],[], Aeq,beq, lb,ub, nonlcon, opts);

    % Desempaquetar solución
    f     = x(1);   fop  = x(2);   p    = x(3);
    Pb    = x(4);   Pbop = x(5);
    YA    = x(6);   YB   = x(7);
    YAop  = x(8);   YBop = x(9);
    k     = x(10);

    f1(t)=f; f1op(t)=fop; p1(t)=p;
    Pb1(t)=Pb; Pb1op(t)=Pbop;
    YA1(t)=YA; YB1(t)=YB; YA1op(t)=YAop; YB1op(t)=YBop;
    bussize(t)=k;

    % Cálculos asociados a la solución (SW, CS, shares, velocidades)
    out_d = demand_nl(f,p,0,Pb, YA,YB,  fop,p,0,Pbop, YAop,YBop, ...
                      k, L,a,Q,abpr,bbpr,tp,tsb,v_walk, ...
                      th_ap,th_bp,th_aop,th_bop,th_nt, ...
                      bet_T_ap,bet_T_bp,bet_T_aop,bet_T_bop, bet_cost, mu, c0c);
    [RevB, RevC] = revenues(Pb, Pbop, 0, 0, L,a,h_peak,h_off, YB,YBop, YA,YAop, eta_rev);
    CostB = bus_oper_cost( ...
                f,p,0,Pb, YA,YB,  fop,p,0,Pbop, YAop,YBop, ...
                k, L,a,Q,abpr,bbpr,tp,tsb, h_peak,h_off, Gb0,Gb1,Gv0,Gv1, cost_amp);
    CS  = consumer_surplus(out_d, Y_grp, bet_cost, mu, w, th_nt);
    SW  = social_welfare(CS, RevB, RevC, CostB, mcpf);

    SW1(t) = SW;   CS1(t) = CS;

    % Shares y velocidades
    Pshare(t)  = (YA + YB)*h_peak / Y_day;
    OPshare(t) = (YAop + YBop)*h_off / Y_day;
    NTshare(t) = 1 - Pshare(t) - OPshare(t);
    CPshare(t) = YA/(YA + YB);    BPshare(t)  = YB/(YA + YB);
    COPshare(t)= YAop/(YAop+YBop); BOPshare(t)= YBop/(YAop+YBop);

    VelAp(t)  = 1/ta_km(f,  p,0,Pb,   YA,  YB,  k, L,a,Q,abpr,bbpr,tp,tsb);
    VelBp(t)  = 1/tb_km(f,  p,0,Pb,   YA,  YB,  k, L,a,Q,abpr,bbpr,tp,tsb);
    VelAop(t) = 1/ta_km(fop,p,0,Pbop, YAop,YBop,k, L,a,Q,abpr,bbpr,tp,tsb);
    VelBop(t) = 1/tb_km(fop,p,0,Pbop, YAop,YBop,k, L,a,Q,abpr,bbpr,tp,tsb);

    % Resumen por nivel de subsidio
    Pb_disp   = max(Pb,   0);
    Pbop_disp = max(Pbop, 0);
    fprintf('{sub=%3.0f%%, SW=%.1f, {f -> %.4f, fop -> %.4f, p -> %.5f, Pb -> %.6f, Pbop -> %.6f, YA -> %.3f, YB -> %.3f, YAop -> %.3f, YBop -> %.3f, k -> %.3f}}\n', ...
            sPerc, SW, f, fop, p, Pb_disp, Pbop_disp, YA, YB, YAop, YBop, k);
end

%% Gráfico ΔSW y tabla de resultados (POST SUBX)
% Requiere que el bucle SUBX haya llenado SW1, CS1, etc.

if exist('subs','var') && ~isempty(subs)
    subs = subs(:);
elseif exist('grilla','var') && ~isempty(grilla)
    subs = grilla(:);
else
    assert(exist('SW1','var')==1,'Falta SW1 (resultado del SUBX).');
    N = numel(SW1);
    subs = linspace(0,100,N).';
end

% Asegurar vectores columna
SW1     = SW1(:);      CS1    = CS1(:);
f1      = f1(:);       f1op   = f1op(:);    p1 = p1(:);
Pb1     = Pb1(:);      Pb1op  = Pb1op(:);
YA1     = YA1(:);      YB1    = YB1(:);
YA1op   = YA1op(:);    YB1op  = YB1op(:);
bussize = bussize(:);

% Variables opcionales
optNames = {'VelAp','VelBp','VelAop','VelBop', ...
            'Pshare','OPshare','NTshare', ...
            'CPshare','BPshare','COPshare','BOPshare'};
for v = optNames
    vn = v{1};
    if exist(vn,'var') && ~isempty(eval(vn))
        assignin('caller', vn, eval([vn '(:)']));
    end
end

% Base para ΔSW (65% o 0%)
USE_BASE_65 = true;
idx65 = find(abs(subs-65)<1e-12, 1);
idx0  = find(abs(subs-0 )<1e-12, 1);
assert(~isempty(idx65) && ~isempty(idx0), 'La grilla debe contener 0%% y 65%%.');

idxBase = idx65;
if ~USE_BASE_65
    idxBase = idx0;
end
SWbase  = SW1(idxBase);
dSW     = SW1 - SWbase;

tol = 1e-6;
fprintf('Check ΔSW(base=%d%%): dSW(base)=%.3e (tol=%g)\n', subs(idxBase), dSW(idxBase), tol);

% Gráfico ΔSW(s)
figure('Color','w');
plot(subs, dSW, '-o', 'LineWidth',1.5, 'MarkerFaceColor',[0.10 0.40 0.80]);
grid on; box on;
xlabel('Subsidio (%)');
ylabel(sprintf('\\Delta SW = SW(s) - SW(%d%%)', subs(idxBase)));
title('Beneficio del subsidio (SUBX)');
yline(0,'--','Color',[0.3 0.3 0.3]);

% Tabla resumen SUBX
if exist('L','var') && ~isempty(L)
    BusFarePeak_perkm  = Pb1  ./ L;
    BusFareOff_perkm   = Pb1op./ L;
else
    BusFarePeak_perkm  = NaN(size(Pb1));
    BusFareOff_perkm   = NaN(size(Pb1op));
end
CarTollPeak = zeros(size(subs));
CarTollOff  = zeros(size(subs));

T = table( ...
    subs, SW1, CS1, dSW, ...
    f1, f1op, p1, ...
    Pb1, Pb1op, BusFarePeak_perkm, BusFareOff_perkm, ...
    YA1, YB1, YA1op, YB1op, bussize, ...
    'VariableNames', { ...
      'Sub','SW','CS','dSW', ...
      'f','fop','p', ...
      'Pb','Pbop','BusFarePeak_perkm','BusFareOff_perkm', ...
      'YA','YB','YAop','YBop','BusSize'});
disp(T)

% Añadir columnas opcionales
if exist('VelAp','var'),   T.VelAp   = VelAp;   end
if exist('VelBp','var'),   T.VelBp   = VelBp;   end
if exist('VelAop','var'),  T.VelAop  = VelAop;  end
if exist('VelBop','var'),  T.VelBop  = VelBop;  end
if exist('Pshare','var'),  T.PeakShare = Pshare;   end
if exist('OPshare','var'), T.OffShare  = OPshare;  end
if exist('NTshare','var'), T.NTshare   = NTshare;  end
if exist('CPshare','var'), T.CarShare_P = CPshare; end
if exist('BPshare','var'), T.BusShare_P = BPshare; end
if exist('COPshare','var'),T.CarShare_O = COPshare;end
if exist('BOPshare','var'),T.BusShare_O = BOPshare;end
T.CarTollPeak = CarTollPeak;
T.CarTollOff  = CarTollOff;

% Reordenar columnas (si existen)
ord = {'Sub','SW','dSW','CS', ...
       'Pb','Pbop','BusFarePeak_perkm','BusFareOff_perkm', ...
       'CarTollPeak','CarTollOff', ...
       'f','fop','BusSize','p', ...
       'VelAp','VelBp','VelAop','VelBop', ...
       'PeakShare','OffShare','NTshare', ...
       'CarShare_P','BusShare_P','CarShare_O','BusShare_O', ...
       'YA','YB','YAop','YBop'};
for k = 1:numel(ord)
    v = ord{k};
    if ismember(v, T.Properties.VariableNames)
        T = movevars(T, v, 'Before', 1);
    end
end

disp(T);

%% DEMANDA PURA — tarifas diferenciadas (grupo×franja, autofinanciamiento)
% Vector de decisión x_pure (18 variables):
%   [f fop p YA YB YAop YBop k PbP_1..PbP_5 PbO_1..PbO_5]

% Semilla coherente con SUBX
x0_pure          = zeros(18,1);
x0_pure(1)       = 33;     % f
x0_pure(2)       = 12;     % fop
x0_pure(3)       = 3.0;    % p
x0_pure(4)       = 500;    % YA
x0_pure(5)       = 600;    % YB
x0_pure(6)       = 200;    % YAop
x0_pure(7)       = 200;    % YBop
x0_pure(8)       = 120;    % k

PbP0             = 0.6;    % tarifa bus peak (semilla)
PbO0             = 0.6;    % tarifa bus off (semilla)
x0_pure(9:13)    = PbP0;
x0_pure(14:18)   = PbO0;

% Cotas coherentes con REF/SUBX
lb_pure = [ ...
      1;   ... f    [bus/h]
      1;   ... fop  [bus/h]
    1.0;  ... p     [stops/km]
    60;   ... YA    [pax/h]
    80;   ... YB
     30;  ... YAop
     30;  ... YBop
      1;  ... k
    zeros(10,1) ...
    ];

ub_pure = [ ...
     70;   ... f
     40;   ... fop
    4.0;   ... p
    900;   ... YA
    900;   ... YB
    600;   ... YAop
    600;   ... YBop
    160;   ... k
    20*ones(10,1) ...
    ];

% Escala típica
typX_pure = max([abs(x0_pure), abs(ub_pure), abs(lb_pure), ...
                 0.5*(ub_pure-lb_pure), ones(size(x0_pure))],[],2);

% Opciones de fmincon
opts_pure = optimoptions('fmincon', ...
    'Algorithm','sqp', ...
    'Display','iter', ...
    'ScaleProblem',true, ...
    'FiniteDifferenceType','central', ...
    'FiniteDifferenceStepSize',1e-4, ...
    'HessianApproximation','lbfgs', ...
    'TypicalX',typX_pure, ...
    'MaxFunctionEvaluations',80000, ...
    'MaxIterations',8000, ...
    'ConstraintTolerance',1e-6, ...
    'OptimalityTolerance',1e-6, ...
    'StepTolerance',1e-10);

% Objetivo DEMANDA PURA: max CS ≡ min(-SW) (sin peajes ni DLcost)
objDP = @(x) obj_pure_demanda( ...
    x, ...
    L,a,Q,abpr,bbpr,tp,tsb,v_walk, ...
    th_ap,th_bp,th_aop,th_bop,th_nt, ...
    bet_T_ap,bet_T_bp,bet_T_aop,bet_T_bop, ...
    bet_cost, mu, c0c, ...
    Y_grp, h_peak, h_off, ...
    Gb0,Gb1,Gv0,Gv1, cost_amp, mcpf, w);

% Restricciones no lineales DEMANDA PURA
nonlconDP = @(x) nonlcon_pure_demanda( ...
    x, ...
    L,a,Q,abpr,bbpr,tp,tsb,v_walk, ...
    th_ap,th_bp,th_aop,th_bop,th_nt, ...
    bet_T_ap,bet_T_bp,bet_T_aop,bet_T_bop, ...
    bet_cost, mu, c0c, ...
    Y_grp, h_peak, h_off, ...
    Gb0,Gb1,Gv0,Gv1, cost_amp);

% Optimización de DEMANDA PURA
[x_pure, fval_pure, exitflag_pure, output_pure] = ...
    fmincon(objDP, x0_pure, [],[], [],[], lb_pure, ub_pure, nonlconDP, opts_pure);

% Reconstruir solución
fDP    = x_pure(1);
fopDP  = x_pure(2);
pDP    = x_pure(3);
YADP   = x_pure(4);
YBDP   = x_pure(5);
YAopDP = x_pure(6);
YBopDP = x_pure(7);
kDP    = x_pure(8);
PbP_DP = (x_pure(9:13)).';
PbO_DP = (x_pure(14:18)).';

% Demanda asociada a DEMANDA PURA
outDP = demand_nl( ...
    fDP,  pDP, 0, PbP_DP,  YADP,   YBDP, ...
    fopDP,pDP, 0, PbO_DP,  YAopDP, YBopDP, ...
    kDP, L,a,Q,abpr,bbpr,tp,tsb,v_walk, ...
    th_ap,th_bp,th_aop,th_bop,th_nt, ...
    bet_T_ap,bet_T_bp,bet_T_aop,bet_T_bop, bet_cost, mu, c0c);

% Flujos diarios de bus por grupo (para RevB)
yPb_by_grp = Y_grp .* outDP.P.nest.P .* outDP.P.cond.P.bus;
yOb_by_grp = Y_grp .* outDP.P.nest.O .* outDP.P.cond.O.bus;

RevB_DP = sum(PbP_DP .* yPb_by_grp) + sum(PbO_DP .* yOb_by_grp);

CostB_DP = bus_oper_cost( ...
    fDP,  pDP, 0, 0, YADP,   YBDP, ...
    fopDP,pDP, 0, 0, YAopDP, YBopDP, ...
    kDP, L,a,Q,abpr,bbpr,tp,tsb, ...
    h_peak, h_off, Gb0,Gb1,Gv0,Gv1, cost_amp);

CS_DP = consumer_surplus(outDP, Y_grp, bet_cost, mu, w, th_nt);

% Velocidades, tamaño operativo y shares
VelAp_DP  = 1 ./ ta_km(fDP,   pDP, 0, 0, YADP,   YBDP,   kDP, L,a,Q,abpr,bbpr,tp,tsb);
VelBp_DP  = 1 ./ tb_km(fDP,   pDP, 0, 0, YADP,   YBDP,   kDP, L,a,Q,abpr,bbpr,tp,tsb);
VelAop_DP = 1 ./ ta_km(fopDP, pDP, 0, 0, YAopDP, YBopDP, kDP, L,a,Q,abpr,bbpr,tp,tsb);
VelBop_DP = 1 ./ tb_km(fopDP, pDP, 0, 0, YAopDP, YBopDP, kDP, L,a,Q,abpr,bbpr,tp,tsb);

k1_DP     = (YBDP * L)/max(fDP,1e-12);

Pshare_DP  = (YADP + YBDP)*h_peak / Y_day;
OPshare_DP = (YAopDP + YBopDP)*h_off / Y_day;
NTshare_DP = 1 - Pshare_DP - OPshare_DP;
CPshare_DP = YADP/(YADP + YBDP);
BPshare_DP = YBDP/(YADP + YBDP);
COPshare_DP= YAopDP/(YAopDP + YBopDP);
BOPshare_DP= YBopDP/(YAopDP + YBopDP);

SW_DP_from_fval_local = -fval_pure;
SW_DP_recalc_local    = social_welfare(CS_DP, RevB_DP, 0, CostB_DP, mcpf, 0);

fprintf('\n=== DEMANDA PURA (tarifas diferenciadas, autofinanciamiento) ===\n');
fprintf(['f=%.4f fop=%.4f p=%.4f YA=%.3f YB=%.3f YAop=%.3f YBop=%.3f ' ...
         'k=%.3f k1=%.3f | exitflag=%d\n'], ...
        fDP, fopDP, pDP, YADP, YBDP, YAopDP, YBopDP, kDP, k1_DP, exitflag_pure);
fprintf('VelAp=%.3f VelBp=%.3f VelAop=%.3f VelBop=%.3f\n', ...
        VelAp_DP, VelBp_DP, VelAop_DP, VelBop_DP);
fprintf(['CS=%.2f | SW(from fval)=%.2f | SW(recalc)=%.2f | RevB=%.2f | ' ...
         'RevC=%.2f | CostB=%.2f | DL=%.2f\n'], ...
        CS_DP, SW_DP_from_fval_local, SW_DP_recalc_local, ...
        RevB_DP, 0, CostB_DP, 0);

% Tabla: tarifas por grupo y flujos diarios de bus
grp = (1:numel(Y_grp)).';
T_pure = table( ...
    grp, ...
    PbP_DP.', PbO_DP.', ...
    yPb_by_grp.', yOb_by_grp.', ...
    'VariableNames', {'Grupo','Pb_Peak','Pb_Off','YBus_Peak_day','YBus_Off_day'});
disp(T_pure);

%% Funciones helper DEMANDA PURA
function val = obj_pure_demanda(x, ...
    L,a,Q,abpr,bbpr,tp,tsb,v_walk, ...
    th_ap,th_bp,th_aop,th_bop,th_nt, ...
    bet_T_ap,bet_T_bp,bet_T_aop,bet_T_bop, ...
    bet_cost, mu, c0c, ...
    Y_grp, h_peak, h_off, ...
    Gb0,Gb1,Gv0,Gv1, cost_amp, mcpf, w)
% OBJ_PURE_DEMANDA: devuelve -SW(x) para DEMANDA PURA

    f    = x(1);   fop  = x(2);   p    = x(3);
    YA   = x(4);   YB   = x(5);
    YAop = x(6);   YBop = x(7);
    k    = x(8);
    PbP  = (x(9:13)).';
    PbO  = (x(14:18)).';

    out = demand_nl( ...
        f,   p, 0, PbP,  YA,   YB, ...
        fop, p, 0, PbO, YAop, YBop, ...
        k, L,a,Q,abpr,bbpr,tp,tsb,v_walk, ...
        th_ap,th_bp,th_aop,th_bop,th_nt, ...
        bet_T_ap,bet_T_bp,bet_T_aop,bet_T_bop, bet_cost, mu, c0c);

    CS = consumer_surplus(out, Y_grp, bet_cost, mu, w, th_nt);

    yPb_by_grp = Y_grp .* out.P.nest.P .* out.P.cond.P.bus;
    yOb_by_grp = Y_grp .* out.P.nest.O .* out.P.cond.O.bus;

    RevB = sum(PbP .* yPb_by_grp) + sum(PbO .* yOb_by_grp);
    RevC = 0;
    DLcost = 0;

    CostB = bus_oper_cost( ...
        f,   p, 0, 0, YA,   YB, ...
        fop, p, 0, 0, YAop, YBop, ...
        k, L,a,Q,abpr,bbpr,tp,tsb, ...
        h_peak, h_off, Gb0,Gb1,Gv0,Gv1, cost_amp);

    SW = social_welfare(CS, RevB, RevC, CostB, mcpf, DLcost);
    val = -SW;
end

function [c,ceq] = nonlcon_pure_demanda(x, ...
    L,a,Q,abpr,bbpr,tp,tsb,v_walk, ...
    th_ap,th_bp,th_aop,th_bop,th_nt, ...
    bet_T_ap,bet_T_bp,bet_T_aop,bet_T_bop, ...
    bet_cost, mu, c0c, ...
    Y_grp, h_peak, h_off, ...
    Gb0,Gb1,Gv0,Gv1, cost_amp)
% Restricciones para DEMANDA PURA:
% capacidad de bus/paradero, consistencia de flujos y autofinanciamiento.

    f    = x(1);   fop  = x(2);   p    = x(3);
    YA   = x(4);   YB   = x(5);
    YAop = x(6);   YBop = x(7);
    k    = x(8);
    PbP  = (x(9:13)).';
    PbO  = (x(14:18)).';

    % Capacidad por tamaño de bus
    c1 = (YB   * L)/max(f,   1e-12) - k;
    c2 = (YBop * L)/max(fop, 1e-12) - k;

    % Capacidad de paradero
    capP = QP(f,   p, 0, 0, YA,   YB);
    capO = QP(fop, p, 0, 0, YAop, YBop);
    c3 = f   - capP;
    c4 = fop - capO;
    c  = [c1; c2; c3; c4];

    % Consistencia de flujos YA/YB con el nested logit
    out = demand_nl( ...
        f,   p, 0, PbP,  YA,   YB, ...
        fop, p, 0, PbO, YAop, YBop, ...
        k, L,a,Q,abpr,bbpr,tp,tsb,v_walk, ...
        th_ap,th_bp,th_aop,th_bop,th_nt, ...
        bet_T_ap,bet_T_bp,bet_T_aop,bet_T_bop, bet_cost, mu, c0c);

    YA_day_P = sum(Y_grp .* out.P.nest.P .* out.P.cond.P.car);
    YB_day_P = sum(Y_grp .* out.P.nest.P .* out.P.cond.P.bus);
    YA_day_O = sum(Y_grp .* out.P.nest.O .* out.P.cond.O.car);
    YB_day_O = sum(Y_grp .* out.P.nest.O .* out.P.cond.O.bus);

    e1 = YA   * h_peak - YA_day_P;
    e2 = YB   * h_peak - YB_day_P;
    e3 = YAop * h_off  - YA_day_O;
    e4 = YBop * h_off  - YB_day_O;

    % Autofinanciamiento: RevB = CostB
    yPb_by_grp = Y_grp .* out.P.nest.P .* out.P.cond.P.bus;
    yOb_by_grp = Y_grp .* out.P.nest.O .* out.P.cond.O.bus;

    RevB = sum(PbP .* yPb_by_grp) + sum(PbO .* yOb_by_grp);

    CostB = bus_oper_cost( ...
        f,   p, 0, 0, YA,   YB, ...
        fop, p, 0, 0, YAop, YBop, ...
        k, L,a,Q,abpr,bbpr,tp,tsb, ...
        h_peak, h_off, Gb0,Gb1,Gv0,Gv1, cost_amp);

    e5 = RevB - CostB;

    ceq = [e1; e2; e3; e4; e5];
end

%% DEMANDA PURA vs REF65: barra de ΔSW

if ~exist('SWref','var')
    error('No existe SWref. Ejecuta primero la sección REF65 (baseline).');
end
if ~exist('fval_pure','var')
    error('No existe fval_pure. Ejecuta primero la sección DEMANDA PURA.');
end

SW_DP_from_fval = -fval_pure;

RevC_DP   = 0;
DLcost_DP = 0;
SW_DP = social_welfare(CS_DP, RevB_DP, RevC_DP, CostB_DP, mcpf, DLcost_DP);

fprintf('\n[CHECK] SW_DP_from_fval = %.4f | SW_DP(recalc) = %.4f | diff = %.2e\n', ...
        SW_DP_from_fval, SW_DP, SW_DP_from_fval - SW_DP);

dSW_DP = SW_DP - SWref;

fprintf('SW_REF65 = %.4f | SW_DemandaPura = %.4f | ΔSW = %.4f\n', ...
        SWref, SW_DP, dSW_DP);

figure('Color','w');
bar(1, dSW_DP, 'FaceColor',[0.85 0.10 0.10]);
hold on;
text(1, dSW_DP, sprintf('  %.2f', dSW_DP), ...
    'VerticalAlignment', ternary(dSW_DP>=0,'bottom','top'), ...
    'HorizontalAlignment','left', ...
    'FontSize',10, ...
    'Color',[0.2 0.2 0.2]);
yline(0,'--','Color',[0.3 0.3 0.3]);
xlim([0.5 1.5]);
set(gca,'XTick',1,'XTickLabel',{'Demanda pura vs REF65'});
ylabel('\Delta SW = SW^{DP} - SW^{REF65}');
title('Cambio en bienestar social respecto al caso base (REF65)');
grid on;
hold off;

ResumenDP = table(SWref, SW_DP, dSW_DP, ...
    'VariableNames', {'SW_REF65','SW_DemandaPura','DeltaSW'});
disp('=== Resumen bienestar: REF65 vs Demanda Pura ===');
disp(ResumenDP);

%% DEMANDA PURA (tabla resumen estilo SUBX, base REF65)

Sub_DP = 0;   % sin subsidio a la oferta

Pb_DP_scalar   = mean(PbP_DP);
Pbop_DP_scalar = mean(PbO_DP);
BusFarePeak_perkm_DP = Pb_DP_scalar   / L;
BusFareOff_perkm_DP  = Pbop_DP_scalar / L;

CarTollPeak_DP = 0;
CarTollOff_DP  = 0;

T_DP = table( ...
    Sub_DP, SW_DP, CS_DP, dSW_DP, ...
    fDP, fopDP, pDP, ...
    Pb_DP_scalar, Pbop_DP_scalar, BusFarePeak_perkm_DP, BusFareOff_perkm_DP, ...
    YADP, YBDP, YAopDP, YBopDP, kDP, ...
    VelAp_DP, VelBp_DP, VelAop_DP, VelBop_DP, ...
    Pshare_DP, OPshare_DP, NTshare_DP, ...
    CPshare_DP, BPshare_DP, COPshare_DP, BOPshare_DP, ...
    CarTollPeak_DP, CarTollOff_DP, ...
    'VariableNames', { ...
      'Sub','SW','CS','dSW', ...
      'f','fop','p', ...
      'Pb','Pbop','BusFarePeak_perkm','BusFareOff_perkm', ...
      'YA','YB','YAop','YBop','BusSize', ...
      'VelAp','VelBp','VelAop','VelBop', ...
      'PeakShare','OffShare','NTshare', ...
      'CarShare_P','BusShare_P','CarShare_O','BusShare_O', ...
      'CarTollPeak','CarTollOff'});

disp('=== DEMANDA PURA (tabla resumen estilo SUBX, respecto a REF65) ===');
disp(T_DP);

function out = ternary(cond, a, b)
    if cond
        out = a;
    else
        out = b;
    end
end

%% MIX — Política mixta (tarifas diferenciadas + subsidio a la oferta)
% Mismo vector de decisión que DEMANDA PURA:
% [f fop p YA YB YAop YBop k PbP_1..PbP_5 PbO_1..PbO_5]

if ~exist('SWref','var')
    error('No existe SWref (REF65). Ejecuta primero la sección REF65.');
end
if ~exist('x_pure','var') || ~exist('CS_DP','var') || ~exist('RevB_DP','var') || ~exist('CostB_DP','var')
    error('Falta solución de DEMANDA PURA. Ejecuta primero la sección DEMANDA PURA.');
end

RevC_DP   = 0;
DLcost_DP = 0;
SW_DP = social_welfare(CS_DP, RevB_DP, RevC_DP, CostB_DP, mcpf, DLcost_DP);
fprintf('\n[INFO MIX] SW_DP (Demanda Pura, s=0) = %.4f\n', SW_DP);

submax_mix = 100;
step_mix   = 5;
subs_mix   = (0:step_mix:submax_mix).';
nMix       = numel(subs_mix);

% Preasignación MIX
f_mix     = nan(nMix,1);
fop_mix   = nan(nMix,1);
p_mix     = nan(nMix,1);
YA_mix    = nan(nMix,1);
YB_mix    = nan(nMix,1);
YAop_mix  = nan(nMix,1);
YBop_mix  = nan(nMix,1);
k_mix     = nan(nMix,1);

SW_mix    = nan(nMix,1);
CS_mix    = nan(nMix,1);
RevB_mix  = nan(nMix,1);
CostB_mix = nan(nMix,1);
exitflag_mix = nan(nMix,1);

PbP_mix   = nan(nMix,5);
PbO_mix   = nan(nMix,5);

x_mix_all = nan(18, nMix);

VelAp_mix_all  = nan(nMix,1);
VelBp_mix_all  = nan(nMix,1);
VelAop_mix_all = nan(nMix,1);
VelBop_mix_all = nan(nMix,1);
Pshare_mix     = nan(nMix,1);
OPshare_mix    = nan(nMix,1);
NTshare_mix    = nan(nMix,1);
CPshare_mix    = nan(nMix,1);
BPshare_mix    = nan(nMix,1);
COPshare_mix   = nan(nMix,1);
BOPshare_mix   = nan(nMix,1);
k1_mix_all     = nan(nMix,1);

lb_mix   = lb_pure;
ub_mix   = ub_pure;
opts_mix = opts_pure;

objMix = objDP;

for t = 1:nMix
    sPerc = subs_mix(t);
    sMix  = sPerc/100;

    fprintf('\n=== MIX: resolviendo para s = %.0f%% ===\n', sPerc);

    if t == 1
        x0_mix = x_pure;
    else
        x0_mix = x_mix_all(:, t-1);
    end

    nonlconMix = @(x) nonlcon_mix_demanda( ...
        x, sMix, ...
        L,a,Q,abpr,bbpr,tp,tsb,v_walk, ...
        th_ap,th_bp,th_aop,th_bop,th_nt, ...
        bet_T_ap,bet_T_bp,bet_T_aop,bet_T_bop, ...
        bet_cost, mu, c0c, ...
        Y_grp, h_peak, h_off, ...
        Gb0,Gb1,Gv0,Gv1, cost_amp);

    [x_opt, fval_opt, exitflag_opt, ~] = ...
        fmincon(objMix, x0_mix, [],[], [],[], lb_mix, ub_mix, nonlconMix, opts_mix);

    x_mix_all(:,t) = x_opt;
    exitflag_mix(t) = exitflag_opt;

    f_mix(t)     = x_opt(1);
    fop_mix(t)   = x_opt(2);
    p_mix(t)     = x_opt(3);
    YA_mix(t)    = x_opt(4);
    YB_mix(t)    = x_opt(5);
    YAop_mix(t)  = x_opt(6);
    YBop_mix(t)  = x_opt(7);
    k_mix(t)     = x_opt(8);
    PbP_mix(t,:) = (x_opt(9:13)).';
    PbO_mix(t,:) = (x_opt(14:18)).';

    if abs(sMix - 1) < 1e-12
        PbP_mix(t,:) = 0;
        PbO_mix(t,:) = 0;
    end

    outM = demand_nl( ...
        f_mix(t),   p_mix(t), 0, PbP_mix(t,:), YA_mix(t),   YB_mix(t), ...
        fop_mix(t), p_mix(t), 0, PbO_mix(t,:), YAop_mix(t), YBop_mix(t), ...
        k_mix(t), L,a,Q,abpr,bbpr,tp,tsb,v_walk, ...
        th_ap,th_bp,th_aop,th_bop,th_nt, ...
        bet_T_ap,bet_T_bp,bet_T_aop,bet_T_bop, bet_cost, mu, c0c);

    CS_mix(t) = consumer_surplus(outM, Y_grp, bet_cost, mu, w, th_nt);

    yPb_by_grp = Y_grp .* outM.P.nest.P .* outM.P.cond.P.bus;
    yOb_by_grp = Y_grp .* outM.P.nest.O .* outM.P.cond.O.bus;

    RevB_mix(t) = sum(PbP_mix(t,:) .* yPb_by_grp) + sum(PbO_mix(t,:) .* yOb_by_grp);

    CostB_mix(t) = bus_oper_cost( ...
        f_mix(t),   p_mix(t), 0, 0, YA_mix(t),   YB_mix(t), ...
        fop_mix(t), p_mix(t), 0, 0, YAop_mix(t), YBop_mix(t), ...
        k_mix(t), L,a,Q,abpr,bbpr,tp,tsb, ...
        h_peak, h_off, Gb0,Gb1,Gv0,Gv1, cost_amp);

    % Shares y velocidades
    Pshare_mix(t)  = (YA_mix(t) + YB_mix(t))*h_peak / Y_day;
    OPshare_mix(t) = (YAop_mix(t) + YBop_mix(t))*h_off / Y_day;
    NTshare_mix(t) = 1 - Pshare_mix(t) - OPshare_mix(t);

    CPshare_mix(t)  = YA_mix(t)/(YA_mix(t)+YB_mix(t));
    BPshare_mix(t)  = YB_mix(t)/(YA_mix(t)+YB_mix(t));
    COPshare_mix(t) = YAop_mix(t)/(YAop_mix(t)+YBop_mix(t));
    BOPshare_mix(t) = YBop_mix(t)/(YAop_mix(t)+YBop_mix(t));

    VelAp_mix_all(t)  = 1 ./ ta_km(f_mix(t),   p_mix(t), 0, 0, YA_mix(t),   YB_mix(t),   k_mix(t), L,a,Q,abpr,bbpr,tp,tsb);
    VelBp_mix_all(t)  = 1 ./ tb_km(f_mix(t),   p_mix(t), 0, 0, YA_mix(t),   YB_mix(t),   k_mix(t), L,a,Q,abpr,bbpr,tp,tsb);
    VelAop_mix_all(t) = 1 ./ ta_km(fop_mix(t), p_mix(t), 0, 0, YAop_mix(t), YBop_mix(t), k_mix(t), L,a,Q,abpr,bbpr,tp,tsb);
    VelBop_mix_all(t) = 1 ./ tb_km(fop_mix(t), p_mix(t), 0, 0, YAop_mix(t), YBop_mix(t), k_mix(t), L,a,Q,abpr,bbpr,tp,tsb);

    k1_mix_all(t) = (YB_mix(t) * L)/max(f_mix(t),1e-12);

    RevC_mix   = 0;
    DLcost_mix = 0;

    SW_from_fval   = -fval_opt;
    SW_recalc      = social_welfare(CS_mix(t), RevB_mix(t), RevC_mix, CostB_mix(t), mcpf, DLcost_mix);
    SW_mix(t)      = SW_recalc;

    VelAp_mix  = VelAp_mix_all(t);
    VelBp_mix  = VelBp_mix_all(t);
    VelAop_mix = VelAop_mix_all(t);
    VelBop_mix = VelBop_mix_all(t);
    k1_mix     = k1_mix_all(t);

    fprintf('\n[MIX s=%3.0f%%]\n', sPerc);
    fprintf(['f=%.4f fop=%.4f p=%.4f YA=%.3f YB=%.3f YAop=%.3f YBop=%.3f ' ...
             'k=%.3f k1=%.3f | exitflag=%d\n'], ...
            f_mix(t), fop_mix(t), p_mix(t), ...
            YA_mix(t), YB_mix(t), YAop_mix(t), YBop_mix(t), ...
            k_mix(t), k1_mix, exitflag_opt);
    fprintf('Tarifas peak PbP_g1..g5 = [%.4f %.4f %.4f %.4f %.4f]\n', ...
            PbP_mix(t,1), PbP_mix(t,2), PbP_mix(t,3), PbP_mix(t,4), PbP_mix(t,5));
    fprintf('Tarifas off  PbO_g1..g5 = [%.4f %.4f %.4f %.4f %.4f]\n', ...
            PbO_mix(t,1), PbO_mix(t,2), PbO_mix(t,3), PbO_mix(t,4), PbO_mix(t,5));
    fprintf('VelAp=%.3f VelBp=%.3f VelAop=%.3f VelBop=%.3f\n', ...
            VelAp_mix, VelBp_mix, VelAop_mix, VelBop_mix);
    fprintf(['CS=%.2f | SW(from fval)=%.2f | SW(recalc)=%.2f | RevB=%.2f | ' ...
             'RevC=%.2f | CostB=%.2f | DL=%.2f\n'], ...
            CS_mix(t), SW_from_fval, SW_recalc, ...
            RevB_mix(t), 0, CostB_mix(t), 0);
end

% ΔSW vs Demanda Pura
dSW_mix_DP = SW_mix - SW_DP;

idx0_mix = find(abs(subs_mix-0)<1e-12, 1);
if ~isempty(idx0_mix)
    fprintf('\n[CHECK MIX] s=0%% → ΔSW = %.3e\n', dSW_mix_DP(idx0_mix));
end

figure('Color','w');
plot(subs_mix, dSW_mix_DP, '-o', ...
    'LineWidth',1.5, ...
    'MarkerFaceColor',[0.10 0.50 0.10], ...
    'Color',[0.10 0.50 0.10]);
grid on; box on;
xlabel('Subsidio a la oferta s (%)');
ylabel('\Delta SW = SW_{MIX}(s) - SW_{DP}(s=0)');
title('Política mixta: bienestar vs. Demanda Pura');
yline(0,'--','Color',[0.3 0.3 0.3]);

T_mix = table( ...
    subs_mix, SW_mix, dSW_mix_DP, CS_mix, RevB_mix, CostB_mix, ...
    f_mix, fop_mix, p_mix, k_mix, ...
    PbP_mix(:,1), PbP_mix(:,2), PbP_mix(:,3), PbP_mix(:,4), PbP_mix(:,5), ...
    PbO_mix(:,1), PbO_mix(:,2), PbO_mix(:,3), PbO_mix(:,4), PbO_mix(:,5), ...
    'VariableNames', { ...
        'Subsidio_s','SW','DeltaSW_vs_DP','CS','RevB','CostB', ...
        'f','fop','p','BusSize_k', ...
        'PbP_g1','PbP_g2','PbP_g3','PbP_g4','PbP_g5', ...
        'PbO_g1','PbO_g2','PbO_g3','PbO_g4','PbO_g5'});

disp('=== Resultados MIX (tarifas diferenciadas + subsidio a la oferta) ===');
disp(T_mix);

%% MIX vs REF65: ΔSW por subsidio

if ~exist('SW_mix','var')
    error('No existe SW_mix. Ejecuta primero el barrido MIX.');
end
if ~exist('SWref','var')
    error('No existe SWref. Ejecuta primero REF65.');
end

if exist('subs_mix','var') && ~isempty(subs_mix)
    subsM = subs_mix(:);
else
    Nmix  = numel(SW_mix);
    subsM = linspace(0,100,Nmix).';
end

SW_mix = SW_mix(:);
if exist('CS_mix','var') && ~isempty(CS_mix)
    CS_mix = CS_mix(:);
else
    CS_mix = nan(size(SW_mix));
end

DeltaSW_mix_REF65 = SW_mix - SWref;

idx65_mix = find(abs(subsM-65)<1e-12, 1);

fprintf('Check MIX vs REF65: SWref = %.4f\n', SWref);
if ~isempty(idx65_mix)
    fprintf('  En MIX, s=65%% → SW_MIX(65)=%.4f | ΔSW=%.4f\n', ...
        SW_mix(idx65_mix), DeltaSW_mix_REF65(idx65_mix));
end

figure('Color','w');
hold on;
plot(subsM, DeltaSW_mix_REF65, '-', ...
    'Color',[0.20 0.40 0.70], ...
    'LineWidth',2.0);
scatter(subsM, DeltaSW_mix_REF65, 36, ...
    'MarkerEdgeColor',[0.20 0.40 0.70], ...
    'MarkerFaceColor',[0.80 0.88 1.00], ...
    'LineWidth',1);
yline(0,'--','Color',[0.4 0.4 0.4],'LineWidth',1);

if ~isempty(idx65_mix)
    x65 = subsM(idx65_mix);
    y65 = DeltaSW_mix_REF65(idx65_mix);
    scatter(x65, y65, 64, [0.90 0.30 0.10], 'filled');
    text(x65, y65, sprintf('  s=65%%, \\DeltaSW=%.1f', y65), ...
        'VerticalAlignment','bottom', ...
        'HorizontalAlignment','left', ...
        'FontSize',10, ...
        'Color',[0.2 0.2 0.2]);
end

xlim([min(subsM)-2, max(subsM)+2]);
yPad = 0.05 * max(1, range(DeltaSW_mix_REF65));
ylim([min(DeltaSW_mix_REF65)-yPad, max(DeltaSW_mix_REF65)+yPad]);

grid on; box on;
set(gca, ...
    'FontName','Helvetica', ...
    'FontSize',11, ...
    'LineWidth',1);

xlabel('Subsidio a la oferta s (%)', 'FontWeight','bold');
ylabel('\Delta SW_{MIX} - SW_{REF65}', ...
       'Interpreter','tex', ...
       'FontWeight','bold');
title('Política mixta vs caso base REF65 (s=65%)', ...
      'FontWeight','bold');
hold off;

CarTollPeak_mix = zeros(size(subsM));
CarTollOff_mix  = zeros(size(subsM));

T_mix_REF65 = table( ...
    subsM, SW_mix, CS_mix, DeltaSW_mix_REF65, ...
    f_mix, fop_mix, p_mix, k_mix, ...
    PbP_mix(:,1), PbP_mix(:,2), PbP_mix(:,3), PbP_mix(:,4), PbP_mix(:,5), ...
    PbO_mix(:,1), PbO_mix(:,2), PbO_mix(:,3), PbO_mix(:,4), PbO_mix(:,5), ...
    VelAp_mix_all, VelBp_mix_all, VelAop_mix_all, VelBop_mix_all, ...
    Pshare_mix, OPshare_mix, NTshare_mix, ...
    CPshare_mix, BPshare_mix, COPshare_mix, BOPshare_mix, ...
    YA_mix, YB_mix, YAop_mix, YBop_mix, ...
    CarTollPeak_mix, CarTollOff_mix, ...
    'VariableNames', { ...
      'Sub','SW','CS','dSW', ...
      'f','fop','p','BusSize', ...
      'PbP_g1','PbP_g2','PbP_g3','PbP_g4','PbP_g5', ...
      'PbO_g1','PbO_g2','PbO_g3','PbO_g4','PbO_g5', ...
      'VelAp','VelBp','VelAop','VelBop', ...
      'PeakShare','OffShare','NTshare', ...
      'CarShare_P','BusShare_P','CarShare_O','BusShare_O', ...
      'YA','YB','YAop','YBop', ...
      'CarTollPeak','CarTollOff'});

disp('=== Resultados política MIX vs caso base REF65 (tabla estilo SUBX) ===');
disp(T_mix_REF65);

%% Diccionario de variables a Excel

filename = 'resultados_tesis.xlsx';

Sheet   = { ...
    'SUBX'; 'SUBX'; 'SUBX'; 'SUBX'; ...
    'MIX_vs_REF65'; 'MIX_vs_REF65'; ...
    'DEMANDA_PURA'; 'DEMANDA_PURA' ...
    };

VarName = { ...
    'Sub'; 'SW'; 'CS'; 'dSW'; ...
    'Sub'; 'dSW'; ...
    'Grupo'; 'Pb_Peak' ...
    };

Unidad = { ...
    '%'; '$/día'; '$/día'; '$/día'; ...
    '%'; '$/día'; ...
    'índice'; '$/viaje' ...
    };

Descripcion = { ...
    'Subsidio a la oferta en porcentaje'; ...
    'Bienestar social diario'; ...
    'Excedente del consumidor diario'; ...
    'ΔSW respecto al caso base de SUBX'; ...
    'Subsidio a la oferta en MIX'; ...
    'ΔSW respecto a REF65'; ...
    'Índice de grupo (1 a 5)'; ...
    'Tarifa bus en punta para el grupo' ...
    };

Dicc = table(Sheet, VarName, Unidad, Descripcion);
writetable(Dicc, filename, 'Sheet','Diccionario');

%% Helper MIX
function [c,ceq] = nonlcon_mix_demanda(x, s, ...
    L,a,Q,abpr,bbpr,tp,tsb,v_walk, ...
    th_ap,th_bp,th_aop,th_bop,th_nt, ...
    bet_T_ap,bet_T_bp,bet_T_aop,bet_T_bop, ...
    bet_cost, mu, c0c, ...
    Y_grp, h_peak, h_off, ...
    Gb0,Gb1,Gv0,Gv1, cost_amp)
% NONLCON_MIX_DEMANDA: igual a NONLCON_PURE_DEMANDA pero con
% RevB = (1-s)*CostB.

    f    = x(1);   fop  = x(2);   p    = x(3);
    YA   = x(4);   YB   = x(5);
    YAop = x(6);   YBop = x(7);
    k    = x(8);
    PbP  = (x(9:13)).';
    PbO  = (x(14:18)).';

    % Capacidad por tamaño de bus
    c1 = (YB   * L)/max(f,   1e-12) - k;
    c2 = (YBop * L)/max(fop, 1e-12) - k;

    % Capacidad de paradero
    capP = QP(f,   p, 0, 0, YA,   YB);
    capO = QP(fop, p, 0, 0, YAop, YBop);
    c3 = f   - capP;
    c4 = fop - capO;

    c  = [c1; c2; c3; c4];

    % Consistencia de flujos
    out = demand_nl( ...
        f,   p, 0, PbP,  YA,   YB, ...
        fop, p, 0, PbO, YAop, YBop, ...
        k, L,a,Q,abpr,bbpr,tp,tsb,v_walk, ...
        th_ap,th_bp,th_aop,th_bop,th_nt, ...
        bet_T_ap,bet_T_bp,bet_T_aop,bet_T_bop, bet_cost, mu, c0c);

    YA_day_P = sum(Y_grp .* out.P.nest.P .* out.P.cond.P.car);
    YB_day_P = sum(Y_grp .* out.P.nest.P .* out.P.cond.P.bus);
    YA_day_O = sum(Y_grp .* out.P.nest.O .* out.P.cond.O.car);
    YB_day_O = sum(Y_grp .* out.P.nest.O .* out.P.cond.O.bus);

    e1 = YA   * h_peak - YA_day_P;
    e2 = YB   * h_peak - YB_day_P;
    e3 = YAop * h_off  - YA_day_O;
    e4 = YBop * h_off  - YB_day_O;

    % Balance financiero con subsidio s
    yPb_by_grp = Y_grp .* out.P.nest.P .* out.P.cond.P.bus;
    yOb_by_grp = Y_grp .* out.P.nest.O .* out.P.cond.O.bus;

    RevB = sum(PbP .* yPb_by_grp) + sum(PbO .* yOb_by_grp);

    CostB = bus_oper_cost( ...
        f,   p, 0, 0, YA,   YB, ...
        fop, p, 0, 0, YAop, YBop, ...
        k, L,a,Q,abpr,bbpr,tp,tsb, ...
        h_peak, h_off, Gb0,Gb1,Gv0,Gv1, cost_amp);

    e5 = RevB - (1 - s) * CostB;

    ceq = [e1; e2; e3; e4; e5];
end

%% Funciones de red y tiempos
function val = QP(f, p, Pa, Pb, YA, YB)
% Capacidad efectiva de paradero (buses/h) en función de carga de pasajeros.
x = YB ./ (f .* p);
x = max(x, 1e-12);
val = ((1.9015 .* x - 70.56) .* log(x)) - 5.5037 .* x + 251.57;
end

function del_int = de(f, p, Pa, Pb, YA, YB)
% Demora interna en parada (h/bus)
fp = f .* p;
fp = max(fp, 1e-12);
x  = YB ./ fp;

cap = QP(f, p, Pa, Pb, YA, YB);
cap = max(cap, 1e-12);
qf  = f ./ cap;

del_int = ( (0.073396 + 0.77422*qf)*x + (0.13713 + 0.17426*qf)*x - 0.0065309*x^2 ) / 3600;
end

function del_queue = dq(f, p, Pa, Pb, YA, YB)
% Demora por cola en la parada (h/bus)
fp  = f .* p;
fp  = max(fp, 1e-12);
x   = YB ./ fp;

cap = QP(f,p,Pa,Pb,YA,YB);
cap = max(cap,1e-12);

del_queue = ((0.76032 + 0.15103*x) .* exp( f .* (4.9222 - 0.037423.*x) ./ cap )) / 3600;
end

function beta = b_equiv(k)
% Equivalencia bus→auto (vehículos equivalentes)
beta = 0.0114 * k + 1.15;
end

function e = epsilon_car(f)
% Fracción de operación de bus que afecta a los autos
e = 1 - 1 / (1.01^f);
end

function tcar_km = ta_km(f, p, Pa, Pb, YA, YB, k, L, a, Q, abpr, bbpr, tp, tsb)
% Tiempo por km en auto (h/km)

load = (L * (YA./a) + b_equiv(k) * f) ./ Q;
road = 0.01666 * (1 + abpr * (load.^bbpr));

del_int = de(f, p, Pa, Pb, YA, YB);
del_q   = dq(f, p, Pa, Pb, YA, YB);

busops = epsilon_car(f) * ( YB .* (tsb ./ max(f,1e-12)) + (tp + del_int + del_q) * p );

tcar_km = road + busops;
end

function tbus_km = tb_km(f, p, Pa, Pb, YA, YB, k, L, a, Q, abpr, bbpr, tp, tsb)
% Tiempo por km en bus (h/km)

load = (L * (YA./a) + b_equiv(k) * f) ./ Q;
road = 0.01666 * (1 + abpr * (load.^bbpr));

del_int = de(f, p, Pa, Pb, YA, YB);
del_q   = dq(f, p, Pa, Pb, YA, YB);

t_stops = ( YB .* (tsb ./ max(f,1e-12)) ) + (tp + del_int + del_q) * p;

tbus_km = road + t_stops;
end

function t = tva_trip(f, p, Pa, Pb, YA, YB, k, L, a, Q, abpr, bbpr, tp, tsb)
% Tiempo por viaje en auto
t = ta_km(f, p, Pa, Pb, YA, YB, k, L, a, Q, abpr,bbpr,tp,tsb) * L;
end

function t = tvb_trip(f, p, Pa, Pb, YA, YB, k, L, a, Q, abpr, bbpr, tp, tsb)
% Tiempo por viaje en bus
t = tb_km(f, p, Pa, Pb, YA, YB, k, L, a, Q, abpr,bbpr,tp,tsb) * L;
end

function tw = teb_wait(f)
% Espera media en parada (h)
tw = 1 ./ (2 * max(f,1e-12));
end

function tw = tab_walk(p, v_walk)
% Tiempo medio de caminata a la parada (h)
tw = 1 ./ (2 * max(p,1e-12) * v_walk);
end

function TG = TG_car(f, p, Pa, Pb, YA, YB, k, L, a, Q, abpr, bbpr, tp, tsb)
% Tiempo generalizado en auto (sólo tiempo en vehículo)
TG = tva_trip(f, p, Pa, Pb, YA, YB, k, L, a, Q, abpr,bbpr,tp,tsb);
end

function TG = TG_bus(f, p, Pa, Pb, YA, YB, k, L, a, Q, abpr, bbpr, tp, tsb, v_walk)
% Tiempo generalizado en bus (viaje + espera + caminata)
TG = tvb_trip(f, p, Pa, Pb, YA, YB, k, L, a, Q, abpr,bbpr,tp,tsb) ...
   + 1.93 * teb_wait(f) ...
   + 3.63 * tab_walk(p, v_walk);
end

function c = c_auto(Pa, L, c0c, a)
% Costo monetario por viaje en auto (por pasajero)
c = ( (Pa + c0c) * L ) / a;
end

function c = c_bus(Pb)
% Costo monetario por viaje en bus (tarifa pagada por pasajero)
c = Pb;
end

function out = demand_nl( ...
    fP,pP,PaP,PbP,YA_P,YB_P, ...
    fO,pO,PaO,PbO,YA_O,YB_O, ...
    k,L,a,Q,abpr,bbpr,tp,tsb,v_walk, ...
    th_ap,th_bp,th_aop,th_bop,th_nt, ...
    bet_T_ap,bet_T_bp,bet_T_aop,bet_T_bop, ...
    bet_cost, mu, c0c)

% Tiempos generalizados
TG_ap  = TG_car(fP,pP,PaP,PbP,YA_P,YB_P,k,L,a,Q,abpr,bbpr,tp,tsb);
TG_bp  = TG_bus(fP,pP,PaP,PbP,YA_P,YB_P,k,L,a,Q,abpr,bbpr,tp,tsb,v_walk);
TG_aop = TG_car(fO,pO,PaO,PbO,YA_O,YB_O,k,L,a,Q,abpr,bbpr,tp,tsb);
TG_bop = TG_bus(fO,pO,PaO,PbO,YA_O,YB_O,k,L,a,Q,abpr,bbpr,tp,tsb,v_walk);

% Costos monetarios
c_autoP  = c_auto(PaP, L, c0c, a);
c_autoO  = c_auto(PaO, L, c0c, a);
c_busP   = c_bus(PbP);
c_busO   = c_bus(PbO);

% Utilidades por grupo
Uap  = th_ap  + bet_T_ap  .* TG_ap  + bet_cost .* c_autoP;
Ubp  = th_bp  + bet_T_bp  .* TG_bp  + bet_cost .* c_busP;
Uaop = th_aop + bet_T_aop .* TG_aop + bet_cost .* c_autoO;
Ubop = th_bop + bet_T_bop .* TG_bop + bet_cost .* c_busO;
UO   = th_nt;

% Probabilidades condicionales modo|periodo
denP  = exp(Uap)  + exp(Ubp);
denOP = exp(Uaop) + exp(Ubop);
P_car_given_P = exp(Uap)  ./ denP;
P_bus_given_P = exp(Ubp)  ./ denP;
P_car_given_O = exp(Uaop) ./ denOP;
P_bus_given_O = exp(Ubop) ./ denOP;

% Logsums
Ap  = log(denP);
Aop = log(denOP);
AO  = UO;

% Probabilidades de nido
numP  = exp(mu * Ap);
numOP = exp(mu * Aop);
numNT = exp(mu * AO);
denN  = numP + numOP + numNT;

P_n_P  = numP  ./ denN;
P_n_O  = numOP ./ denN;
P_n_NT = numNT ./ denN;

% Salida
out.U.Uap  = Uap;   out.U.Ubp  = Ubp;   out.U.Uaop = Uaop;  out.U.Ubop = Ubop; out.U.UO = UO;
out.P.cond.P.car = P_car_given_P;  out.P.cond.P.bus = P_bus_given_P;
out.P.cond.O.car = P_car_given_O;  out.P.cond.O.bus = P_bus_given_O;
out.logsums.Ap = Ap; out.logsums.Aop = Aop; out.logsums.AO = AO;
out.P.nest.P  = P_n_P;  out.P.nest.O = P_n_O;  out.P.nest.NT = P_n_NT;
end

function flows = probs_to_flows(out, Y_grp)
% Conversión de probabilidades a flujos por modo/periodo

YA_P_by = Y_grp .* out.P.nest.P .* out.P.cond.P.car;
YB_P_by = Y_grp .* out.P.nest.P .* out.P.cond.P.bus;
YA_O_by = Y_grp .* out.P.nest.O .* out.P.cond.O.car;
YB_O_by = Y_grp .* out.P.nest.O .* out.P.cond.O.bus;
Y_NT_by = Y_grp .* out.P.nest.NT;

flows.YA_P = sum(YA_P_by);  flows.YB_P = sum(YB_P_by);
flows.YA_O = sum(YA_O_by);  flows.YB_O = sum(YB_O_by);
flows.Y_NT = sum(Y_NT_by);

flows.bygrp.YA_P = YA_P_by; flows.bygrp.YB_P = YB_P_by;
flows.bygrp.YA_O = YA_O_by; flows.bygrp.YB_O = YB_O_by;
flows.bygrp.Y_NT = Y_NT_by;

flows.Pp  = out.P.nest.P;
flows.Pop = out.P.nest.O;
flows.Pnt = out.P.nest.NT;
end

function CostB = bus_oper_cost( ...
    fP,pP,PaP,PbP,YA_P,YB_P, ...
    fO,pO,PaO,PbO,YA_O,YB_O, ...
    k,L,a,Q,abpr,bbpr,tp,tsb, ...
    hP,hO, Gb0,Gb1,Gv0,Gv1, cost_amp)

% Costo de operación del bus (día)
tbP_km = tb_km(fP,pP,PaP,PbP,YA_P,YB_P,k,L,a,Q,abpr,bbpr,tp,tsb);

Gb = Gb0 + Gb1*k;
Gv = Gv0 + Gv1*k;

CostB = cost_amp * ( ...
           fP * tbP_km * Gb * (hP + hO) + ...
           Gv * ( fP*hP + fO*hO ) );
end

function [RevB, RevC] = revenues(PbP,PbO, PaP,PaO, L,a, hP,hO, ...
                                 YB_P,YB_O, YA_P,YA_O, eta_rev)
% Ingresos de bus y peajes

RevB = PbP*YB_P*hP + PbO*YB_O*hO;

RevC_gross = (PaP*(L/a)*YA_P*hP) + (PaO*(L/a)*YA_O*hO);
RevC = (1 - eta_rev) * RevC_gross;
end

function CS = consumer_surplus(out, Y_grp, bet_cost, mu, w, th_nt)
% Excedente del consumidor con logsum nested logit

Ap  = out.logsums.Ap;   Aop = out.logsums.Aop;   AO = th_nt;
LS  = log( exp(mu.*Ap) + exp(mu.*Aop) + exp(mu.*AO) );
CS_i = (Y_grp ./ (-bet_cost)) .* (LS./mu);
CS   = sum(w .* CS_i);
end

function SW = social_welfare(CS, RevB, RevC, CostB, mcpf, DLcost)
% Bienestar social total
if nargin < 6, DLcost = 0; end
SW = CS + mcpf*( (RevB + RevC) - (CostB + DLcost) );
end

%% Funciones de REF65 (obj y restricciones)
function val = obj_ref(x, L,a,Q,abpr,bbpr,tp,tsb,v_walk, ...
    th_ap,th_bp,th_aop,th_bop,th_nt, ...
    bet_T_ap,bet_T_bp,bet_T_aop,bet_T_bop, bet_cost, mu, c0c, ...
    Y_grp, h_peak, h_off, ...
    Gb0,Gb1,Gv0,Gv1, cost_amp, eta_rev, mcpf, dl_cost, w)

    f=x(1); fop=x(2); p=x(3); Pb=x(4); Pbop=x(5);
    YA=x(6); YB=x(7); YAop=x(8); YBop=x(9); k=x(10);

    out = demand_nl(f,p,0,Pb, YA,YB,  fop,p,0,Pbop, YAop,YBop, ...
                k, L,a,Q,abpr,bbpr,tp,tsb,v_walk, ...
                th_ap,th_bp,th_aop,th_bop,th_nt, ...
                bet_T_ap,bet_T_bp,bet_T_aop,bet_T_bop, bet_cost, mu, c0c);

    [RevB, RevC] = revenues(Pb, Pbop, 0, 0, L,a, h_peak,h_off, YB,YBop, YA,YAop, eta_rev);
    CostB = bus_oper_cost( ...
                f,p,0,Pb, YA,YB,  fop,p,0,Pbop, YAop,YBop, ...
                k, L,a,Q,abpr,bbpr,tp,tsb, h_peak,h_off, Gb0,Gb1,Gv0,Gv1, cost_amp);

    CS = consumer_surplus(out, Y_grp, bet_cost, mu, w, th_nt);
    SW = social_welfare(CS, RevB, RevC, CostB, mcpf, dl_cost);

    val = -SW;
end

function [c,ceq] = nonlcon_ref(x, L,a,Q,abpr,bbpr,tp,tsb,v_walk, ...
    th_ap,th_bp,th_aop,th_bop,th_nt, ...
    bet_T_ap,bet_T_bp,bet_T_aop,bet_T_bop, bet_cost, mu, c0c, ...
    Y_grp, h_peak, h_off, ...
    Gb0,Gb1,Gv0,Gv1, cost_amp)

    f=x(1); fop=x(2); p=x(3); Pb=x(4); Pbop=x(5);
    YA=x(6); YB=x(7); YAop=x(8); YBop=x(9); k=x(10);

    % Capacidad por bus
    c1 = (YB  * L)/max(f,   1e-12) - k;
    c2 = (YBop* L)/max(fop, 1e-12) - k;

    % Capacidad de paradero
    capP = QP(f,   p, 0, Pb,   YA,   YB);
    capO = QP(fop, p, 0, Pbop, YAop, YBop);
    c3 = f   - capP;
    c4 = fop - capO;

    c = [c1; c2; c3; c4];

    % Flujos y balance operador sin subsidio
    out = demand_nl(f,p,0,Pb, YA,YB,  fop,p,0,Pbop, YAop,YBop, ...
                k, L,a,Q,abpr,bbpr,tp,tsb,v_walk, ...
                th_ap,th_bp,th_aop,th_bop,th_nt, ...
                bet_T_ap,bet_T_bp,bet_T_aop,bet_T_bop, bet_cost, mu, c0c);

    e1 = YA   * h_peak - sum(Y_grp .* out.P.nest.P .* out.P.cond.P.car);
    e2 = YB   * h_peak - sum(Y_grp .* out.P.nest.P .* out.P.cond.P.bus);
    e3 = YAop * h_off  - sum(Y_grp .* out.P.nest.O .* out.P.cond.O.car);
    e4 = YBop * h_off  - sum(Y_grp .* out.P.nest.O .* out.P.cond.O.bus);

    RevB = Pb*YB*h_peak + Pbop*YBop*h_off;
    CostB = bus_oper_cost( ...
                f,p,0,Pb, YA,YB,  fop,p,0,Pbop, YAop,YBop, ...
                k, L,a,Q,abpr,bbpr,tp,tsb, h_peak,h_off, Gb0,Gb1,Gv0,Gv1, cost_amp);
    e5 = RevB - CostB;

    ceq = [e1; e2; e3; e4; e5];
end
