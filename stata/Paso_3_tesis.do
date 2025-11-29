* ================== SCRIPT 3: dimensiones por grupo/modo/periodo (HORAS & USD) ==================
clear all
set more off
cd "C:\Users\drake\OneDrive\Escritorio\tesis\Nueva carpeta"
cap mkdir "outputs"

* 0) Cargar viajes enriquecidos (trae mins + horas)
use "outputs/EMD18_trip_enriched.dta", clear

* 1) Particiones por (g,m,t) y versión con NT=10%
preserve
    keep if !missing(grupo, periodo, modo2)
    contract grupo periodo modo2                   // _freq por (g,t,m)
    bys grupo periodo: egen N_gp = total(_freq)
    gen double s_within = _freq / N_gp
    bys grupo: egen N_g = total(_freq)
    gen double P_t = N_gp / N_g
    gen double p = s_within * P_t

    * Escalamos al 90% "viaja" y reservamos 10% "no viaja"
    gen double p_final = 0.90 * p

    * Mapeo de alternativas
    gen byte alt = .
    replace alt = 1 if modo2==2 & periodo==1   // Auto-peak
    replace alt = 2 if modo2==2 & periodo==2   // Auto-off
    replace alt = 3 if modo2==1 & periodo==1   // TP-peak
    replace alt = 4 if modo2==1 & periodo==2   // TP-off

    * (a) versión LONG
    save "outputs/particiones_finales_con_NT_long.dta", replace

    * (b) versión WIDE + NT
    keep grupo alt p p_final
    reshape wide p p_final, i(grupo) j(alt)
    gen double p_final_nt = 0.10
    order grupo p_final1 p_final2 p_final3 p_final4 p_final_nt
    save "outputs/particiones_finales_con_NT.dta", replace
    export delimited using "outputs/particiones_finales_con_NT.csv", replace
restore

* 2) Tiempo generalizado promedio por (g,m,t) EN HORAS
local phiW  = 1.93
local phiWK = 3.63
capture drop G_h
gen double G_h = T_vehicle_h ///
               + `phiW'  * cond(modo2==1, W_h, 0) ///
               + `phiWK' * cond(modo2==1, walk_h, 0)

preserve
    keep if inlist(modo2,1,2) & !missing(grupo, periodo, G_h)
    collapse (mean) G = G_h, by(grupo modo2 periodo)   // G en HORAS
    save "outputs/Gmt_by_g_m_t_hours.dta", replace
    export delimited using "outputs/Gmt_by_g_m_t_hours.csv", replace
restore

* 3) Costos por (g,m,t) — placeholder (0) por ahora
use "outputs/Gmt_by_g_m_t_hours.dta", clear
gen double C = 0
keep grupo modo2 periodo C
save "outputs/costs_by_g_m_t.dta", replace
export delimited using "outputs/costs_by_g_m_t.csv", replace

* 4) β_T por HORA: β_T,g = β_cost,g * VOT_h,g  (VOT_h,g en USD/h)
preserve
    use "outputs/EMD18_model_base_TP48_AUTO5.dta", clear
    keep grupo ing_grupo_medio
    bys grupo: keep if _n==1
    sort grupo

    * VOT por hora (50% del salario horario; 208 h/mes). ing_grupo_medio ya en USD.
    gen double w_hora_g = ing_grupo_medio / 208
    gen double VOT_h_g  = 0.5 * w_hora_g

    * β_cost por grupo (1..5) — tus valores:
 gen double bet_cost = .
    replace bet_cost = -1.70 if grupo==1
    replace bet_cost = -1.40 if grupo==2
    replace bet_cost = -1.20 if grupo==3
    replace bet_cost = -1.00 if grupo==4
    replace bet_cost = -0.72 if grupo==5
    * β_T por hora
    gen double bet_T_h = bet_cost * VOT_h_g

    * Guardar y exportar
    keep grupo bet_cost VOT_h_g bet_T_h
    save "outputs/betaT_hours_by_group.dta", replace
    export delimited using "outputs/betaT_hours_by_group.csv", replace

    * Línea lista para MATLAB
    levelsof grupo, local(grs)
    local vecstr ""
    foreach g of local grs {
        quietly su bet_T_h if grupo==`g', meanonly
        local vecstr `"`vecstr' `=string(r(mean),"%8.4f")'"'
    }
    di as res "MATLAB -> bet_T_h = [ `vecstr' ];   % (USD/h * util/$ = util/h)"
restore


* Más simple: rehacer vector y exportar limpio
preserve
    use "outputs/EMD18_model_base_TP48_AUTO5.dta", clear
    keep grupo ing_grupo_medio
    bys grupo: keep if _n==1
    sort grupo
    gen double w_hora_g = ing_grupo_medio / 208
    gen double VOT_h_g  = 0.5 * w_hora_g

    * asignar β_cost por grupo
  gen double bet_cost = .
    replace bet_cost = -1.70 if grupo==1
    replace bet_cost = -1.40 if grupo==2
    replace bet_cost = -1.20 if grupo==3
    replace bet_cost = -1.00 if grupo==4
    replace bet_cost = -0.72 if grupo==5
    * β_T por hora

    gen double bet_T_h = bet_cost * VOT_h_g

    * Guardar a disco
    keep grupo bet_cost VOT_h_g bet_T_h
    save "outputs/betaT_hours_by_group.dta", replace
    export delimited using "outputs/betaT_hours_by_group.csv", replace

    * Imprimir línea MATLAB
    levelsof grupo, local(grs)
    local vecstr ""
    foreach g of local grs {
        quietly su bet_T_h if grupo==`g', meanonly
        local bt = r(mean)
        local vecstr `"`vecstr' `=string(`bt',"%8.4f")'"'
    }
    di as res "MATLAB -> bet_T_h = [ `vecstr' ];   % (USD/h * util/$ = util/h)"
capture confirm variable bet_T_h
if _rc {
    di as err "No encuentro 'bet_T_h' en memoria. Corre primero el paso 4)."
    exit 459
}

* Factores relativos a Auto-peak (Santiago, BS2014)
local f_ap  = 1
local f_aop = 0.35/0.74
local f_bp  = 1.31/0.74
local f_bop = 0.86/0.74

* β_T por modo×periodo (util/h)
gen double bet_T_ap  = `f_ap'  * bet_T_h
gen double bet_T_aop = `f_aop' * bet_T_h
gen double bet_T_bp  = `f_bp'  * bet_T_h
gen double bet_T_bop = `f_bop' * bet_T_h

* Imprimir vectores en orden de grupo, estilo MATLAB
sort grupo
levelsof grupo, local(grs)

local vec_ap  ""
local vec_aop ""
local vec_bp  ""
local vec_bop ""

foreach g of local grs {
    quietly su bet_T_ap  if grupo==`g', meanonly
    local vec_ap  `"`vec_ap'  `=string(r(mean),"%8.4f")'"'
    quietly su bet_T_aop if grupo==`g', meanonly
    local vec_aop `"`vec_aop' `=string(r(mean),"%8.4f")'"'
    quietly su bet_T_bp  if grupo==`g', meanonly
    local vec_bp  `"`vec_bp'  `=string(r(mean),"%8.4f")'"'
    quietly su bet_T_bop if grupo==`g', meanonly
    local vec_bop `"`vec_bop' `=string(r(mean),"%8.4f")'"'
}

di as res "MATLAB -> bet_T_ap  = [ `vec_ap'  ];  % util/h"
di as res "MATLAB -> bet_T_aop = [ `vec_aop' ];  % util/h"
di as res "MATLAB -> bet_T_bp  = [ `vec_bp'  ];  % util/h"
di as res "MATLAB -> bet_T_bop = [ `vec_bop'  ];  % util/h"

* Limpieza visual (opcional)
drop bet_T_ap bet_T_aop bet_T_bp bet_T_bop


