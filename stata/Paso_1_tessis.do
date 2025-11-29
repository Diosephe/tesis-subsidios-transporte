*--------------------------------------------------
* Proyecto: Tesis
* TAREA: ENMODO-2018 (Hogar) → deciles NSE (1..10) → quintiles (1..5)
* AUTOR: Yoseph Barrera
*--------------------------------------------------

clear all
macro drop _all
capture log close
set more off

* Directorio de trabajo 
cd "C:\Users\drake\OneDrive\Escritorio\tesis\Nueva carpeta"

* (Opcional) log
cap mkdir "outputs"
log using "outputs\01_nse_quint.log", text replace

*--------------------------------------------------
	* ) IMPORTAR EXCEL (Hogar)
	import excel "Hogar_ENMODO18.xlsx", firstrow case(lower) clear
	keep id_hogar provincia_hogar partido_hogar radio_censal_hogar nse
	save "outputs/HOGAR_limpio.dta", replace

	* ) Viajes
	import excel "Viajes_ENMODO18.xlsx", firstrow case(lower) clear
	keep id_hogar id_persona id_viaje jur_origen jur_destino hora_inicio hora_final modo_des menos_2_cuadras cantidad_etapas tipo_viaje distancia_viaje_km duracion_minutos
	save "outputs/VIAJES_limpio.dta", replace

	* ) Etapas
	import excel "Etapas_ENMODO18.xlsx", firstrow case(lower) clear
	keep id_hogar id_persona id_viaje id_etapa ///
		 vii_10_medio_transporte vii_11_cuadras_caminadas ///
		 vii_26_cuanto_espero_el_transpor ///
		 vii_28_cuadras_caminadas_hasta_e ///
		 vii_30_cuadras_caminadas_al_baja
	save "outputs/ETAPAS_limpio.dta", replace

* ============================================
* Base final con TODO junto (conserva TODOS los hogares)
* ============================================
use "outputs/HOGAR_limpio.dta", clear
bys id_hogar: keep if _n==1

merge 1:m id_hogar using "outputs/VIAJES_limpio.dta", gen(_m_hv) keep(master match)
merge 1:m id_hogar id_persona id_viaje using "outputs/ETAPAS_limpio.dta", gen(_m_hve) keep(master match)

order id_hogar id_persona id_viaje id_etapa
save "outputs/EMD18_TODO_JUNTO.dta", replace

* (Chequeo rápido de hogares únicos)
egen tag_h = tag(id_hogar)
capture drop _m_hv _m_hve   // borrar indicadores de merge
count if tag_h
di as txt "Hogares únicos en EMD18_TODO_JUNTO: " r(N)
drop tag_h




*--------------------------------------------------
* 0) IMPORTAR 

use "outputs/EMD18_TODO_JUNTO.dta", clear
egen tag_h = tag(id_hogar)
count if tag_h
di as txt "Hogares únicos: " r(N)
drop tag_h

*--------------------------------------------------
* 1) LIMPIEZA (si nse llegó como string)
capture confirm string variable nse
if _rc==0 {
    replace nse = strtrim(nse)
    replace nse = subinstr(nse," ","",.)
    replace nse = "" if inlist(upper(nse),"NSE","NA","N/A","NAN",".","MISS","NULL")
    di as txt ">> Revisión tras limpieza:"
    count if missing(nse) | nse==""
    tab nse if missing(real(nse)) & nse!="" , miss
}
d

* Convertir NSE a numérico (crea nse_dec; conserva nse original)
capture drop nse_dec

capture confirm string variable nse
if _rc==0 {
    * nse es string: limpieza mínima y conversión
    replace nse = strtrim(nse)
    replace nse = subinstr(nse," ","",.)
    replace nse = "" if inlist(upper(nse),"NA","N/A","NAN",".","MISS","NULL","NSE")
    gen byte nse_dec = real(nse)
}
else {
    * nse ya es numérica: solo copiar al nuevo nombre (byte)
    gen byte nse_dec = nse
}

label var nse_dec "NSE (deciles 1..10, numérico)"
assert inrange(nse_dec,1,10) | missing(nse_dec)

* (opcional) chequeo rápido
describe nse nse_dec
tab nse_dec, miss
*------------------------------------------------------
*Grupos

* Mapear decil -> ingreso medio (misma unidad que tu fuente)
recode nse_dec ///
    (1 = 4279)   ///
    (2 = 8560)   ///
    (3 = 11571)  ///
    (4 = 14413)  ///
    (5 = 17932)  ///
    (6 = 21597)  ///
    (7 = 26961)  ///
    (8 = 33621)  ///
    (9 = 44240)  ///
    (10 = 87821), gen(ing_dec_medio)

label var ing_dec_medio "Ingreso medio por decil NSE"
format ing_dec_medio %12.3f

* Chequeo rápido
tab nse_dec, sum(ing_dec_medio)

*--------------------------------------------------
* Grupos 
*--------------------------------------------------
* (A) Construir el mapa grupo -> ingreso medio ponderado (USANDO HOGARES ÚNICOS)
preserve
keep id_hogar nse_dec ing_dec_medio
bys id_hogar: keep if _n==1                // una fila por hogar
keep if !missing(nse_dec, ing_dec_medio)

collapse (count) Ndec = id_hogar (mean) mean_dec = ing_dec_medio, by(nse_dec)

* Mapear decil -> grupo
gen byte grupo = .
replace grupo = 1 if inlist(nse_dec,1,2)
replace grupo = 2 if inlist(nse_dec,3,4)
replace grupo = 3 if inlist(nse_dec,5,6,7)
replace grupo = 4 if inlist(nse_dec,8,9)
replace grupo = 5 if nse_dec==10
label define G 1 "G1(1-2)" 2 "G2(3-4)" 3 "G3(5-7)" 4 "G4(8-9)" 5 "G5(10)", replace
label values grupo G

* Promedio ponderado por # de hogares dentro del grupo
gen double num = mean_dec * Ndec
collapse (sum) num (sum) den = Ndec, by(grupo)
gen double y_g = num/den
format y_g %9.3f

* (opcional) mirar frecuencias por grupo
egen long N_total = total(den)
gen double pct_grupo = 100*den/N_total
format pct_grupo %6.2f
list grupo y_g den pct_grupo, noobs

tempfile map_g
save `map_g', replace
restore

* (B) En tu base original: asignar grupo a cada observación
capture drop grupo
gen byte grupo = .
replace grupo = 1 if inlist(nse_dec,1,2)
replace grupo = 2 if inlist(nse_dec,3,4)
replace grupo = 3 if inlist(nse_dec,5,6,7)   // <- no olvides G3
replace grupo = 4 if inlist(nse_dec,8,9)
replace grupo = 5 if nse_dec==10
label values grupo G      // G ya quedó definida en (A)

* (C) Traer el ingreso medio del grupo (ponderado por hogares únicos) a cada fila
merge m:1 grupo using `map_g', nogen keep(match master)

rename y_g ing_grupo_medio
label var ing_grupo_medio "Ingreso medio del grupo (ponderado por hogares únicos)"
format ing_grupo_medio %9.3f
replace ing_grupo_medio = ing_grupo_medio/39


*----------------------------------------
preserve
keep id_hogar provincia_hogar partido_hogar radio_censal_hogar ///
     nse nse_dec ing_dec_medio grupo num den ing_grupo_medio
describe
restore

* --- Parámetros ---
local Hmes    = 208      // 48 h/sem (baseline legal). Cambia a 160 si quieres.
local fracVOT = 0.5      // VOT = 50% del salario

* (Opcional) si ing_grupo_medio está en miles de ARS, descomenta:
* gen double ing_grupo_medio_ars = 1000*ing_grupo_medio
* replace ing_grupo_medio = ing_grupo_medio_ars

* Chequeos
assert !missing(grupo, ing_grupo_medio)

* 1) Salario por hora y VOT por grupo (quedan en cada fila)
gen double w_hora_g   = ing_grupo_medio / `Hmes'
gen double VOT_h_g    = `fracVOT' * w_hora_g
gen double VOT_min_g  = VOT_h_g / 60

* (Opcional) betas si tus tiempos están en minutos y normalizas beta_C=-1
*gen double betaT_min_g = - VOT_min_g

* 2) Resumen único por grupo
preserve
collapse (mean) ing_grupo_medio w_hora_g VOT_h_g VOT_min_g , by(grupo)
format ing_grupo_medio w_hora_g VOT_h_g VOT_min_g %9.3f
list, noobs sepby(grupo)

* === Guardar VOT por grupo en USD/h para Script 3 ===
save "outputs/VOTg_USD_by_group.dta", replace
export delimited using "outputs/VOTg_USD_by_group.csv", replace
restore


*******----------- Sacar Frecuencias.
preserve
keep id_hogar grupo
bys id_hogar: keep if _n==1
tab grupo, missing   // muestra Frecuencia y Percent (sobre hogares)
restore
* -------------------------- Normalizar y parsear la hora de inicio ------------------------
* ------------------ Normalizar y clasificar usando SOLO hora_inicio ------------------
* 1) Elegir variable (solo hora_inicio)
local hvar "hora_inicio"
capture confirm variable `hvar'
if _rc {
    di as err "No existe la variable hora_inicio."
    exit 498
}

* 2) Pasar a numérico (horas decimales tipo 7.5, 14.75, etc.)
capture confirm numeric variable `hvar'
if _rc {
    gen strL _hi = `hvar'
    replace _hi = trim(_hi)
    replace _hi = subinstr(_hi, ",", ".", .)   // por si viene con coma decimal
    destring _hi, gen(__hi_num) force
    drop `hvar'
    rename __hi_num `hvar'
    drop _hi
}

* (sanity check) rango 0..24; 24 ≡ 0
assert `hvar'>=0 & `hvar'<=24 if !missing(`hvar')
replace `hvar' = mod(`hvar',24) if !missing(`hvar')

* 3) Clasificar periodos con horas decimales
capture drop periodo
gen byte periodo = .
replace periodo = 1 if (`hvar'>=7  & `hvar'<11) | (`hvar'>=16 & `hvar'<20)
replace periodo = 2 if missing(periodo) & !missing(`hvar')
label define P 1 "peak" 2 "offpeak", replace
label values periodo P

* 4) Chequeo rápido
count if missing(`hvar')
di as txt "Filas sin hora_inicio: " r(N)
tab periodo, missing
*==================== Preparar 'modo_des' y mapeos (TP=4,8; Auto=5; 7 excluido) ====================*

* Limpieza previa por si existen versiones anteriores
capture drop modo3 modo2

* Modo agregado a 3 categorías: TP / Auto / Otros (excluye 7 de TP)
gen byte modo3 = .
replace modo3 = 1 if inlist(modo_des,4,8)                 // Transporte público: SOLO códigos 4 y 8
replace modo3 = 2 if inlist(modo_des,5)                    // Auto (solo 5)
replace modo3 = 3 if inlist(modo_des,1,2,3,6,7,9)          // Otros (incluye 7 explícitamente)
label define M3 1 "TP(4,8)" 2 "Auto(5)" 3 "Otros(1,2,3,6,7,9)", replace
label values modo3 M3

* Modo binario para el modelo (TP vs Auto); deja missing a los 'Otros'
gen byte modo2 = .
replace modo2 = 1 if modo3==1      // TP
replace modo2 = 2 if modo3==2      // Auto
label define M2 1 "TP" 2 "Auto", replace
label values modo2 M2

* (Diagnóstico rápido)
tab modo_des modo3, m
tab modo2, m
count if modo_des==7
di as txt "Excluidos (modo==7): " %9.0f r(N)


*==================== Guardar la base hasta aquí ====================*
save "outputs/EMD18_model_base_TP48_AUTO5.dta", replace
export delimited using "outputs/EMD18_model_base_TP48_AUTO5.csv", replace

*========================================================
