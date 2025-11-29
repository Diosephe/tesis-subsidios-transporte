

clear all
macro drop _all
capture log close
set more off

* Directorio de trabajo 
cd "C:\Users\drake\OneDrive\Escritorio\tesis\Nueva carpeta"


*========================================================
* Base ya preparada (TP=4,8; Auto=5; peak/off definidos)
*========================================================
use "outputs/EMD18_model_base_TP48_AUTO5.dta", clear

*========================================================
* (A) Particiones por GRUPO — HOGARES y VIAJES (solo mostrar)
*========================================================
* Hogares por grupo (hogar único con id_hogar)
preserve
keep id_hogar grupo
bys id_hogar: keep if _n==1
keep if !missing(grupo)
contract grupo                          // _freq = N hogares en g
egen long HH_total = total(_freq)
gen double share_hh = _freq / HH_total
gen double pct_hh   = 100*share_hh
format share_hh %6.4f
format pct_hh   %6.2f
di as txt "=== Partición de hogares por grupo ==="
list grupo _freq HH_total share_hh pct_hh, noobs
restore

* Viajes por grupo (nivel viaje)
preserve
keep if !missing(grupo)
contract grupo                          // _freq = N viajes en g
egen long Y_day = total(_freq)
gen double share_trip = _freq / Y_day
gen double pct_trip   = 100*share_trip
gen double Y_grp      = share_trip * Y_day
format share_trip %6.4f
format pct_trip   %6.2f
di as txt "=== Partición de viajes por grupo ==="
list grupo _freq Y_day share_trip pct_trip Y_grp, noobs
restore

*========================================================
* (C) Trip-level limpio: dur_min (desde string) + espera (TP) + caminar  => T_vehicle_min
*========================================================

* -- 0) DURACIÓN TOTAL (min) desde 'duracion_minutos' (string) --
capture drop dur_min
capture confirm numeric variable duracion_minutos
if _rc {
    gen strL __s_dur = trim(duracion_minutos)
    replace __s_dur = subinstr(__s_dur, ",", ".", .)
    gen double dur_min = .
    replace dur_min = real(__s_dur) if regexm(__s_dur, "^[0-9]+(\\.[0-9]+)?$")
    replace dur_min = 60*real(regexs(1)) + real(regexs(2)) if ///
        regexm(__s_dur, "^([0-9]{1,2}):([0-9]{2})$")
    drop __s_dur
}
else {
    gen double dur_min = duracion_minutos
}
* Asegurar una única duración por viaje (evita duplicados por etapas)
bys id_hogar id_persona id_viaje: egen double dur_trip = max(dur_min)
replace dur_min = dur_trip
drop dur_trip

label var dur_min "Duración total del viaje (min)"
* -- 0.b) DISTANCIA a numérico (km) si existe -> dist_km
capture drop dist_km
local distvar ""
foreach cand in distancia_viaje_km distancia_km distancia_viaje distancia_kms distancia {
    capture confirm variable `cand'
    if !_rc & "`distvar'"=="" local distvar "`cand'"
}
if "`distvar'" != "" {
    capture confirm numeric variable `distvar'
    if _rc {
        gen strL __s_dist = trim(`distvar')
        replace __s_dist = subinstr(__s_dist, ",", ".", .)
        destring __s_dist, gen(dist_km) force
        drop __s_dist
    }
    else gen double dist_km = `distvar'
}


* -- 1) VIEW TRIP-LEVEL (una fila por viaje) con variables de viaje --
preserve
keep id_hogar id_persona id_viaje grupo periodo modo2 dur_min dist_km
bys id_hogar id_persona id_viaje: keep if _n==1
tempfile trips
save `trips', replace
restore


* -- 2) ESPERA (TP) desde ETAPAS (variable de espera en string -> num) y agregar por viaje --
preserve
keep if inlist(vii_10_medio_transporte,4,8)          // solo ETAPAS TP
* convertir espera por ETAPA a numérico (min)
tempvar w_etapa
capture confirm numeric variable vii_26_cuanto_espero_el_transpor
if _rc {
    gen strL __s_w = trim(vii_26_cuanto_espero_el_transpor)
    replace __s_w = subinstr(__s_w, ",", ".", .)
    gen double `w_etapa' = .
    replace `w_etapa' = real(__s_w) if regexm(__s_w, "^[0-9]+(\\.[0-9]+)?$")
    replace `w_etapa' = 60*real(regexs(1)) + real(regexs(2)) if ///
        regexm(__s_w, "^([0-9]{1,2}):([0-9]{2})$")
    drop __s_w
}
else {
    gen double `w_etapa' = vii_26_cuanto_espero_el_transpor
}
* promedio (o suma) de espera por VIAJE; usa mean para ser conservador
collapse (sum)  W_min = `w_etapa', by(id_hogar id_persona id_viaje)
tempfile waitsTP
save `waitsTP', replace
restore

* -- 3) CAMINAR (ACCESO/EGRESO) sin sobrecontar por etapas
preserve
* Solo etapas de TP y solo las dos variables de acceso/egreso
keep if inlist(vii_10_medio_transporte,4,8)
keep id_hogar id_persona id_viaje ///
     vii_28_cuadras_caminadas_hasta_e ///
     vii_30_cuadras_caminadas_al_baja

* Pasar a numérico de forma robusta
tempvar n28 n30
capture confirm numeric variable vii_28_cuadras_caminadas_hasta_e
if _rc {
    gen strL __s28 = trim(vii_28_cuadras_caminadas_hasta_e)
    replace __s28 = subinstr(__s28, ",", ".", .)
    gen double `n28' = real(__s28)
    drop __s28
}
else gen double `n28' = vii_28_cuadras_caminadas_hasta_e

capture confirm numeric variable vii_30_cuadras_caminadas_al_baja
if _rc {
    gen strL __s30 = trim(vii_30_cuadras_caminadas_al_baja)
    replace __s30 = subinstr(__s30, ",", ".", .)
    gen double `n30' = real(__s30)
    drop __s30
}
else gen double `n30' = vii_30_cuadras_caminadas_al_baja

replace `n28' = 0 if missing(`n28')
replace `n30' = 0 if missing(`n30')

* CLAVE: tomar el MÁXIMO por viaje (no sumar por cada etapa)
bys id_hogar id_persona id_viaje: egen double w2_max = max(`n28')
bys id_hogar id_persona id_viaje: egen double w3_max = max(`n30')
bys id_hogar id_persona id_viaje: keep if _n==1

* Bloques -> minutos (0.1 km por cuadra; 3.6 km/h)
gen double walk_min = (w2_max + w3_max) * (60*(0.1/4.5))

keep id_hogar id_persona id_viaje walk_min
tempfile walkTrip
save `walkTrip', replace
restore


* -- 4) RECONSTRUIR TRIP-LEVEL en packTrip y PEGAR a la base completa --
preserve
use `trips', clear                        // 1 fila por viaje
merge m:1 id_hogar id_persona id_viaje using `waitsTP', nogen keep(master match)
merge m:1 id_hogar id_persona id_viaje using `walkTrip', nogen keep(master match)

* rellenar faltantes (idéntico a tu lógica original)
replace W_min   = 0 if missing(W_min) | modo2==2     // Auto sin espera
replace walk_min = 0 if missing(walk_min)

* tiempo en vehículo puro
capture drop T_vehicle_min
gen double T_vehicle_min = dur_min - W_min - walk_min
replace T_vehicle_min = max(T_vehicle_min, 0)
* === Versión en horas (coherente con beta_T en 1/h y MATLAB) ===
gen double dur_h       = dur_min/60
gen double W_h         = W_min/60
gen double walk_h      = walk_min/60
gen double T_vehicle_h = T_vehicle_min/60

label var dur_h        "Duración total (h)"
label var W_h          "Espera TP (h)"
label var walk_h       "Caminar (h)"
label var T_vehicle_h  "Tiempo en vehículo (h)"


* empaquetar solo lo trip-level para pegar a la base "grande"
keep id_hogar id_persona id_viaje grupo periodo modo2 dur_min W_min walk_min T_vehicle_min dist_km
tempfile packTrip
save `packTrip', replace
restore

* ¡ahora sí! pegamos a tu base COMPLETA (no perdemos nada)
merge m:1 id_hogar id_persona id_viaje using `packTrip', nogen

* -- 5) T_en_vehículo = duración - espera - caminar (truncado a >=0) --
capture drop T_vehicle_min
gen double T_vehicle_min = dur_min - W_min - walk_min

count if T_vehicle_min < 0
if r(N)>0 di as err "Aviso: " r(N) " viajes con T_vehicle_min < 0; se truncarán a 0."
replace T_vehicle_min = max(dur_min - W_min - walk_min, 0)


* -- 6) Mostrar promedios (consola) --
preserve
keep if inlist(modo2,1,2) & !missing(T_vehicle_min)
collapse (mean) T_min=T_vehicle_min, by(modo2)
format T_min %9.2f
di as txt "=== T_en_vehículo promedio (min) por modo (trip-level) ==="
list modo2 T_min, noobs
restore

capture confirm numeric variable periodo
if !_rc {
    preserve
    keep if inlist(modo2,1,2) & !missing(periodo, T_vehicle_min)
    collapse (mean) T_min=T_vehicle_min, by(modo2 periodo)
    format T_min %9.2f
    di as txt "=== T_en_vehículo promedio (min) por modo y periodo (trip-level) ==="
    list modo2 periodo T_min, noobs
    restore

    preserve
    keep if modo2==1 & !missing(periodo, W_min)
    collapse (mean) W_min, by(periodo)
    format W_min %9.2f
    di as txt "=== Espera TP promedio (min) por periodo (trip-level) ==="
    list periodo W_min, noobs
    restore
}

*========================================================
* (D) L (km por viaje) — SOLO mostrar
*========================================================
* (D) L (km por viaje) — DEFINICIÓN CANÓNICA (trip-level)
*========================================================
* OJO: medimos L en la base trip-level para evitar sobreponderar por etapas.

preserve
    use `trips', clear                      // 1 fila por viaje (ya creada arriba)
    keep if inlist(modo2,1,2) & inlist(periodo,1,2)
    quietly summ dist_km if !missing(dist_km) & dist_km>0, meanonly
    di as res "L_trip (modelo, sin rangos) = " %6.3f r(mean)
restore
* === Guardar tabla de viajes ENRIQUECIDA (1 fila por viaje) ===
* === Guardar tabla de viajes ENRIQUECIDA (1 fila por viaje) ===
preserve
    * Asegurar versión en horas ANTES de filtrar variables
    capture drop dur_h W_h walk_h T_vehicle_h
    gen double dur_h       = dur_min/60
    gen double W_h         = W_min/60
    gen double walk_h      = walk_min/60
    gen double T_vehicle_h = T_vehicle_min/60

    keep id_hogar id_persona id_viaje ///
         grupo periodo modo2 ///
         dur_min W_min walk_min T_vehicle_min dist_km ///
         dur_h   W_h   walk_h   T_vehicle_h
    bys id_hogar id_persona id_viaje: keep if _n==1
    save "outputs/EMD18_trip_enriched.dta", replace
    export delimited using "outputs/EMD18_trip_enriched.csv", replace
restore
di as res ">> OK: outputs/EMD18_trip_enriched.dta (mins + horas) creado."
*==================== (E7R) Velocidad simple con base refinada ====================*
clear all
set more off
cd "C:\Users\drake\OneDrive\Escritorio\tesis\Nueva carpeta"

use "outputs/EMD18_trip_enriched.dta", clear

* Chequeos mínimos
foreach v in dist_km dur_min periodo modo2 {
    capture confirm numeric variable `v'
    if _rc {
        di as err "Falta `v' en EMD18_trip_enriched.dta"; exit 605
    }
}


* Velocidad puerta-a-puerta (km/h) usando SOLO distancia y duración total (min→h)
capture drop kmh_simple
gen double kmh_simple = dist_km/(dur_min/60) if dur_min>0 & !missing(dist_km)
describe kmh_simple
* Nos quedamos con los modos/periodos del modelo
keep if inlist(modo2,1,2) & inlist(periodo,1,2)

* Tabla: promedio y mediana por periodo × modo
preserve
    collapse (mean)   mean_kmh = kmh_simple ///
             (median) med_kmh  = kmh_simple ///
             (count)  N        = kmh_simple, by(periodo modo2)

    label define P 1 "peak" 2 "off", replace
    label define M 1 "TP"   2 "Auto", replace
    label values periodo P
    label values modo2 M

    gen str5 periodo_lbl = cond(periodo==1,"peak","off")
    gen str3 modo_lbl    = cond(modo2==1,"TP","Auto")
    order periodo_lbl modo_lbl mean_kmh med_kmh N

    di as res "=== Velocidad (km/h) por periodo × modo — media y mediana ==="
    list, noobs

    export delimited using "outputs/VELsimple_refined_byperiod_mode_mean_median.csv", replace
    save "outputs/VELsimple_refined_byperiod_mode_mean_median.dta", replace
restore
*==================== FIN (E7R) ===================================================*
use "outputs/EMD18_trip_enriched.dta", clear

* 1) Una fila por viaje (evita sobreponderar viajes con más etapas)
bys id_hogar id_persona id_viaje: keep if _n==1

* 2) Mismas celdas del modelo
keep if inlist(modo2,1,2) & inlist(periodo,1,2)



* 4) L global = media de dist_km en esta muestra "canónica"
quietly summ dist_km, meanonly
local L_global = r(mean)
di as res "L_global (km por viaje, media, trip-level, modelo+filtros) = " %6.3f `L_global'
