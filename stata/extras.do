*==================== Velocidades: door-to-door vs en-vehículo ====================*
clear all
set more off
cd "C:\Users\drake\OneDrive\Escritorio\tesis\Nueva carpeta"

use "outputs/EMD18_trip_enriched.dta", clear

* 0) Una fila por viaje y celdas del modelo
bys id_hogar id_persona id_viaje: keep if _n==1
keep if inlist(modo2,1,2) & inlist(periodo,1,2)

* 1) Filtros robustos (evitar extremos)
keep if inrange(dur_min, 5, 120) & inrange(dist_km, 0.5, 180)

* 2) Definiciones de velocidad
capture drop kmh_door kmh_inveh
gen double kmh_door = dist_km / (dur_min/60)

capture drop T_vehicle_min
gen double T_vehicle_min = max(dur_min - W_min - walk_min, 1)   // evitar 0
gen double kmh_inveh     = dist_km / (T_vehicle_min/60)

label var kmh_door  "Velocidad puerta-a-puerta (km/h)"
label var kmh_inveh "Velocidad en vehículo (km/h)"

* 3) Resumen robusto: mediana y rango intercuartílico por periodo × modo
preserve
    keep if !missing(kmh_door, kmh_inveh, periodo, modo2)

    bys periodo modo2: egen double p25_door = pctile(kmh_door),  p(25)
    bys periodo modo2: egen double p50_door = pctile(kmh_door),  p(50)
    bys periodo modo2: egen double p75_door = pctile(kmh_door),  p(75)

    bys periodo modo2: egen double p25_in   = pctile(kmh_inveh), p(25)
    bys periodo modo2: egen double p50_in   = pctile(kmh_inveh), p(50)
    bys periodo modo2: egen double p75_in   = pctile(kmh_inveh), p(75)

    bys periodo modo2: egen long   N = count(kmh_door)
    gen double IQR_door = p75_door - p25_door
    gen double IQR_in   = p75_in   - p25_in

    bys periodo modo2: keep if _n==1
    label define P 1 "peak" 2 "off", replace
    label define M 1 "TP"   2 "Auto", replace
    label values periodo P
    label values modo2 M

    order periodo modo2 p50_door p25_door p75_door IQR_door p50_in p25_in p75_in IQR_in N
    rename (p50_door p25_door p75_door p50_in p25_in p75_in) ///
           (med_door  q25_door q75_door med_in  q25_in  q75_in)

    di as res "=== Velocidades (km/h) por periodo × modo ==="
    di as res "Puerta-a-puerta: mediana [q25–q75]  |  En-vehículo: mediana [q25–q75]"
    list periodo modo2 med_door q25_door q75_door IQR_door med_in q25_in q75_in IQR_in N, noobs

    export delimited using "outputs/VEL_door_vs_inveh_byperiod_mode.csv", replace
    save "outputs/VEL_door_vs_inveh_byperiod_mode.dta", replace
restore
*==============================================================================*

*==================== Velocidad por GRUPO con L global vs L_g ====================*
clear all
set more off
cd "C:\Users\drake\OneDrive\Escritorio\tesis\Nueva carpeta"

use "outputs/EMD18_trip_enriched.dta", clear
bys id_hogar id_persona id_viaje: keep if _n==1
keep if inlist(modo2,1,2) & inlist(periodo,1,2)

* Filtros robustos (ajusta si quieres)
keep if inrange(dur_min,5,120) & inrange(dist_km,0.5,80)

* Asegurar T_vehicle_min
capture confirm variable T_vehicle_min
if _rc {
    replace W_min   = 0 if missing(W_min)   | modo2==2
    replace walk_min= 0 if missing(walk_min)
    gen double T_vehicle_min = max(dur_min - W_min - walk_min, 1)
}
else replace T_vehicle_min = max(T_vehicle_min, 1)

* L global (media de dist_km)
quietly summ dist_km
local L_global = r(mean)

preserve
    * Medias por grupo
    collapse (mean) dur_min T_vehicle_min dist_km, by(grupo)
    format dur_min T_vehicle_min dist_km %9.2f

    * Velocidades ratio-de-medias con L global y con L_g
    gen double vdoor_Lglob  = `L_global' / (dur_min/60)
    gen double vinveh_Lglob = `L_global' / (T_vehicle_min/60)
    gen double vdoor_Lg     = dist_km     / (dur_min/60)
    gen double vinveh_Lg    = dist_km     / (T_vehicle_min/60)

    label var vdoor_Lglob  "km/h door (L global)"
    label var vinveh_Lglob "km/h in-veh (L global)"
    label var vdoor_Lg     "km/h door (L_g)"
    label var vinveh_Lg    "km/h in-veh (L_g)"

    order grupo dist_km dur_min T_vehicle_min vdoor_Lglob vinveh_Lglob vdoor_Lg vinveh_Lg
    format vdoor_* vinveh_* %9.2f

    di as res "=== Referencia: L_global (km por viaje, media) ==="
    di %6.3f `L_global'
    di as res "=== Velocidades por GRUPO (ratio de medias) ==="
    list, noobs sepby(grupo)
restore

* (Opcional) chequeo de mediana de velocidades por grupo (para comparar)
preserve
    gen double kmh_door_trip  = dist_km/(dur_min/60)
    gen double kmh_inveh_trip = dist_km/(T_vehicle_min/60)
    collapse (median) med_door=kmh_door_trip (median) med_in=kmh_inveh_trip, by(grupo)
    format med_* %9.2f
    di as res "=== Mediana de velocidades por GRUPO (solo referencia) ==="
    list, noobs
restore
*==============================================================================*
*========================================================
* SCRIPT: vel_prom_modo_periodo.do
* Objetivo: Velocidad promedio por modo × periodo (peak/off)
*           v = L_global / mean(dur_h)  [y análogo con T_vehicle_h si existe]
* No guarda archivos; solo imprime tablas.
*========================================================
clear all
set more off
cd "C:\Users\drake\OneDrive\Escritorio\tesis\Nueva carpeta"

use "outputs/EMD18_trip_enriched.dta", clear

* --- Trip-level canónico: 1 fila por viaje ---
bys id_hogar id_persona id_viaje: keep if _n==1

* --- Celdas del modelo y requisitos mínimos ---
keep if inlist(modo2,1,2) & inlist(periodo,1,2)
keep if !missing(dist_km, dur_min)

* --- L_global (km por viaje, media trip-level) ---
quietly summ dist_km if dist_km>0, meanonly
local L_global = r(mean)

* --- Duración promedio por modo × periodo (min → h) ---
preserve
    collapse (mean) dur_min = dur_min, by(modo2 periodo)
    gen double dur_h = dur_min/60

    * Velocidad puerta-a-puerta (km/h): ratio-de-medias
    gen double v_door = `L_global' / dur_h

    * Etiquetas legibles
    label define P 1 "peak" 2 "off", replace
    label define M 1 "TP"   2 "Auto", replace
    label values periodo P
    label values modo2 M

    order periodo modo2 dur_min dur_h v_door
    format dur_min dur_h v_door %9.2f

    di as res "L_global (km) = " %6.3f `L_global'
    di as res "=== Velocidad promedio (km/h) — puerta-a-puerta ==="
    list periodo modo2 dur_min dur_h v_door, noobs
restore

* --- (Opcional rápido) En-vehículo si ya está en la base ---
capture confirm variable T_vehicle_min
if !_rc {
    preserve
        keep if !missing(T_vehicle_min)
        collapse (mean) Tv_min = T_vehicle_min, by(modo2 periodo)
        gen double Tv_h = Tv_min/60
        gen double v_inveh = `L_global' / Tv_h

        label define P 1 "peak" 2 "off", replace
        label define M 1 "TP"   2 "Auto", replace
        label values periodo P
        label values modo2 M

        order periodo modo2 Tv_min Tv_h v_inveh
        format Tv_min Tv_h v_inveh %9.2f

        di as res "=== Velocidad promedio (km/h) — en vehículo ==="
        list periodo modo2 Tv_min Tv_h v_inveh, noobs
    restore
}


