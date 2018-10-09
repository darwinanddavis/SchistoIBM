; script-file DEB-INF-IBM
; Authors: Dave Civitello (dcivite@emory.edu), adapted from Ben Martin's DEB-IBM

; implementation of the DEB-infection equationsof Civitello et al. XXXX in an IBM
; check ... for the user-manual and the ODD
; published in "...", ... 20XX
; ==========================================================================================================================================
; ========================== DEFINITION OF PARAMETERS AND STATE VARIABLES ==================================================================
; ==========================================================================================================================================
breed [snails snail]
breed [eggs egg]


; global parameters: are accessible for patches and turtles
globals[
  ;e0    ; t L^2, initial reserves of the embryos at the start of the simulation
  L_0      ; cm, initial structural volume

]
; ------------------------------------------------------------------------------------------------------------------------------------------
; parameters for the environment: here only prey density

patches-own[
  F        ; # / m^2, prey density
  M        ; # / m^2 miracidial density
  Z        ; # / m^2 cercarial density
  d_F      ; change of prey density in time
  d_M      ; change in miracidial density in time
  d_Z      ; change in cercarial density in time
]
; ------------------------------------------------------------------------------------------------------------------------------------------

; definition of parameters for the individuals:
; the notation follows the DEBtool-notation as far as possible
; deviation: rates are indicated with "_rate" rather than a dot
; each individual(turtle) in the model has the following parameters
snails-own[
  ; - - - - - - - - - - - - - - - STATE VARIABLES - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ; --- Host state variables ---
  L           ; cm, structural length
  dL          ; change of structural length in time
  e_H         ; scaled reserve density of hosts
  de_H        ; change of scaled reserve density of hosts in time
  D           ; host development
  dD          ; change of host development in time
  R_H         ; energy in reproduction buffer
  dR_H        ; change of energy in reproduction buffer (reproduction rate)
      ; --- Parasite state variables ---
  P           ; parasite biomass
  dP          ; change of parasite biomass in time
  R_P         ; energy in parasite reproduction
  dR_P        ; change of energy in parasite reproduction

  ; - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ; - - - - - - - - - - - - - - - FLUXES and COMPOUND PARAMETERS (used by several submodels) - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  a_M         ; assimilation flux
  C           ; mobilization flux
  g           ; costs of growth
  f_H         ; host scaled functional response
  m_v         ; somatic maintenance rate
  m_D         ; development maintenance rate
  Pdens       ; parasite biomass density
  repro_P     ; parasite biomass invested in reproduction
  f_P         ; parasite scaled functional response
  kstar       ; realized kappa
  DELTA       ; mobilized reserve shortfall
  ingestion   ; scaled ingestion rate
  biomass     ; allometric estimate of biomass
  d_EX        ; exposure to miracidia
  d_EGGS      ; number of eggs produced
  d_CERCS     ; number of cercariae released
  ; - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ; - - - - - - - - - - - - - - -  DEB PARAMETERS (with dimension and name) - - - - - - - - - - - - - - - - - - - - - - - - - - - -
]

eggs-own[
  ; - - - - - - - - - - - - - - - STATE VARIABLES - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  age         ; age of egg, hatches at 7 d
  ; - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ; - - - - - - - - - - - - - - -  DEB PARAMETERS (with dimension and name) - - - - - - - - - - - - - - - - - - - - - - - - - - - -
]
; ==========================================================================================================================================
; ========================== SETUP PROCEDURE: SETTING INITIAL CONDITIONS ===================================================================
; ==========================================================================================================================================

to setup
  clear-all

 create-snails 60                   ; 10 turtles are created in the beginning
 ask  snails  [
  set L random-float 2 + 4
  set e_H 0.9
  set D 0
  set R_H 0
  set P 0
  set R_P 0
    ;individual-variability  ; first their individual variability in the parameter is set
  ;calc-embryo-reserve-investment     ; then the initial energy is calculated for each
 ]

 ask patches [ set F 1 ]; set initial value of prey to their carrying capacity
  reset-ticks
end

; ==========================================================================================================================================
; ========================== GO PROCEDURE: RUNNING THE MODEL ===============================================================================
; ==========================================================================================================================================
; the go statement below is the order in which all procedures are run each timestep

to go
  ask snails[calc-DEB]                       ; first all individuals calculate the change in their state variables based on the current conditions
  ask eggs[calc-age]                         ; increment age of all developing eggs
  ask patches [calc-env]        ; calculate changes in food density

   update                           ; the the state variables of the individuals and prey are updated based on the delta value

   ask snails with [R_H >= 0.015] ; how many eggs will each host lay
   [lay-eggs]

  tick
  do-plots                          ; then the plots are updated
  if count turtles = 0 [stop]

end

; ==========================================================================================================================================
; ========================== SUBMODELS =====================================================================================================
; ==========================================================================================================================================



; ------------------------------------------------------------------------------------------------------------------------------------------
; ----------------- DEB DYNAMICS -------------------------------------------------------------------------------------------------------
; ------------------------------------------------------------------------------------------------------------------------------------------
; Calculates changes in DEB state variables for a host feeding on a logistically growing food source


to calc-DEB
  ifelse D <= D_B
    [set f_H 0]
    [set f_H F / (F_50 + F)]
  ; set compound parameters to keep calculating code clean
  set a_M i_M * y_EF
  set kstar min list (k + alpha * P) 1
  set g 1 / (kstar * y_VE * E_M)
  set Pdens P / (Chi * L ^ 3)
  set repro_P Pdens ^ 2 / (p_h ^ 2 + Pdens ^ 2)
  set m_V k * a_M / (Chi * L_M)
  set m_D mu_D * m_V
  set C (g * e_H)/( g + e_H ) * (a_M * L ^ 2 + m_V * Chi * L ^ 3 / (kstar * g))
  set f_P e_H / (e_50 + e_H)
  set DELTA m_V * Chi * L ^ 3 - kstar * C
  set ingestion i_M * f_H * L ^ 2
  set biomass 0.0096 * L ^ 3
  ; calculate rates of change for state variables
  set dL max list (g * y_VE / ( 3 * Chi) * (kstar * a_M * e_H - m_V * Chi * L)/(g + e_H)) 0
  ;set dL max list (y_VE /(3 * Chi) * (kstar * a_M * e - m_V * Chi * L) / (1 + y_VE * kstar * E_M * e)) 0
  set de_H a_M / (Chi * E_M * L) * (f_H - e_H) - i_PM * f_P * Pdens / E_M
  ifelse D < D_R ; maturity and repro for non- and mild starvation
    [set dD (1 - kstar) * C - m_D * D - (min list DELTA 0) ]
    [set dD min list ((1 - kstar) * C - m_D * D - (min list DELTA 0)) 0]
  ifelse D < D_R ; they only invest into maturity until they reach puberty
    [set dR_H 0]
    [set dR_H max list ((1 - kstar) * C - m_D * D - (min list DELTA 0)) 0]
  if m_V * Chi * L ^ 3 - C > 0 ; if extreme starvation
    [set dL 0
     set de_H (a_M * L ^ 2 * f_H - (m_V * Chi * L ^ 3 + i_PM * f_P * P))/(E_M * Chi * L ^ 3)
     set dD (- m_D * D)
     set dR_H 0]
  set dP (Y_PE * i_PM * f_P - m_P - repro_P) * P
  set dR_P y_RP * repro_P * P
  set d_EX random-poisson (epsilon * M * (1 / volume)) ; number of parasites depleted from env
  set P P + 0.00001425 * random-poisson (susceptibility * d_EX)
  set d_EGGS min list 120 random-poisson (R_H / 0.015)
  set d_CERCS random-poisson (R_P / 0.00004)
  set dR_P dR_P - d_CERCS * 0.00004
  if e_H <= 0 [die]
end


; ------------------------------------------------------------------------------------------------------------------------------------------
; ----------------- Egg development -------------------------------------------------------------------------------------------------------
; ------------------------------------------------------------------------------------------------------------------------------------------
; increments age of developing eggs


to calc-age
  set age age + 1 / timestep
  if age >= 10 [
    if random-float 1 > 0.75 [die]
    hatch-snails 1
      [set L 0.75
      set e_H 0.9
      set D 0
      set R_H 0
      set P  0
      set R_P 0]
    die
  ]
end
;-------------------------------------------------------------------------------------------------------------------------------------------
;--------- LAY EGGS ------------------------------------------------------------------------------------------------------------------------
;-------------------------------------------------------------------------------------------------------------------------------------------
;the following procedure is run for mature individuals which have enough energy to reproduce
; they create 1 offspring and give it the following state variables and DEB parameters
;the initial reserves is set to the value determined by the bisection method in "calc_egg_size"

to lay-eggs
  set R_H R_H - d_EGGS * 0.015
  hatch-eggs d_EGGS
    [set age 0]
end


; ------------------------------------------------------------------------------------------------------------------------------------------
; ----------------- LOGISTIC PREY ----------------------------------------------------------------------------------------------------------
; ------------------------------------------------------------------------------------------------------------------------------------------
 ;the following procedure calculates change in prey density this procedure is only run when prey dynamics are set to "logistic" in the user interface

to calc-env
   set d_F (r_F) * F * (1 - (F / K_F))   - sum [ ingestion ] of snails-here / volume ; units may be off for volume
  set d_M (- mu_M * M) - sum [d_EX] of snails-here / volume; background loss rate; depletion through exposure needs added
   set d_Z (- mu_Z * Z) + sum [d_CERCS] of snails-here / volume; background death rate + production from infected hosts
end


; ------------------------------------------------------------------------------------------------------------------------------------------
; ----------------- UPDATE -----------------------------------------------------------------------------------------------------------------
; ------------------------------------------------------------------------------------------------------------------------------------------

to update
; individuals update their state variables based on the calc_state variable proccesses
  ask snails
  [
    set L L + dL    / timestep
    set e_H e_H + de_H    / timestep
    set D D + dD    / timestep
    set R_H R_H + dR_H    / timestep
    set P P + dP    / timestep
    set R_P R_P + dR_P    / timestep
   ;if ticks mod timestep = age-day [if random-float 1 < background-mortality [die]]
 ]
  ask patches [
    set F F + d_F / timestep
    set M (M + d_M / timestep)  + ifelse-value member? (ticks / timestep) [1 15 29 43] [2500 / volume][0]
    set Z Z + d_Z / timestep
  ]
end

; ------------------------------------------------------------------------------------------------------------------------------------------
; ----------------- PLOT -------------------------------------------------------------------------------------------------------------------
; ------------------------------------------------------------------------------------------------------------------------------------------

to do-plots
 set-current-plot "stage class density"
 set-current-plot-pen "embryo"
 set-plot-pen-interval 1 / timestep
    ifelse any? eggs  [plot count eggs / volume]
    [plot 0]
   set-current-plot-pen "juvenile"
   set-plot-pen-interval 1 / timestep
  ifelse any? snails with [D  < D_R] [plot count snails with [D  < D_R] / volume]
  [plot 0]
  set-current-plot-pen "adult"
  set-plot-pen-interval 1 / timestep
  ifelse any? snails with [D  > D_R] [plot count snails with [D  > D_R] / volume]
  [plot 0]




  set-current-plot "Environment"
  set-current-plot-pen "food density"
  set-plot-pen-interval 1 / timestep
  plot mean [F] of patches
  set-current-plot-pen "miracidia"
  set-plot-pen-interval 1 / timestep
  plot mean [M] of patches
    set-current-plot-pen "cercariae"
  set-plot-pen-interval 1 / timestep
  ifelse (any? snails with [d_CERCS > 0]) [plot  mean [d_CERCS] of snails with [d_CERCS > 0]] [plot 0] ; [Z] of patches
  ;plot  sum [Z] of patches
;
;  set-current-plot "population density"
;  set-plot-pen-interval 1 / timestep
; plot count snails with [D > 0]



  set-current-plot "size distribution"
  histogram [L] of snails with [D > D_B]

    set-current-plot "juv e distribution"
  histogram [e_H] of snails with [D  < D_R]

    set-current-plot "adult e distribution"
  histogram [e_H] of snails with [D  > D_R]

    set-current-plot "parasite distribution"
  histogram [d_CERCS * timestep] of snails with [d_CERCS > 0]
end
@#$#@#$#@
GRAPHICS-WINDOW
-2
11
31
45
-1
-1
25.0
1
10
1
1
1
0
1
1
1
0
0
0
0
1
1
1
ticks
30.0

BUTTON
24
21
87
54
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
25
58
88
91
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

INPUTBOX
148
98
198
158
L_M
34.0
1
0
Number

INPUTBOX
208
98
258
158
Chi
8.2E-5
1
0
Number

INPUTBOX
321
98
371
158
y_EF
0.16
1
0
Number

INPUTBOX
12
367
62
427
i_M
0.027
1
0
Number

INPUTBOX
67
366
117
426
F_50
0.5
1
0
Number

INPUTBOX
266
98
316
158
E_M
74.2
1
0
Number

INPUTBOX
15
141
65
201
k
0.816
1
0
Number

INPUTBOX
16
206
66
266
mu_D
0.02
1
0
Number

INPUTBOX
73
206
123
266
D_R
0.385
1
0
Number

INPUTBOX
71
141
121
201
D_B
0.0
1
0
Number

INPUTBOX
238
180
288
240
ALPHA
0.1
1
0
Number

INPUTBOX
184
244
234
304
p_h
9.92
1
0
Number

INPUTBOX
239
244
289
304
e_50
0.04
1
0
Number

INPUTBOX
183
181
233
241
i_PM
0.316
1
0
Number

INPUTBOX
183
382
233
442
y_PE
0.99
1
0
Number

INPUTBOX
237
317
287
377
m_P
0.006
1
0
Number

INPUTBOX
181
314
231
374
y_RP
0.024
1
0
Number

INPUTBOX
10
300
60
360
r_F
0.9
1
0
Number

INPUTBOX
67
301
117
361
K_F
5.0
1
0
Number

INPUTBOX
218
10
268
70
Volume
500.0
1
0
Number

INPUTBOX
155
10
211
70
timestep
20.0
1
0
Number

INPUTBOX
376
98
426
158
y_VE
0.8
1
0
Number

PLOT
649
11
1278
161
stage class density
time
density
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"embryo" 1.0 0 -7500403 true "" ""
"juvenile" 1.0 0 -14439633 true "" ""
"adult" 1.0 0 -7858858 true "" ""

PLOT
651
334
851
484
size distribution
size
frequency
0.0
33.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 1 -7500403 true "" ""

PLOT
863
332
1063
482
juv e distribution
NIL
NIL
0.0
1.0
0.0
10.0
true
false
"" ""
PENS
"default" 0.1 1 -16777216 true "" ""

PLOT
1074
335
1274
485
adult e distribution
NIL
NIL
0.0
1.0
0.0
10.0
true
false
"" ""
PENS
"default" 0.1 1 -16777216 true "" ""

PLOT
649
168
1301
318
Environment
NIL
NIL
0.0
10.0
0.0
1.0
true
true
"" ""
PENS
"food density" 1.0 0 -16777216 true "" ""
"miracidia" 1.0 0 -7500403 true "" ""
"cercariae" 1.0 0 -2674135 true "" ""

INPUTBOX
14
435
64
495
mu_M
1.0
1
0
Number

INPUTBOX
238
383
288
443
epsilon
17.0
1
0
Number

INPUTBOX
71
435
121
495
susceptibility
0.05
1
0
Number

INPUTBOX
127
442
177
502
mu_Z
1.0
1
0
Number

PLOT
444
333
644
483
parasite distribution
NIL
NIL
0.0
10000.0
0.0
10.0
true
false
"" ""
PENS
"default" 0.1 1 -16777216 true "" ""

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.0.4
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
