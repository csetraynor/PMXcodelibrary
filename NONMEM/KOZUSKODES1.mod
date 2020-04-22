;; 1. Based on: 
;; 2. Description: 3CMT+2CMT PK+PD
;; x1. Author: kllt565

$PROBLEM TWO-SENS-GROW

;; 4. Date: 23.04.2020
;; 5. Version: 1
;; 6. Label:
;; Basic model
;; 7. Structural model:
;; Five compartmental PKPD
;; 8. Covariate model:
;; No covariates
;; 9. Inter-individual variability:
;; KG KDQ 
;; 10. Inter-occasion variability:
;; No IOV
;; 11. Residual variability:
;; Additive + Proportional
;; 12. Estimation:
;; IMP

$INPUT ID AMT TIME DV WT
$DATA ../DerivedData/COLO205.csv IGNORE=@ 
$SUB ADVAN2

$PK

TVKG=EXP(THETA(1))
MU_1=LOG(TVKG)
KG = EXP(MU_1+ETA(1))

TVKMAX=EXP(THETA(2))
MU_2=LOG(TVKMAX)
KMAX = EXP(MU_2+ETA(2))

TVP0=EXP(THETA(3))
MU_3=LOG(TVP0)
P0 = EXP(MU_3+ETA(3))

TV

TVec50 = 0.1 //  level of compound to have half effect
gamma = 1

// PD parameters
TVkg   = 0.366 // net proliferation rate

TVkmax = 0.419 // max kill rate for the compound

TVkdq  = 0.0154 // destruction rate of quiescent cells
TVkneg = 0.00909 // retardation of growth constant in Gompertz model, microenvironment constant
TVkdn  = 0.0349

CL = K*V

S2 = V

$ERROR 

IPRED=F
W=SQRT(THETA(4)**2+THETA(5)**2*IPRED*IPRED)  ; proportional + additive error
IRES=DV-IPRED
IWRES=IRES/W
Y=IPRED+W*EPS(1)


$THETA
1             	; KA ; h-1 ; LOG
-2.5            ; K  ; h-1 ; LOG
-0.5          	; V  ; L ; LOG
0.1           	; add error
0.1           	; prop error

$OMEGA
0.1			; IIV_KA ; LOG
0.1			; IIV_K ; LOG
0.1			; IIV_V ; LOG


$SIGMA
1 FIX

; Parameter estimation - FOCE
;$EST METHOD=1 INTER NOABORT MAXEVAL=9999 PRINT=1 NSIG=3 SIGL=9

; Parameter estimation - IMP
$EST METHOD=IMP ISAMPLE=300 NITER=300 RANMETHOD=3S2
CTYPE=3 CITER=10 CALPHA=0.05 CINTERVAL=3
PRINT=1 NOABORT INTERACTION GRD=DDDSS

; Parameer estimation - SAEM
;$EST METHOD=SAEM ISAMPLE=2 NBURN=1000 NITER=500 RANMETHOD=3S2
;CTYPE=3 CITER=10 CALPHA=0.05 CINTERVAL=10
;PRINT=1 NOABORT INTERACTION GRD=DDDSS

; Objective function and covariance evaluation
$EST METHOD=IMP INTER EONLY= 1 MAPITER=0 ISAMPLE = 2000 NITER = 10 RANMETHOD=3S2 NOABORT PRINT=1 NSIG=3 SIGL=9 GRD=DDDSS

$COV PRINT=E UNCONDITIONAL SIGL=10

$TABLE ID TIME IPRED IWRES IRES CWRES NPDE
FILE=sdtab1 NOPRINT ONEHEADER FORMAT=tF13.4
$TABLE ID ETAS(1:LAST); individual parameters
FILE=patab1 NOPRINT ONEHEADER FORMAT=tF13.4
$TABLE ID ; continuous covariates
FILE=cotab1 NOPRINT ONEHEADER FORMAT=tF13.4
$TABLE ID ; categorical covariates
FILE=catab1 NOPRINT ONEHEADER FORMAT=tF13.4
