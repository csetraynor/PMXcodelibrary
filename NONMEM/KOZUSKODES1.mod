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
$SUB ADVAN13 TOL=12
$MODEL
	COMP = (IRNC)
	COMP = (IRNP)
	COMP = (SN38C)
	COMP = (P)
	COMP = (Q)
	COMP = (N)
	COMP = (A)

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

TVKDQ=EXP(THETA(4))
MU_4=LOG(TVKDQ)
KDQ = EXP(MU_4+ETA(4))

TVKNEG=EXP(THETA(5))
MU_5=LOG(TVKNEG)
KNEG=EXP(MU_5+ETA(5))

TVKDN=EXP(THETA(6))
MU_6=LOG(TVKDN)
KDN=EXP(MU_6+ETA(6))

M = KG + KDQ
PVINF = KDQ / M

$DES

;;variables
IRNCP = A(1)/VIRN
Fm = IRNCP/ (KM+IRNCP)
SN38CP = A(3)/VSN38
Pmod = A(4)+P0
Vt = Pmod + A(5) + A(6) + A(7)
PVratio = Pmod / Vt
Psi = -m*PVratio*PVratio + (m+kneg)*PVratio - kneg*PVINF

;;PKMODEL
DADT(1)  =  - Fm - A(1) * k0i + k21*A(2) - k12*A(1)
DADT(2)  =  - k21*A(2) + k12*A(1)
DADT(3)  =   Fm - SN38_C * k0s
;;PDMODEL
DADT(4) = kg * Pmod - Vt * Psi - kmax*effects*Pmod
DADT(5) = Vt * Psi - kdq * A(5)
DADT(6) = kdq * A(5) - kdn * A(6)
DADT(7) = kmax*effects*Pmod - kdn * A(7)

$ERROR 

IPRED=Vt
W=SQRT(SIGMA(1,1)*IPRED**2+SIGMA(2,2))  ; proportional + additive error
IRES=DV-IPRED
IWRES=IRES/W
Y=IPRED+W*EPS(1)


$THETA
-1             	; KG ; days-1 ; LOG
-2.5            ; KMAX  ; days-1 ; LOG
-0.5          	; P0  ; ml ; LOG
-1		; KDQ ; days-1 ; LOG
-1		; KNEG ; days-1 ; LOG
-1		; KDN ; days-1 ; LOG


$OMEGA
0.1			; IIV_KG ; LOG
0.1			; IIV_KMAX ; LOG
0.1			; IIV_P0 ; LOG
0.1			; IIV_KDQ ; LOG
0.1			; IIV_KNEG ; LOG
0.1			; IIV_KN ; LOG


$SIGMA
0.1		; prop error
......		; add error

; Parameter estimation - FOCE
;$EST METHOD=1 INTER NOABORT MAXEVAL=9999 PRINT=1 NSIG=3 SIGL=9

; Parameter estimation - IMP
$EST METHOD=IMP ISAMPLE=300 NITER=300 RANMETHOD=3S2
CTYPE=3 CITER=10 CALPHA=0.05 CINTERVAL=3
PRINT=1 NOABORT INTERACTION

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
