// Created: 2020-01-06 15:13:31 UTC
[GLOBAL]

//-- #include "global.h"
//--  Y = (max(0,x).^k)./(tau^k + max(0,x).^k);

double kmEQ(double c, double k, double vmax) {
  double a = vmax * c;
  return a/(k + c);
}

double HillEQ(double x, double k, double tau) {
  double a = pow(std::max(x,0.0),k);
  return a/(pow(tau,k) + a);
}

[PARAM] 
// PK parameters
TVVirn    = 0.224      // Volume central compartment
TVVsn     = 0.0349      // Volume SN-38
TVVmax  = 0.78   // Max metabolism rate
TVKm    = 0.285    // Amount of half Vmax
TVk0i   = 10.56      // elimination rate constant
TVk0s   = 53.76      // elimination rate constant
TVk12   = 24*0.169
TVk21   = 24*0.121

// effect parameters
TVec50 = 0.1 //  level of compound to have half effect
gamma = 1

// PD parameters
TVkmax = 0.419 // max kill rate for the compound
TVkg   = 0.366 // net proliferation rate
TVkdq  = 0.0154 // destruction rate of quiescent cells
TVkneg = 0.00909 // retardation of growth constant in Gompertz model, microenvironment constant
TVkdn  = 0.0349

// modelled initial values
TVV0 = 0.27

TVkdummy = 0 // for sensitivity analysis

[OMEGA] 0 0 0 0 0 0 0 0 0 0 0 0 0

[SIGMA] @labels PROP ADD
0 0

[INIT]
// Initial conditions (5)
IRN_C = 0      // specie 1
IRN_P = 0
SN38_C = 0     // specie 2
Q = 0
P = 0
N = 0
A = 0
AUC = 0

[MAIN]
double Virn     = TVVirn*exp(ETA(1));
double Vsn      = TVVsn*exp(ETA(2));
double Vmax   = TVVmax*exp(ETA(3));
double Km     = TVKm*exp(ETA(4));
double k0i    = TVk0i*exp(ETA(5));
double k0s    = TVk0s*exp(ETA(6));
double k12    = TVk12*exp(ETA(7));
double k21    = TVk21*exp(ETA(8));
double ec50   = TVec50*exp(ETA(9));
double kmax   = TVkmax*exp(ETA(10));
double kg     = TVkg*exp(ETA(11));
double kdq    = TVkdq*exp(ETA(11));
double kneg   = TVkneg*exp(ETA(12));
double kdn    = TVkdn*exp(ETA(13));
double P0     = TVV0;
double Q0     = 0;

[ODE]

double IRN_c    = IRN_C/ Virn;
double Fm       = kmEQ(IRN_c, Km, Vmax);
double SN38_cp  = SN38_C / Vsn;
double effects  = HillEQ(SN38_cp, gamma, ec50);
double m = kg + kdq;
double P_mod = P + P0;   /// modelling initial conditions to match Stan
double Q_mod = Q + Q0;
double Vt = P_mod + Q_mod + N + A;
double PVratio = P_mod/Vt;
double PVratioinf = kdq / m;
double Psi = -m*PVratio*PVratio + (m + kneg)*PVratio - kneg*PVratioinf;

// PK system
dxdt_IRN_C    =  - Fm - IRN_C * k0i + k21*IRN_P - k12*IRN_C;
dxdt_IRN_P    =  - k21*IRN_P + k12*IRN_C; 
dxdt_SN38_C   =  Fm - SN38_C * k0s;

// PD model
dxdt_P = kg * P_mod - Vt * Psi - kmax*effects*P_mod;
dxdt_Q = Vt * Psi - kdq * Q_mod;
dxdt_N = kdq * Q_mod - kdn * N;
dxdt_A = kmax*effects*P_mod - kdn * A;

dxdt_AUC = Vt;

[TABLE] 
double IPRED = Vt;
double TPratio = P/Vt;
double DV = IPRED*(1+PROP)+ADD;
int i = 0;
while(DV < 0 && i < 100){
  simeps();
  DV = IPRED*(1+PROP)+ADD;
  ++i;
}

$CAPTURE TPratio Vt effects DV P_mod Q_mod
  
  