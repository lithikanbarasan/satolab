// Juan G Restrepo Model with ICaL coupled gating

#ifndef ___SUBCELL_H
#define ___SUBCELL_H

#include <iostream>
using namespace std;
#define _USE_MATH_DEFINES
#include <cfloat>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <fstream>

class CSubcell {
 private:
  // cell parameters
  int nxny;

  int finemesh3;
  int nnx;  // the number of finemesh (x)
  int nny;  // the number of finemesh (y)
  int nnz;  // the number of finemesh (z)
  int nnxnny;

  double vi;
  double vs;
  static const double vp_ave;

  double tausT;
  double tausL;
  double taumninv;
  double tauiT;  // D/dx/dx*dt
  double tauiL;
  double taunsrT;
  double taunsrL;

  void computeIci(void);
  void computeIcnsr(void);
  double *Ici, *Icnsr;
#ifdef ___NO_CS_BUFFER
  double *csmn;
  void computecsmn(void);
#else
  void computeIcs(void);
  double *Ics;
  double *Idps;
#endif
  int bino(double num, double p, int ii);
  double calcvp(double mean, double std, double lim1, double lim2, int ii);

  double dt;
  double vup, kup, KNSR;
  double vnaca;
  double Jmax;
  double taup;
  double tausi;
  double tautr;

  double gca;  // scale 1

  double cao;

  double gleak;
  double gcabk;
  double qslcap;
  double BCSQN;

  //  RyR parameters
  double Kc;
  double nM;
  double nD;
  double hh;
  double KK;
  double rhoinf;

  double Ku;
  double Kb;
  double tauu;
  double taub;
  double tauc1;
  double tauc2;
  double BCSQN0;

  int seed;
  bool initialized;
  int bc;

  // #ifdef ___SIGMOID
  double Kcp;
  double pedk12;
  double pedk43;
  // #endif
  int NoCaL;  // # of CaL channel per CaRU

  double MaxSR;
  double MinSR;
  double ec50SR;
  double hkosrca;

  void delarray(void);

#ifdef ___EGTA
  double BEGTA;
#endif

 public:
  unsigned int *xsx, *xsy, *xsz, *xsw;
  int nx;        // the number of CRU (x)
  int ny;        // the number of CRU (y)
  int nz;        // the number of CRU (z)
  int finemesh;  // X times fine mesh (must be odd number , must be >3)
  double xi;

  int n;   // the number of CRU x*y*z
  int nn;  // the number of fine mesh

  static const double vjsr;
  double vnsr;

  CSubcell(int sizex = 65, int sizey = 27, int sizez = 11, int fmesh = 1,
           double xii = 1);
  void init(double initci = 0.1, double initcj = 750);
  virtual ~CSubcell();

  //  double *ci,*cs,*cp,*cjsr,*cnsr,*cati,*cats;
  double *ci, *cs, *cp, *cjsr, *cnsr, *cati;
  double *cscp1, *cscp2, *Itr;
#ifdef ___DETERMINISTIC
  double *c1, *c2, *i1ca, *i1ba, *i2ca, *i2ba, *fryr1, *fryr2, *fryr3;
#else
  int *y, *ryr1, *ryr2, *ryr3, *nryr;
#endif
  double *vp;
  double *Jmaxx;

  int *crupos;

  void pace(double v, double nai);

  double iupave, icaave, incxave, irave, ileakave, icabkave, islcapave;
  double computeaveci(void);
  double computeavecs(void);
  double computeavecnsr(void);
  CSubcell &operator=(const CSubcell &sc);

  void setdt(double newdt) { dt = newdt; };
  double getdt(void) { return dt; };
  void setgca(double newgca) { gca = newgca; };
  double getgca(void) { return gca; };
  void setgncx(double newgncx) { vnaca = newgncx; };
  double getgncx(void) { return vnaca; };
  void setJmax(double newJmax);
  double getJmax(void) { return Jmax; };
  void setvup(double newvup) { vup = newvup; };
  double getvup(void) { return vup; };
  void setkup(double newkup) { kup = newkup; };
  double getkup(void) { return kup; };
  void setKNSR(double newKNSR) { KNSR = newKNSR; };
  double getKNSR(void) { return KNSR; };

  void setgleak(double newgleak) { gleak = newgleak; };
  double getgleak(void) { return gleak; };
  void setcao(double newcao) { cao = newcao; };
  double getcao(void) { return cao; };
  void setBCSQN(double newBCSQN) { BCSQN = newBCSQN; };
  double getBCSQN(void) { return BCSQN; };

  void setgcabk(double newgcabk) { gcabk = newgcabk; };
  double getgcabk(void) { return gcabk; };
  void setqslcap(double newqslcap) { qslcap = newqslcap; };
  double getqslcap(void) { return qslcap; };

  void setKc(double newKc) { Kc = newKc; };
  double getKc(void) { return Kc; };
  void setnM(double newnM) { nM = newnM; };
  double getnM(void) { return nM; };
  void setnD(double newnD) { nD = newnD; };
  double getnD(void) { return nD; };
  void sethh(double newhh) { hh = newhh; };
  double gethh(void) { return hh; };
  void setKK(double newKK) { KK = newKK; };
  double getKK(void) { return KK; };
  void setrhoinf(double newrhoinf) { rhoinf = newrhoinf; };
  double getrhoinf(void) { return rhoinf; };
  void setKu(double newKu) { Ku = newKu; };
  double getKu(void) { return Ku; };
  void setKb(double newKb) { Kb = newKb; };
  double getKb(void) { return Kb; };
  void settauu(double newtauu) { tauu = newtauu; };
  double gettauu(void) { return tauu; };
  void settaub(double newtaub) { taub = newtaub; };
  double gettaub(void) { return taub; };
  void settauc1(double newtauc1) { tauc1 = newtauc1; };
  double gettauc1(void) { return tauc1; };
  void settauc2(double newtauc2) { tauc2 = newtauc2; };
  double gettauc2(void) { return tauc2; };
  void setBCSQN0(double newBCSQN0) { BCSQN0 = newBCSQN0; };
  double getBCSQN0(void) { return BCSQN0; };

  void setvi(double newvi) { vi = newvi; }
  double getvi(void) { return vi; }
  void setvs(double newvs) {
    vs = newvs;
    for (int id = 0; id < n; id++) {
      cscp2[crupos[id]] = 1 / taup * vp[id] / vs;
    }
  }
  double getvs(void) { return vs; }

  double gettausi(void) { return tausi; }
  double settausi(double newval) {
    tausi = newval;
    return tausi;
  }

  double gettautr(void) { return tautr; }
  double settautr(double newval) {
    tautr = newval;
    return tautr;
  }

  // #ifdef ___SIGMOID
  void setKcp(double newKcp) { Kcp = newKcp; };
  double getKcp(void) { return Kcp; };
  void setpedk12(double newpedk12) { pedk12 = newpedk12; };
  double getpedk12(void) { return pedk12; };
  void setpedk43(double newpedk43) { pedk43 = newpedk43; };
  double getpedk43(void) { return pedk43; };
  // #endif

  double getMaxSR(void) { return MaxSR; };
  void setMaxSR(double newMaxSR) { MaxSR = newMaxSR; };
  double getMinSR(void) { return MinSR; };
  void setMinSR(double newMinSR) { MinSR = newMinSR; };
  double getec50SR(void) { return ec50SR; };
  void setec50SR(double newec50SR) { ec50SR = newec50SR; };
  double gethkosrca(void) { return hkosrca; };
  void sethkosrca(double newhkosrca) { hkosrca = newhkosrca; };

  // functions: diffusion
  double gettausT(void) { return tausT; }
  double settausT(double newval) {
    tausT = newval;
    taumninv = 2 / tausL + 4 / tausT;
    return tausT;
  }
  double gettausL(void) { return tausL; }
  double settausL(double newval) {
    tausL = newval;
    taumninv = 2 / tausL + 4 / tausT;
    return tausL;
  }
  double gettauiT(void) { return tauiT; }
  double settauiT(double newval) {
    tauiT = newval;
    return tauiT;
  }
  double gettauiL(void) { return tauiL; }
  double settauiL(double newval) {
    tauiL = newval;
    return tauiL;
  }
  double gettaunsrT(void) { return taunsrT; }
  double settaunsrT(double newval) {
    taunsrT = newval;
    return taunsrT;
  }
  double gettaunsrL(void) { return taunsrL; }
  double settaunsrL(double newval) {
    taunsrL = newval;
    return taunsrL;
  }

  double calcjcalave(void) { return icaave; };
  double calcjnacaave(void) { return incxave; };
  double calcicabkave(void) { return icabkave; };
  double calcislcapave(void) { return islcapave; };
  void srand(int sed = 0);
  void setboundary(int bcc = 3);
  void resetBuffer(void);

  int getNoCaL(void) { return NoCaL; }
  int setNoCaL(int newval) {
    NoCaL = newval;
    return NoCaL;
  }

  double couplingstrength1;
  double couplingstrength2;
  double csx;
  double pox;

#ifdef ___EGTA
  double *caEGTAi, *caEGTAs;
  double setBEGTA(double newBEGTA) {
    BEGTA = newBEGTA;
    resetBuffer();
    return BEGTA;
  }
  double getBEGTA(void) { return BEGTA; }
#endif

#ifdef ___NCX
  double NCXalpha;
#endif

#ifdef ___DEBUG
  bool bSTOP;
#endif
};
#endif
