#ifdef _OPENMP
#include <omp.h>
#endif

#include <fstream>
#include <iostream>
#include <sstream>
using namespace std;
#include "ap.h"
#include "log.h"
#include "recsubcell.h"
#include "subcell.h"

// for open and close functions for /dev/urandom
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

int main(int argc, char *argv[]) {
  const int NoCaL = 10;

  const int fmesh = 1;
  const int nx = 65;
  const int ny = 27;
  const int nz = 11;
  //  const int nx=5;
  //  const int ny=5;
  //  const int nz=5;
  const int nn = nx * ny * nz;

  Logging *llog = new Logging(LOG_UNIQUE);

  int randomData = open("/dev/random", O_RDONLY);
  int r;
  ssize_t ret = read(randomData, &r, sizeof r);
  close(randomData);
  r = abs(r);

  double cps1 = 502.127;
  double cps2 = 5.02111;

  char tmpchar[255];
  sprintf(tmpchar, "seed=%d, cps1=%f, cps2=%f", r, cps1, cps2);
  llog->Note(tmpchar);

  CSubcell sc(nx, ny, nz, fmesh, 1);
  sc.srand(r);

  sc.setdt(0.01);
  double dt = sc.getdt();

  sc.setNoCaL(NoCaL);
  sc.init(0.1, 700);

  sc.couplingstrength1 = cps1;
  sc.couplingstrength2 = cps2;
  sc.pox = 0.1;
  sc.csx = 10;

  CActionPotential ap;
  ap.v = -80;
  ap.nai = 8;

  RECSUBCELLPARAM param;

  CRecSubcell rec(&sc, &ap, &param);

  int Tn = 1000 / dt;
  for (int tn = 0; tn < Tn; tn++) {
    sc.pace(ap.v, ap.nai);
  }

  for (int i = 0; i < nn; i++) {
    sc.ci[i] = sc.cs[i] = sc.cp[i] = 0.1;
    sc.cjsr[i] = sc.cnsr[i] = 700;
  }

  CSubcell sc2(nx, ny, nz, fmesh, 1);
  sc2 = sc;
  double activation[200];
  double peakicalvalue[200];
  double releasevalue[200];
  for (int i = 0; i < 200; i++) {
    activation[i] = peakicalvalue[i] = releasevalue[i] = 0;
  }

  char fname[255];
  for (int v = -80; v <= 80; v += 5) {
    sprintf(fname, "ical%d.txt", v);
    ofstream osical(fname);
    sc = sc2;
    double peak = 0;
    double peakical = 0;
    double peakical2 = 0;
    double ir = 0;
    ap.v = v;
    Tn = 300 / dt;
    for (int tn = 0; tn < Tn; tn++) {
      if (tn % 100 == 0)
        osical << tn * dt << "\t" << sc.calcjcalave() << "\t"
               << sc.computeaveci() << endl;
      sc.pace(ap.v, ap.nai);
      int c0 = 0;
      for (int id = 0; id < nn; id++) {
        for (int j = 0; j < NoCaL; j++) {
          switch (sc.y[id * NoCaL + j]) {
            case 0:
              c0++;
              break;
          }
        }
      }
      double o = c0 / double(nn * NoCaL);
      if (o > peak) {
        peak = o;
      }

      if (fabs(sc.calcjcalave()) > peakical) {
        peakical = fabs(sc.calcjcalave());
        peakical2 = sc.calcjcalave();
      }
      ir += sc.irave;
    }
    activation[v + 80] = peak;
    peakicalvalue[v + 80] = peakical2;
    releasevalue[v + 80] = ir * dt;
  }

  sprintf(fname, "%d_activation.txt", r);
  ofstream osact(fname);
  sprintf(fname, "%d_iv.txt", r);
  ofstream osiv(fname);

  for (int v = -80; v <= 80; v += 5) {
    osact << v << "\t" << activation[v + 80] << endl;
    osiv << v << "\t" << peakicalvalue[v + 80] << "\t" << releasevalue[v + 80]
         << endl;
  }

  sprintf(fname, "%d_count.txt", r);
  ofstream oscnt(fname);
  for (int v = -70; v <= 40; v += 10) {
    sprintf(fname, "%d_opentime%d.txt", r, v);
    ofstream osopentime(fname);
    sc = sc2;
    int cnt[NoCaL + 1];
    for (int i = 0; i <= NoCaL; i++) cnt[i] = 0;
    Tn = 10 / dt;

    bool openstate[nn];
    int maxcount[nn];
    int starttime[nn];
    for (int id = 0; id < nn; id++) {
      openstate[id] = false;
      maxcount[id] = 0;
    }

    for (int tn = 0; tn < Tn; tn++) {
      //      if (tn%1000==0)
      //        cout<<tn*dt<<endl;
      ap.v = v;
      sc.pace(ap.v, ap.nai);
      //      if (tn%tnrec==0)
      for (int id = 0; id < nn; id++) {
        int counter = 0;
        for (int j = 0; j < NoCaL; j++) {
          if (sc.y[id * NoCaL + j] == 0) counter++;
        }
        if (openstate[id] == false && counter > 0) {
          openstate[id] = true;
          maxcount[id] = counter;
          starttime[id] = tn;
        } else if (openstate[id] == true && counter == 0) {
          osopentime << (tn - starttime[id]) * dt << "\t" << maxcount[id]
                     << endl;
          openstate[id] = false;
          cnt[maxcount[id]]++;
          maxcount[id] = 0;
        } else if (openstate[id] == true && counter > 0) {
          if (maxcount[id] < counter) maxcount[id] = counter;
        }
      }
    }
    oscnt << v;
    for (int i = 1; i <= NoCaL; i++) oscnt << "\t" << cnt[i];
    oscnt << endl;
  }

  delete llog;
  return 0;
}
