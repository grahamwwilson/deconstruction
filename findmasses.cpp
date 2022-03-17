#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TRandom1.h"
using namespace std;

// electroweakino.f90
extern"C" {
void electroweakino_(double *M1, double *M2, double *mu, double *tanb);
}

// struct for mval common block
extern "C" {
   extern struct{
      double n1,n2,n3,n4;
      double c1,c2;
      double v1[4],v2[4],v3[4],v4[4];
      double U[4], V[4];
      double sn1,sn2,sn3,sn4;
   } mval_;                                // Filled using mval common block
}

void dumpNMatrix(){
  cout << "N1j " << mval_.v1[0] << " " << mval_.v1[1] << " " 
                << mval_.v1[2] << " " << mval_.v1[3] << endl;
  cout << "N2j " << mval_.v2[0] << " " << mval_.v2[1] << " " 
                << mval_.v2[2] << " " << mval_.v2[3] << endl;
  cout << "N3j " << mval_.v3[0] << " " << mval_.v3[1] << " " 
                << mval_.v3[2] << " " << mval_.v3[3] << endl;
  cout << "N4j " << mval_.v4[0] << " " << mval_.v4[1] << " " 
                << mval_.v4[2] << " " << mval_.v4[3] << endl;
}
void dumpUMatrix(){
  cout << "U " << mval_.U[0] << " " << mval_.U[1] << endl;
  cout << "  " << mval_.U[2] << " " << mval_.U[3] << endl;
}
void dumpVMatrix(){
  cout << "V " << mval_.V[0] << " " << mval_.V[1] << endl;
  cout << "  " << mval_.V[2] << " " << mval_.V[3] << endl;
}

int main()
{

 const float xdal = 1.0/sqrt(3.0);
 TFile* f;
 f = (TFile*) new TFile("histos_ewinos_Point.root","recreate");

 TH1D* hn1type = (TH1D*) new TH1D("hn1type","; N1 type ", 3, 0.5, 3.5);
 TH1D* hn2type = (TH1D*) new TH1D("hn2type","; N2 type ", 3, 0.5, 3.5);
 TH1D* hn3type = (TH1D*) new TH1D("hn3type","; N3 type ", 3, 0.5, 3.5);
 TH1D* hn4type = (TH1D*) new TH1D("hn4type","; N4 type ", 3, 0.5, 3.5);

 TH1D* hn1 = (TH1D*) new TH1D("hn1","; Mass(N1) (GeV) ", 100, 0.0, 1000.0);
 TH1D* hn2 = (TH1D*) new TH1D("hn2","; Mass(N2) (GeV) ", 100, 0.0, 1000.0);
 TH1D* hn3 = (TH1D*) new TH1D("hn3","; Mass(N3) (GeV) ", 100, 0.0, 1000.0);
 TH1D* hn4 = (TH1D*) new TH1D("hn4","; Mass(N4) (GeV) ", 100, 0.0, 1000.0);
 TH1D* hc1 = (TH1D*) new TH1D("hc1","; Mass(C1) (GeV) ", 100, 0.0, 1000.0);
 TH1D* hc2 = (TH1D*) new TH1D("hc2","; Mass(C2) (GeV) ", 100, 0.0, 1000.0);

 TH1D* hdm = (TH1D*) new TH1D("hdm","; DeltaM (N2-N1) (GeV)",100, -10.0, 190.0);
 TH1D* hdc = (TH1D*) new TH1D("hdc","; DeltaM (C1-N1) (GeV)",100, -10.0, 190.0);
 TH1D* hdn = (TH1D*) new TH1D("hdn","; DeltaM (N2-C1) (GeV)",100, -10.0, 190.0);

 TH2D* hdcn = (TH2D*) new TH2D("hdcn","; Mass (C1) (GeV); Mass (N1) (GeV)", 
                               100, 100.0, 750.0, 100, 0.0, 600.0);

 TH2D* hdal1 = (TH2D*) new TH2D("hdal1","LSP Mixing (B-W-H)", 
                               120, -0.6, 0.6, 105, -0.35, 0.70);
 TH2D* hdal2 = (TH2D*) new TH2D("hdal2","N2 Mixing (B-W-H)", 
                               120, -0.6, 0.6, 105, -0.35, 0.70);
 TH2D* hdal3 = (TH2D*) new TH2D("hdal3","N3 Mixing (B-W-H)", 
                               120, -0.6, 0.6, 105, -0.35, 0.70);
 TH2D* hdal4 = (TH2D*) new TH2D("hdal4","N4 Mixing (B-W-H)", 
                               120, -0.6, 0.6, 105, -0.35, 0.70);

 TRandom1* rg = (TRandom1*) new TRandom1(4358,3);

 unsigned long int ntrials = 0;
 int nthrown = 0;
 
 while (nthrown < 100){

      ntrials++;

 // Maybe should sample uniformly in logM 
 // Should have a function that returns suitable (M1, M2, mu, tanb) points
      double M1 = double(rg->Uniform(50.0,400.0));
      double M2 = double(rg->Uniform(50.0,400.0));
      double mu = double(rg->Uniform(50.0,400.0));
      double test1 = rg->Uniform(1.0);
      double test2 = rg->Uniform(1.0);
      double test3 = rg->Uniform(1.0);
      double test4 = rg->Uniform(1.0);
      if(test1<0.5)M1=-M1;                       // Sample +/- M1 uniformly
      if(test2<0.5)M2=-M2;                       // Sample +/- M2 uniformly
      if(test3<0.5)mu=-mu;                       // Sample +/- mu uniformly
      double tanb = double(rg->Uniform(1.0,10.0));
      if(test4<0.5)tanb = 1.0/tanb;
 //     if(abs(mu) > M1 || abs(mu) > M2) continue;
 
 // Calculate the corresponding mass spectrum and couplings
      electroweakino_(&M1, &M2, &mu, &tanb);     // Fortran subroutine
/*      cout << "N1 " << mval_.n1 << " N2 " << mval_.n2 
           << " N3 " << mval_.n3 << " N4 " << mval_.n4 << endl;
      cout << "C1 " << mval_.c1 << " C2 " << mval_.c2 << endl;      */

//      if( abs(mval_.n1 - 100.0) > 5.0 || abs(mval_.c1 - 105.0) > 5.0 || abs(mval_.c2 - 200.0) > 5.0) continue;
      if( abs(mval_.n1 - 220.0) > 5.0 || abs(mval_.c1 - 235.0) > 5.0 || abs(mval_.n2 - 250.0) > 5.0) continue;

      float fb1 = pow(mval_.v1[0],2);
      float fw1 = pow(mval_.v1[1],2);
      float fh1 = pow(mval_.v1[2],2) + pow(mval_.v1[3],2);
      
      cout << "SELECTED! " << nthrown << endl; 
      cout << "M1 " << M1 << " M2 " << M2 
           << " mu " << mu << " tanb " << tanb << endl;
      cout << "N1 " << mval_.n1 << " N2 " << mval_.n2 
           << " N3 " << mval_.n3 << " N4 " << mval_.n4 << endl;
      cout << "C1 " << mval_.c1 << " C2 " << mval_.c2 << endl;
      cout << "LSP " << fb1 << " " << fw1 << " "  << fh1 << " (B, W, H) " << endl;
      cout << "signs " << mval_.sn1/mval_.n1 << " " 
                       << mval_.sn2/mval_.n2 << " " 
                       << mval_.sn3/mval_.n3 << " " 
                       << mval_.sn4/mval_.n4 << " " << endl;

      dumpNMatrix();
      dumpUMatrix();
      dumpVMatrix();

      float binoLSP = fb1 - max(fw1, fh1);
      float winoLSP = fw1 - max(fb1, fh1);
      float hinoLSP = fh1 - max(fb1, fw1);
      float checksum = fb1 + fw1 + fh1;
      cout << "Checksum = " << checksum << endl;
 // Dalitz plot coordinates. Draw equilateral triangle with height 1 with 
 // origin corresponding to fb=fw=fh=1/3 with the 
 // coordinate satisfying fb+fw+fh = 1.
      float xval1 = (fb1-fw1)/sqrt(3.0);
      float yval1 = fh1 - (1.0/3.0);
      hdal1->Fill(xval1,yval1);
      cout << binoLSP << " " << winoLSP << " " << hinoLSP << endl;
      cout << " " << endl;

      int type1 = 0;
      if (fb1 > max(fw1, fh1) ) type1 = 1;
      if (fw1 > max(fb1, fh1) ) type1 = 2;
      if (fh1 > max(fb1, fw1) ) type1 = 3;
      hn1type->Fill(type1);

// N2 mixing
      float fb2 = pow(mval_.v2[0],2);
      float fw2 = pow(mval_.v2[1],2);
      float fh2 = pow(mval_.v2[2],2) + pow(mval_.v2[3],2);
      float xval2 = (fb2-fw2)/sqrt(3.0);
      float yval2 = fh2 - (1.0/3.0);
      hdal2->Fill(xval2,yval2);
      int type2 = 0;
      if (fb2 > max(fw2, fh2) ) type2 = 1;
      if (fw2 > max(fb2, fh2) ) type2 = 2;
      if (fh2 > max(fb2, fw2) ) type2 = 3;
      hn2type->Fill(type2);

// N3 mixing
      float fb3 = pow(mval_.v3[0],2);
      float fw3 = pow(mval_.v3[1],2);
      float fh3 = pow(mval_.v3[2],2) + pow(mval_.v3[3],2);
      float xval3 = (fb3-fw3)/sqrt(3.0);
      float yval3 = fh3 - (1.0/3.0);
      hdal3->Fill(xval3,yval3);
      int type3 = 0;
      if (fb3 > max(fw3, fh3) ) type3 = 1;
      if (fw3 > max(fb3, fh3) ) type3 = 2;
      if (fh3 > max(fb3, fw3) ) type3 = 3;
      hn3type->Fill(type3);

// N4 mixing
      float fb4 = pow(mval_.v4[0],2);
      float fw4 = pow(mval_.v4[1],2);
      float fh4 = pow(mval_.v4[2],2) + pow(mval_.v4[3],2);
      float xval4 = (fb4-fw4)/sqrt(3.0);
      float yval4 = fh4 - (1.0/3.0);
      hdal4->Fill(xval4,yval4);
      int type4 = 0;
      if (fb4 > max(fw4, fh4) ) type4 = 1;
      if (fw4 > max(fb4, fh4) ) type4 = 2;
      if (fh4 > max(fb4, fw4) ) type4 = 3;
      hn4type->Fill(type4);

      hn1->Fill(mval_.n1);
      hn2->Fill(mval_.n2);
      hn3->Fill(mval_.n3);
      hn4->Fill(mval_.n4);
      hc1->Fill(mval_.c1);
      hc2->Fill(mval_.c2);

      hdm->Fill(mval_.n2-mval_.n1);
      hdc->Fill(mval_.c1-mval_.n1);
      hdn->Fill(mval_.n2-mval_.c1);     
      hdcn->Fill(mval_.c1,mval_.n1);
      nthrown++;

 }
 f->Write();
 cout << "Ntrials " << ntrials << endl;

 return 0;
}
