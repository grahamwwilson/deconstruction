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
      double M1 = 254.0;
      double M2 = 480.0;
      double mu = 768.0;
      double tanb = 5.0;
 
 // Calculate the corresponding mass spectrum and couplings
      electroweakino_(&M1, &M2, &mu, &tanb);    // Fortran subroutine

      float fb1 = pow(mval_.v1[0],2);
      float fw1 = pow(mval_.v1[1],2);
      float fh1 = pow(mval_.v1[2],2) + pow(mval_.v1[3],2);
       
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
//      hdal1->Fill(xval1,yval1);
      cout << binoLSP << " " << winoLSP << " " << hinoLSP << endl;
      cout << " " << endl;

      int type1 = 0;
      if (fb1 > max(fw1, fh1) ) type1 = 1;
      if (fw1 > max(fb1, fh1) ) type1 = 2;
      if (fh1 > max(fb1, fw1) ) type1 = 3;
//      hn1type->Fill(type1);

// N2 mixing
      float fb2 = pow(mval_.v2[0],2);
      float fw2 = pow(mval_.v2[1],2);
      float fh2 = pow(mval_.v2[2],2) + pow(mval_.v2[3],2);
      float xval2 = (fb2-fw2)/sqrt(3.0);
      float yval2 = fh2 - (1.0/3.0);
//      hdal2->Fill(xval2,yval2);
      int type2 = 0;
      if (fb2 > max(fw2, fh2) ) type2 = 1;
      if (fw2 > max(fb2, fh2) ) type2 = 2;
      if (fh2 > max(fb2, fw2) ) type2 = 3;
//      hn2type->Fill(type2);

// N3 mixing
      float fb3 = pow(mval_.v3[0],2);
      float fw3 = pow(mval_.v3[1],2);
      float fh3 = pow(mval_.v3[2],2) + pow(mval_.v3[3],2);
      float xval3 = (fb3-fw3)/sqrt(3.0);
      float yval3 = fh3 - (1.0/3.0);
//      hdal3->Fill(xval3,yval3);
      int type3 = 0;
      if (fb3 > max(fw3, fh3) ) type3 = 1;
      if (fw3 > max(fb3, fh3) ) type3 = 2;
      if (fh3 > max(fb3, fw3) ) type3 = 3;
//      hn3type->Fill(type3);

// N4 mixing
      float fb4 = pow(mval_.v4[0],2);
      float fw4 = pow(mval_.v4[1],2);
      float fh4 = pow(mval_.v4[2],2) + pow(mval_.v4[3],2);
      float xval4 = (fb4-fw4)/sqrt(3.0);
      float yval4 = fh4 - (1.0/3.0);
//      hdal4->Fill(xval4,yval4);
      int type4 = 0;
      if (fb4 > max(fw4, fh4) ) type4 = 1;
      if (fw4 > max(fb4, fh4) ) type4 = 2;
      if (fh4 > max(fb4, fw4) ) type4 = 3;
//      hn4type->Fill(type4);

/*
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
*/
//      nthrown++;

// f->Write();
// cout << "Ntrials " << ntrials << endl;

 return 0;
}
