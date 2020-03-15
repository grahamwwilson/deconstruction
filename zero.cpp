#include "Math/Functor.h"
#include "Math/RootFinder.h"
#include <TMath.h>
#include <iostream>
#include <iomanip>

using namespace std;

double mZ=91.1876; double s2W=0.232;
double M1=254.0; double M2=480.0; double mu=768.0; double tanb=5.0;

#include "MyFunction.C"
 
int main() {
 
   // RootFinder with Base Functions
   ROOT::Math::Functor1D f1D(&MyFunction);
 
   // Create the RootFinder and find the root
   ROOT::Math::RootFinder rf(ROOT::Math::RootFinder::kGSL_BRENT);
   rf.SetFunction(f1D, 200.0, 300.0); 
   rf.Solve(1000,0.0,1.0e-14);
   int status = rf.Status();
   int niter = rf.Iterations();

   cout << "root finding status : " << status << " niterations = " << niter << endl; 

// Use FORTRAN compatibility output using examples from 
// page 29 of Barton and Nackmann, Scientific and Engineering C++
   cout << setiosflags(ios::showpoint | ios::uppercase);
   cout << setiosflags(ios::fixed);
   cout.precision(16);

   cout << rf.Root() << endl;

   cout << rf.Name() << endl;
//   cout << rf.XLower() << endl;
//   cout << rf.XUpper() << endl;

   // Find the solution vector from the root position
   double m = rf.Root();
   double result;
   result = EvaluateFunction(m);
   cout << " f(m) = " << result  << endl;

   return 0;

}
