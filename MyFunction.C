double MyFunction(double m){

// Set up quartic equation

   double A,B,C,D;
   double value;
   double beta,sin2b;

   beta = atan2(tanb,1.0);
   sin2b = sin(2.0*beta);

   A = M1 + M2;
   B = M1*M2 - mu*mu - mZ*mZ;
   C = -A*mu*mu + mu*mZ*mZ*sin2b - mZ*mZ*(M1+(M2-M1)*s2W);
   D = -M1*M2*mu*mu + mu*mZ*mZ*sin2b*(M1+(M2-M1)*s2W);
  
   value = m*m*m*m - A*m*m*m + B*m*m - C*m + D; 

   return value;
}

double EvaluateFunction(double m){

// Set up quartic equation

   double A,B,C,D;
   double value;
   double beta,sin2b;
   double Ap,Bp,Cp,Dp;
   double value2;
   double mp;
   double C1,C2,C3;

   beta = atan2(tanb,1.0);
   sin2b = sin(2.0*beta);

   A = M1 + M2;
   B = M1*M2 - mu*mu - mZ*mZ;
   C = -A*mu*mu + mu*mZ*mZ*sin2b - mZ*mZ*(M1+(M2-M1)*s2W);
   D = -M1*M2*mu*mu + mu*mZ*mZ*sin2b*(M1+(M2-M1)*s2W);

   Ap = A/100.0;
   Bp = B/(100.0*100.0);
   Cp = C/(100.0*100.0*100.0);
   Dp = D/(100.0*100.0*100.0*100.0);
  
   cout << "A,B,C,D:    " << A << " " << B << " " << C << " " << D << endl;
   cout << "A,B,C,D('): " << Ap << " " << Bp << " " << Cp << " " << Dp << endl;

   value = m*m*m*m - A*m*m*m + B*m*m - C*m + D; 

   mp = m/100.0;
   value2 = mp*mp*mp*mp - Ap*mp*mp*mp + Bp*mp*mp - Cp*mp + Dp;
   cout << "value2 " << value2 << endl;

// Diagnostics 
   C1 = -A*mu*mu;
   C2 =  mu*mZ*mZ*sin2b;
   C3 = - mZ*mZ*(M1+(M2-M1)*s2W);
   cout << "C terms: " << C1 << " " << C2 << " " << C3 << endl;

   return value;
}
