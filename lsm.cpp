#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "utils.hpp"
#include "lsm.hpp"
//#include "omp.h"

std::vector<double>
LSM::calcEfuncs(const double q) {
    // Computes the "extra" functions of Q which come in beyond the
    // Zel'dovich approximation.
    const int Nk=kLin.size();
    const double kmax=exp(kLin[Nk-1]);
    int Nint=(int)(8*kmax*q+512);
    if (Nint>=20000) Nint=20000;
    const double hh=(kLin[Nk-1]-kLin[0])/Nint;
    double sum0=0,sum1=0,sum2=0,sum3=0,sum4=0,sum5=0,sum6=0,sum7=0,
           sum8=0,sum9=0,sum10=0,sum11=0,sum12=0,sumS=0;
    //std::ofstream files("letR1",std::ios::trunc);
#pragma omp parallel for reduction(+:sum0,sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,sum9,sum10,sum11,sum12,sumS)
    for (int i=1; i<Nint; ++i) {
      double xx = kLin[0]+i*hh;
      double kk = exp(xx);
      double k2 = kk*kk;
      double kq = kk*q;
      double R1 = R1spl(kk);
      double R2 = R2spl(kk);
      double Q1 = Q1spl(kk);
      double Q2 = Q2spl(kk);
      double Q5 = Q5spl(kk);
      double Q8 = Q8spl(kk);
      double Qs = Qsspl(kk);
      double QI = QIspl(kk);
      double R12 = R12spl(kk);
      double RI = RIspl(kk);
      //extra Q splines to store
     // double Q3 = Q3spl(kk);
     // double Q4 = Q4spl(kk);
     // double Q6 = Q6spl(kk);
     // double Q7 = Q7spl(kk);
     // double Q9 = Q9spl(kk);
     // double Q10 = Q10spl(kk);
     // double Q11 = Q11spl(kk);
     // double Q12 = Q12spl(kk);
     // double Q13 = Q13spl(kk);
      //std::cout << exp(xx) << std::endl;
      //files << kLin[0] << std::endl;
      //files.close();      
      double ap = cos(M_PI/2.*kk/kmax);
      std::vector<double> jl=sphBess(kk*q);
      double j1 = kq*jl[1];
      double j2,j3;
      if (kq<0.1) {
        j2 = pow(kq,2.0)/15.  - pow(kq,4.0)/210.;
        j3 = pow(kq,3.0)/105. - pow(kq,5.0)/1890.;
      }
      else {
        j2 = 3.*jl[1]-jl[0];
        j3 = 5.*j2/(kq)-(kq)*jl[1];
      }
      int wt= 2+2*(i%2);
      //Changes by George to account for MG correlators
      sum0 += kk*wt*(9./98.*Q1*(2./3.-2*jl[1]));		// X^{(22)}
      sum1 += kk*wt*(5./21.*R1*(2./3.-2*jl[1]));		// X^{(13)}
      sum2 += kk*wt*(9./98.*Q1*(-2*jl[0]+6*jl[1]));		// Y^{(22)}
      sum3 += kk*wt*(5./21.*R1*(-2*jl[0]+6*jl[1]));		// Y^{(13)}
      //sum4 += kk*wt*(2*(R1-R2)+3*R1*jl[0]-3*(3*R1+4*R2+2*Q5)*jl[1])/14.;//X1210
      //sum5 += kk*wt*(3*R1+4*R2+2*Q5)*(jl[0]-3*jl[1])*(-3./14.);// Y_{10}^{(12)}
      sum4 += kk*wt*(2*(RI-R2)+3*RI*jl[0]-3*(RI+2*R2+2*R12+2*Q5)*jl[1])/14.;//X1210
      sum5 += kk*wt*(RI+2*R2+2*R12+2*Q5)*(jl[0]-3*jl[1])*(-3./14.);// Y_{10}^{(12)}
      //sum6 +=    wt*(R1*j1)*(-3./7.);			// V_1^{(112)}
      sum6 +=    wt*(RI*j1)*(-3./7.);			// V_1^{(112)}
      //sum7 +=    wt*(Q1*j1)*(-3./7.);			// V_3^{(112)}
      sum7 +=    wt*(QI*j1)*(-3./7.);			// V_3^{(112)}
      //sumS +=    wt*(2*R1+4*R2+Q1+2*Q2)*(3./7.*j2/(kk*q));	// S^{(112)}
      sumS +=    wt*(2*RI+4*R2+QI+2*Q2)*(3./7.*j2/(kk*q));	// S^{(112)}
      //sum8 +=    wt*(2*R1+4*R2+Q1+2*Q2)*j3*(-3./7.);	// T^{(112)}
      sum8 +=    wt*(2*RI+4*R2+QI+2*Q2)*j3*(-3./7.);	// T^{(112)}
      sum9 += k2*wt*(R1*j1)*(-5./21.);			// U^{(3)}
      sum10+= k2*wt*(Q8*j1)*(-3./7.);			// U_{20}^{(2)}
      //sum11+= k2*wt*((R1+R2)*j1)*(-6./7.);		// U_{11}^{(2)}
      sum11+= k2*wt*((R12)*j1)*(-6./7.);		// U_{11}^{(2)}
      sum12+= k2*wt*(Qs*j1)*(-2./7.)*ap;		// Shear term
    }
    //std::cout  << sum0 << std::endl;
    //files.close();
    sum6 += sumS;
    sum7 += sumS;
    std::vector<double> sum(16);
    sum[ 1] = sum0 * hh/3.0/(2*M_PI*M_PI);
    sum[ 2] = sum1 * hh/3.0/(2*M_PI*M_PI);
    sum[ 4] = sum2 * hh/3.0/(2*M_PI*M_PI);
    sum[ 5] = sum3 * hh/3.0/(2*M_PI*M_PI);
    sum[ 6] = sum4 * hh/3.0/(2*M_PI*M_PI);
    sum[ 7] = sum5 * hh/3.0/(2*M_PI*M_PI);
    sum[ 8] = sum6 * hh/3.0/(2*M_PI*M_PI);
    sum[ 9] = sum7 * hh/3.0/(2*M_PI*M_PI);	
    sum[10] = sum8 * hh/3.0/(2*M_PI*M_PI);	
    sum[12] = sum9 * hh/3.0/(2*M_PI*M_PI);	
    sum[13] = sum10* hh/3.0/(2*M_PI*M_PI);	
    sum[14] = sum11* hh/3.0/(2*M_PI*M_PI);	
    sum[15] = sum12* hh/3.0/(2*M_PI*M_PI);	
    // Now tabulate the pieces going as Plin.
    sum0=sum1=sum2=0;
#pragma omp parallel for reduction(+:sum0,sum1,sum2)
    for (int i=1; i<Nint; ++i) {
      double xx = kLin[0]+i*hh;
      double kk = exp(xx);
      double k2 = kk*kk;
      double kq = kk*q;
      int    jj = (int)(i*hh*dkinv);
      if (jj>=pLin.size()-2) jj=pLin.size()-2;
      double pk = exp(pLin[jj]+(xx-kLin[jj])*
                     (pLin[jj+1]-pLin[jj])/(kLin[jj+1]-kLin[jj]));
      std::vector<double> jl=sphBess(kk*q);
      double j1 = kq*jl[1];
      int wt= 2+2*(i%2);
      sum0 += kk*wt*pk*(2./3.-2*jl[1]);		// X^{(11)}
      sum1 += kk*wt*pk*(-2.*jl[0]+6*jl[1]);	// Y^{(11)}
      sum2 += k2*wt*pk*(-j1);			// U^{(1)}
    }
    sum[ 0] = sum0 * hh/3.0/(2*M_PI*M_PI);
    sum[ 3] = sum1 * hh/3.0/(2*M_PI*M_PI);
    sum[11] = sum2 * hh/3.0/(2*M_PI*M_PI);
    return(sum);
}

void
LSM::tabulateEfuncs() {
    // Tabulate the "extra" functions.
    // First compute them on a coarse grid.
    std::vector<double> qq;
    const int Nsample=150;
    try {
      qq.resize(Nsample);
      X11.resize(Nsample);
      X22.resize(Nsample);
      X13.resize(Nsample);
      Y11.resize(Nsample);
      Y22.resize(Nsample);
      Y13.resize(Nsample);
      X1210.resize(Nsample);
      Y1210.resize(Nsample);
      V1.resize(Nsample);
      V3.resize(Nsample);
      TT.resize(Nsample);
      U1.resize(Nsample);
      U3.resize(Nsample);
      U220.resize(Nsample);
      U211.resize(Nsample);
      S2D.resize(Nsample);
    } catch(std::exception& e) {myexception(e);}
    delta=(qmax-qmin)/(Nsample-1);
    for (int i=0; i<Nsample; ++i) {
      double q = qmin+i*delta;
      std::vector<double> ef=calcEfuncs(q);
      qq[i]   = q;
      X11[i]  = ef[ 0];
      X22[i]  = ef[ 1];
      X13[i]  = ef[ 2];
      Y11[i]  = ef[ 3];
      Y22[i]  = ef[ 4];
      Y13[i]  = ef[ 5];
      X1210[i]= ef[ 6];
      Y1210[i]= ef[ 7];
      V1[i]   = ef[ 8];
      V3[i]   = ef[ 9];
      TT[i]   = ef[10];
      U1[i]   = ef[11];
      U3[i]   = ef[12];
      U220[i] = ef[13];
      U211[i] = ef[14];
      S2D[i]  = ef[15];
    }
    // then fit splines and retabulate it onto a finer grid.
    Spline X11Spline(qq,X11);
    Spline X22Spline(qq,X22);
    Spline X13Spline(qq,X13);
    Spline Y11Spline(qq,Y11);
    Spline Y22Spline(qq,Y22);
    Spline Y13Spline(qq,Y13);
    Spline X1210Spline(qq,X1210);
    Spline Y1210Spline(qq,Y1210);
    Spline V1Spline(qq,V1);
    Spline V3Spline(qq,V3);
    Spline TTSpline(qq,TT);
    Spline U1Spline(qq,U1);
    Spline U3Spline(qq,U3);
    Spline U220Spline(qq,U220);
    Spline U211Spline(qq,U211);
    Spline S2DSpline(qq,S2D);
    try {
      X11.resize(NqTable);
      X22.resize(NqTable);
      X13.resize(NqTable);
      Y11.resize(NqTable);
      Y22.resize(NqTable);
      Y13.resize(NqTable);
      X1210.resize(NqTable);
      Y1210.resize(NqTable);
      V1.resize(NqTable);
      V3.resize(NqTable);
      TT.resize(NqTable);
      U1.resize(NqTable);
      U3.resize(NqTable);
      U220.resize(NqTable);
      U211.resize(NqTable);
      S2D.resize(NqTable);
    } catch(std::exception& e) {myexception(e);}
    delta=(qmax-qmin)/(NqTable-1);
    for (int i=0; i<NqTable; ++i) {
      double q  = qmin+i*delta;
      X11[i]  = X11Spline(q);
      X22[i]  = X22Spline(q);
      X13[i]  = X13Spline(q);
      Y11[i]  = Y11Spline(q);
      Y22[i]  = Y22Spline(q);
      Y13[i]  = Y13Spline(q);
      X1210[i]= X1210Spline(q);
      Y1210[i]= Y1210Spline(q);
      V1[i]   = V1Spline(q);
      V3[i]   = V3Spline(q);
      TT[i]   = TTSpline(q);
      U1[i]   = U1Spline(q);
      U3[i]   = U3Spline(q);
      U220[i] = U220Spline(q);
      U211[i] = U211Spline(q);
      S2D[i]  = S2DSpline(q);
    }
}

void
LSM::writeSaveFile(const char fname[]) {
    // Save the XX, YY, etc. arrays to a file.
    std::ofstream fs(fname,std::ios::trunc);
    if (!fs) {
      std::cerr<<"Unable to open "<<fname<<" for writing."<<std::endl;
      myexit(1);
    }
    for (int i=0; i<NqTable; ++i)

      fs << std::scientific << std::setw(20) << std::setprecision(9) << X11[i]
         << std::scientific << std::setw(20) << std::setprecision(9) << X22[i]
         << std::scientific << std::setw(20) << std::setprecision(9) << X13[i]
         << std::scientific << std::setw(20) << std::setprecision(9) << Y11[i]
         << std::scientific << std::setw(20) << std::setprecision(9) << Y22[i]
         << std::scientific << std::setw(20) << std::setprecision(9) << Y13[i]
         << std::scientific << std::setw(20) << std::setprecision(9) << X1210[i]
         << std::scientific << std::setw(20) << std::setprecision(9) << Y1210[i]
         << std::scientific << std::setw(20) << std::setprecision(9) << V1[i]
         << std::scientific << std::setw(20) << std::setprecision(9) << V3[i]
         << std::scientific << std::setw(20) << std::setprecision(9) << TT[i]
         << std::scientific << std::setw(20) << std::setprecision(9) << U1[i]
         << std::scientific << std::setw(20) << std::setprecision(9) << U3[i]
         << std::scientific << std::setw(20) << std::setprecision(9) << U220[i]
         << std::scientific << std::setw(20) << std::setprecision(9) << U211[i]
         << std::scientific << std::setw(20) << std::setprecision(9) << S2D[i]
         << std::endl;
    fs.close();

}

void
LSM::readSaveFile(const char fname[]) {
    // Read the XX, YY, etc. arrays from a file.
    try {
      X11.resize(NqTable);
      X22.resize(NqTable);
      X13.resize(NqTable);
      Y11.resize(NqTable);
      Y22.resize(NqTable);
      Y13.resize(NqTable);
      X1210.resize(NqTable);
      Y1210.resize(NqTable);
      V1.resize(NqTable);
      V3.resize(NqTable);
      TT.resize(NqTable);
      U1.resize(NqTable);
      U3.resize(NqTable);
      U220.resize(NqTable);
      U211.resize(NqTable);
      S2D.resize(NqTable);
    } catch(std::exception& e) {myexception(e);}
    std::ifstream fs(fname);
    if (!fs) {
      std::cerr<<"Unable to open "<<fname<<" for reading."<<std::endl;
      myexit(1);
    }
    for (int i=0; i<NqTable; ++i) {
      std::string ss;
      getline(fs,ss);
      if (fs.fail()) {
        std::cerr<<"Error reading line "<<i<<" of "<<fname<<std::endl;
      }
      std::istringstream(ss) >> X11[i] >> X22[i] >> X13[i]
                             >> Y11[i] >> Y22[i] >> Y13[i]
                             >> X1210[i] >> Y1210[i]
                             >> V1[i] >> V3[i] >> TT[i] >> U1[i] >> U3[i]
                             >> U220[i] >> U211[i] >> S2D[i];
    }
    fs.close();
}

std::vector<double>
LSM::interpEfuncs(const double q) {
    // Does a linear interpolation to return the "extra" functions.
    std::vector<double> ef(16);
    int k=(NqTable-1)*(q-qmin)/(qmax-qmin);
    if (q>qmin && q<qmax) {
      double dq = (q-(qmin+k*delta))/delta;
      ef[ 0]=X11[k]+dq*(X11[k+1]-X11[k]);
      ef[ 1]=X22[k]+dq*(X22[k+1]-X22[k]);
      ef[ 2]=X13[k]+dq*(X13[k+1]-X13[k]);
      ef[ 3]=Y11[k]+dq*(Y11[k+1]-Y11[k]);
      ef[ 4]=Y22[k]+dq*(Y22[k+1]-Y22[k]);
      ef[ 5]=Y13[k]+dq*(Y13[k+1]-Y13[k]);
      ef[ 6]=X1210[k]+dq*(X1210[k+1]-X1210[k]);
      ef[ 7]=Y1210[k]+dq*(Y1210[k+1]-Y1210[k]);
      ef[ 8]=V1[k]+dq*(V1[k+1]-V1[k]);
      ef[ 9]=V3[k]+dq*(V3[k+1]-V3[k]);
      ef[10]=TT[k]+dq*(TT[k+1]-TT[k]);
      ef[11]=U1[k]+dq*(U1[k+1]-U1[k]);
      ef[12]=U3[k]+dq*(U3[k+1]-U3[k]);
      ef[13]=U220[k]+dq*(U220[k+1]-U220[k]);
      ef[14]=U211[k]+dq*(U211[k+1]-U211[k]);
      ef[15]=S2D[k]+dq*(S2D[k+1]-S2D[k]);
    }
    else {
      if (q<qmin) {
        ef[0]=ef[1]=ef[2]=ef[3]=ef[4]=ef[5]=ef[6]=ef[7]=ef[8]=ef[9]
             =ef[10]=ef[11]=ef[12]=ef[13]=ef[14]=ef[15]=0;
      }
      if (q>qmax) {
        ef[ 0]=X11[NqTable-1];
        ef[ 1]=X22[NqTable-1];
        ef[ 2]=X13[NqTable-1];
        ef[ 3]=Y11[NqTable-1];
        ef[ 4]=Y22[NqTable-1];
        ef[ 5]=Y13[NqTable-1];
        ef[ 6]=X1210[NqTable-1];
        ef[ 7]=Y1210[NqTable-1];
        ef[ 8]=V1[NqTable-1];
        ef[ 9]=V3[NqTable-1];
        ef[10]=TT[NqTable-1];
        ef[11]=U1[NqTable-1];
        ef[12]=U3[NqTable-1];
        ef[13]=U220[NqTable-1];
        ef[14]=U211[NqTable-1];
        ef[15]=S2D[NqTable-1];
      }
    }
    return(ef);
}
//
std::vector<double>
LSM::calcEfuncsMG(const double q) {
    // Computes the time derivatives of the "extra" functions of Q which come in beyond the
    // Zel'dovich approximation. These will be necessary for claulcating direct RSD in MG. Added by George. ' will indicate time derivatives.
    const int Nk=kLin.size();
    const double kmax=exp(kLin[Nk-1]);
    int Nint=(int)(8*kmax*q+512);
    if (Nint>=20000) Nint=20000;
    const double hh=(kLin[Nk-1]-kLin[0])/Nint;
    double sum0=0,sum1=0,sum2=0,sum3=0,sum4=0,sum5=0,sum6=0,sum7=0,sumS=0,sum13=0,sum14=0,sum15=0,
           sum8=0,sum9=0,sum10=0,sum11=0,sum12=0,sum16=0,sum17=0,sum18=0,sum19=0,sum20=0,sum21=0,sum22=0,sum23=0,sum24=0,sum25=0,sum26=0,sum27=0,sum28=0,sum29=0,sum30=0,
           sum31=0,sum32=0,sum33=0,sum34=0,sum35=0,sum36=0,sum37=0,sum38=0,sum39=0,sum40=0, sumS11p2=0, sumS1p1p2=0, sumS11p2p=0,sumS1p1p2p=0;
#pragma omp parallel for reduction(+:sum0,sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,sum9,sum10,sum11,sum12,sumS,sum13,sum14,sum15,sum16,sum17,sum18,sum19,sum20,sum21,sum22,sum23,sum24,sum25,sum26,sum27,sum28,sum29,sum30,sum31,sum32,sum33,sum34,sum35,sum36,sum37)
//#pragma omp parallel for reduction(+:sum0,sum1,sum2,sum3,sum4,sum5,sum6,sum7)
    for (int i=1; i<Nint; ++i) {
      double xx = kLin[0]+i*hh;
      double kk = exp(xx);
      double k2 = kk*kk;
      double kq = kk*q;
      //std::cout << kk << std::endl;
      double R1 = R1spl(kk);
      double R2 = R2spl(kk);
      double Q1 = Q1spl(kk);
      double Q2 = Q2spl(kk);
      double Q5 = Q5spl(kk);
      double Q8 = Q8spl(kk);
      double Qs = Qsspl(kk);
      double QI = QIspl(kk);
      double R12 = R12spl(kk);
      double RI = RIspl(kk);
      double R1prime = R1primespl(kk);
      double R2prime = R2primespl(kk);
      double Q1prime = Q1primespl(kk);
      double Q2prime = Q2primespl(kk);
      double Q5prime = Q5primespl(kk);
      double Q8prime = Q8primespl(kk);
      double Qsprime = Qprimesspl(kk);
      double QIprime = QIprimespl(kk);
      double R12prime = R12primespl(kk);
      double RIprime = RIprimespl(kk);
      double Q1halfprime = Q1halfprimespl(kk);
      //growth factor in MG
      double fgrowth = fgrowthspl(kk);
      //std::cout << "hey" << std::endl;
      //extra Q splines to store
     // double Q3 = Q3spl(kk);
     // double Q4 = Q4spl(kk);
     // double Q6 = Q6spl(kk);
     // double Q7 = Q7spl(kk);
     // double Q9 = Q9spl(kk);
     // double Q10 = Q10spl(kk);
     // double Q11 = Q11spl(kk);
     // double Q12 = Q12spl(kk);
     // double Q13 = Q13spl(kk);
    
      double ap = cos(M_PI/2.*kk/kmax);
      std::vector<double> jl=sphBess(kk*q);
      double j1 = kq*jl[1];
      double j2,j3;
      if (kq<0.1) {
        j2 = pow(kq,2.0)/15.  - pow(kq,4.0)/210.;
        j3 = pow(kq,3.0)/105. - pow(kq,5.0)/1890.;
      }
      else {
        j2 = 3.*jl[1]-jl[0];
        j3 = 5.*j2/(kq)-(kq)*jl[1];
      }
      int wt= 2+2*(i%2);
      //Changes by George to account for MG correlators
      sum0 += kk*wt*(9./98.*Q1halfprime*(2./3.-2*jl[1]));		// X^{(22')}
      sum1 += kk*wt*(5./21.*R1prime*(2./3.-2*jl[1]));		// X^{(13')}
      sum2 += kk*wt*(9./98.*Q1halfprime*(-2*jl[0]+6*jl[1]));		// Y^{(22')}
      sum3 += kk*wt*(5./21.*R1prime*(-2*jl[0]+6*jl[1]));		// Y^{(13')}
      sum4 += kk*wt*(2*(RIprime-R2prime)+3*RIprime*jl[0]-3*(RIprime+2*R2prime+2*R12prime+2*Q5prime)*jl[1])/14.;//X12'10
      sum5 += kk*wt*(RIprime+2*R2prime+2*R12prime+2*Q5prime)*(jl[0]-3*jl[1])*(-3./14.);// Y_{10}^{(12')}
      sum6 +=    wt*(RIprime*j1)*(-3./7.);			// V_1^{(112')}
      sum7 +=    wt*(QIprime*j1)*(-3./7.);			// V_3^{(112')}
      sumS +=    wt*(2*RIprime+4*R2prime+QIprime+2*Q2prime)*(3./7.*j2/(kk*q));	// S^{(112')}
      sum8 +=    wt*(2*RIprime+4*R2prime+QIprime+2*Q2prime)*j3*(-3./7.);	// T^{(112')}
      sum9 += k2*wt*(R1prime*j1)*(-5./21.);			// U^{(3')}
      sum10+= k2*wt*(Q8prime*j1)*(-3./7.);			// U_{20}^{(2')}
      sum11+= k2*wt*((R12prime)*j1)*(-6./7.);		// U_{11}^{(2')}
      sum12+= k2*wt*(Qsprime*j1)*(-2./7.)*ap;		// Shear term

      sum16 += kk*wt*(9./98.*Q1prime*(2./3.-2*jl[1]));		// X^{(2'2')}
      sum17 += kk*wt*fgrowth*(5./21.*R1*(2./3.-2*jl[1]));		// X^{(1'3)}
      sum18 += kk*wt*(9./98.*Q1prime*(-2*jl[0]+6*jl[1]));		// Y^{(2'2')}
      sum19 += kk*wt*fgrowth*(5./21.*R1*(-2*jl[0]+6*jl[1]));		// Y^{(1'3)}
      sum20 += kk*wt*fgrowth*(5./21.*R1prime*(2./3.-2*jl[1]));		// X^{(1'3')}
      sum21 += kk*wt*fgrowth*(5./21.*R1prime*(-2*jl[0]+6*jl[1]));		// Y^{(1'3')}

      sum22 += kk*wt*fgrowth*(2*(RI-R2)+3*RI*jl[0]-3*(RI+2*R2+2*R12+2*Q5)*jl[1])/14.;//X1'210
      sum23 += kk*wt*fgrowth*(RI+2*R2+2*R12+2*Q5)*(jl[0]-3*jl[1])*(-3./14.);// Y_{10}^{(1'2)}
      sum24 += kk*wt*fgrowth*(2*(RIprime-R2prime)+3*RIprime*jl[0]-3*(RIprime+2*R2prime+2*R12prime+2*Q5prime)*jl[1])/14.;//X1'2'_{10}
      sum25 += kk*wt*fgrowth*(RIprime+2*R2prime+2*R12prime+2*Q5prime)*(jl[0]-3*jl[1])*(-3./14.);// Y_{10}^{(1'2')}

      sum26 +=    wt*fgrowth*(RI*j1)*(-3./7.);			// V_1^{(11'2)}
      sum27 +=    wt*fgrowth*(QI*j1)*(-3./7.);			// V_3^{(11'2)}
      sumS11p2 +=  wt*fgrowth*(2*RI+4*R2+QI+2*Q2)*(3./7.*j2/(kk*q));	// S^{(11'2)}
      sum28 +=    wt*fgrowth*(2*RI+4*R2+QI+2*Q2)*j3*(-3./7.);	// T^{(11'2)}

      sum29 +=    wt*fgrowth*fgrowth*(RI*j1)*(-3./7.);			// V_1^{(1'1'2)}
      sum30 +=    wt*fgrowth*fgrowth*(QI*j1)*(-3./7.);			// V_3^{(1'1'2)}
      sumS1p1p2 +=  wt*fgrowth*fgrowth*(2*RI+4*R2+QI+2*Q2)*(3./7.*j2/(kk*q));	// S^{(1'1'2)}
      sum31 +=    wt*fgrowth*fgrowth*(2*RI+4*R2+QI+2*Q2)*j3*(-3./7.);	// T^{(1'1'2)}

      sum32 +=    wt*fgrowth*(RIprime*j1)*(-3./7.);			// V_1^{(11'2')}
      sum33 +=    wt*fgrowth*(QIprime*j1)*(-3./7.);			// V_3^{(11'2')}
      sumS11p2p +=    wt*fgrowth*(2*RIprime+4*R2prime+QIprime+2*Q2prime)*(3./7.*j2/(kk*q));	// S^{(11'2')}
      sum34 +=    wt*fgrowth*(2*RIprime+4*R2prime+QIprime+2*Q2prime)*j3*(-3./7.);	// T^{(11'2')}

      sum35 +=    wt*fgrowth*fgrowth*(RIprime*j1)*(-3./7.);			// V_1^{(1'1'2')}
      sum36 +=    wt*fgrowth*fgrowth*(QIprime*j1)*(-3./7.);			// V_3^{(1'1'2')}
      sumS1p1p2p +=    wt*fgrowth*fgrowth*(2*RIprime+4*R2prime+QIprime+2*Q2prime)*(3./7.*j2/(kk*q));	// S^{(1'1'2')}
      sum37 +=    wt*fgrowth*fgrowth*(2*RIprime+4*R2prime+QIprime+2*Q2prime)*j3*(-3./7.);	// T^{(1'1'2')}
    }
    //files.close();
    sum6 += sumS;
    sum7 += sumS;
    sum26 += sumS11p2;
    sum27 += sumS11p2;
    sum29 += sumS1p1p2;
    sum30 += sumS1p1p2;
    sum32 += sumS11p2p;
    sum33 += sumS11p2p;
    sum35 += sumS1p1p2p;
    sum36 += sumS1p1p2p;
    //std::vector<double> sum(16);
    std::vector<double> sum(43);
    sum[ 1] = sum0 * hh/3.0/(2*M_PI*M_PI);
    sum[ 2] = sum1 * hh/3.0/(2*M_PI*M_PI);
    sum[ 4] = sum2 * hh/3.0/(2*M_PI*M_PI);
    sum[ 5] = sum3 * hh/3.0/(2*M_PI*M_PI);
    sum[ 6] = sum4 * hh/3.0/(2*M_PI*M_PI);
    sum[ 7] = sum5 * hh/3.0/(2*M_PI*M_PI);
    sum[ 8] = sum6 * hh/3.0/(2*M_PI*M_PI);
    sum[ 9] = sum7 * hh/3.0/(2*M_PI*M_PI);	
    sum[10] = sum8 * hh/3.0/(2*M_PI*M_PI);	
    sum[12] = sum9 * hh/3.0/(2*M_PI*M_PI);	
    sum[13] = sum10* hh/3.0/(2*M_PI*M_PI);	
    sum[14] = sum11* hh/3.0/(2*M_PI*M_PI);	
    sum[15] = sum12* hh/3.0/(2*M_PI*M_PI);	
    sum[16] = sum16* hh/3.0/(2*M_PI*M_PI);	
    sum[17] = sum17* hh/3.0/(2*M_PI*M_PI);	
    sum[18] = sum18* hh/3.0/(2*M_PI*M_PI);
    sum[19] = sum19* hh/3.0/(2*M_PI*M_PI);	
    sum[20] = sum20* hh/3.0/(2*M_PI*M_PI);	
    sum[21] = sum21* hh/3.0/(2*M_PI*M_PI);
    sum[22] = sum22* hh/3.0/(2*M_PI*M_PI);
    sum[23] = sum23* hh/3.0/(2*M_PI*M_PI);
    sum[24] = sum24* hh/3.0/(2*M_PI*M_PI);
    sum[25] = sum25* hh/3.0/(2*M_PI*M_PI);
    sum[26] = sum26* hh/3.0/(2*M_PI*M_PI);
    sum[27] = sum27* hh/3.0/(2*M_PI*M_PI);
    sum[28] = sum28* hh/3.0/(2*M_PI*M_PI);
    sum[29] = sum29 * hh/3.0/(2*M_PI*M_PI);	
    sum[30] = sum30 * hh/3.0/(2*M_PI*M_PI);	
    sum[31] = sum31 * hh/3.0/(2*M_PI*M_PI);	
    sum[32] = sum32* hh/3.0/(2*M_PI*M_PI);	
    sum[33] = sum33* hh/3.0/(2*M_PI*M_PI);	
    sum[34] = sum34* hh/3.0/(2*M_PI*M_PI);	
    sum[35] = sum35* hh/3.0/(2*M_PI*M_PI);
    sum[36] = sum36* hh/3.0/(2*M_PI*M_PI);	
    sum[37] = sum37* hh/3.0/(2*M_PI*M_PI);
    // Now tabulate the pieces going as Plin.
    sum0=sum1=sum2=sum3=sum4=sum5=sum6=sum7=0;
#pragma omp parallel for reduction(+:sum0,sum1,sum2,sum3,sum4,sum5,sum6,sum7)
    for (int i=1; i<Nint; ++i) {
      double xx = kLin[0]+i*hh;
      double kk = exp(xx);
      double k2 = kk*kk;
      double kq = kk*q;
      double fgrowth = fgrowthspl(kk);
      int    jj = (int)(i*hh*dkinv);
      if (jj>=pLin.size()-2) jj=pLin.size()-2;
      double pk = exp(pLin[jj]+(xx-kLin[jj])*
                     (pLin[jj+1]-pLin[jj])/(kLin[jj+1]-kLin[jj]));
      std::vector<double> jl=sphBess(kk*q);
      double j1 = kq*jl[1];
      int wt= 2+2*(i%2);
      sum0 += kk*wt*pk*(2./3.-2*jl[1]);		// X^{(11)}
      sum1 += kk*wt*pk*(-2.*jl[0]+6*jl[1]);	// Y^{(11)}
      sum2 += k2*wt*pk*(-j1);			// U^{(1)}
      sum3 += kk*wt*pk*fgrowth*(2./3.-2*jl[1]);		// X^{(11')}
      sum4 += kk*wt*pk*fgrowth*(-2.*jl[0]+6*jl[1]);	// Y^{(11')}
      sum5 += k2*wt*pk*fgrowth*(-j1);			// U^{(1')}
      sum6 += kk*wt*pk*fgrowth*fgrowth*(2./3.-2*jl[1]);		// X^{(1'1')}
      sum7 += kk*wt*pk*fgrowth*fgrowth*(-2.*jl[0]+6*jl[1]);	// Y^{(1'1')}
    }
    sum[ 0] = sum0 * hh/3.0/(2*M_PI*M_PI);
    sum[ 3] = sum1 * hh/3.0/(2*M_PI*M_PI);
    sum[11] = sum2 * hh/3.0/(2*M_PI*M_PI);
    sum[38] = sum3 * hh/3.0/(2*M_PI*M_PI);
    sum[39] = sum4 * hh/3.0/(2*M_PI*M_PI);
    sum[40] = sum5 * hh/3.0/(2*M_PI*M_PI);
    sum[41] = sum6 * hh/3.0/(2*M_PI*M_PI);
    sum[42] = sum7 * hh/3.0/(2*M_PI*M_PI);
    //std::cout  << sum[ 0] << std::endl;
    return(sum);
}

void
LSM::tabulateEfuncsMG() {
    // Tabulate the "extra" functions.
    // First compute them on a coarse grid.
    std::vector<double> qq;
    const int Nsample=150;
    try {
      qq.resize(Nsample);
      X11MG.resize(Nsample);
      X22p.resize(Nsample);
      X13p.resize(Nsample);
      Y11MG.resize(Nsample);
      Y22p.resize(Nsample);
      Y13p.resize(Nsample);
      X12p10.resize(Nsample);
      Y12p10.resize(Nsample);
      V1112p.resize(Nsample);
      V3112p.resize(Nsample);
      TT112p.resize(Nsample);
      U1MG.resize(Nsample);
      U3p.resize(Nsample);
      U220p.resize(Nsample);
      U211p.resize(Nsample);
      S2Dp.resize(Nsample);
      X2p2p.resize(Nsample);
      X1p3.resize(Nsample);
      Y2p2p.resize(Nsample);
      Y1p3.resize(Nsample);
      X1p3p.resize(Nsample);
      Y1p3p.resize(Nsample);
      X1p210.resize(Nsample);
      Y1p210.resize(Nsample);
      X1p2p10.resize(Nsample);
      Y1p2p10.resize(Nsample);
      V111p2.resize(Nsample);
      V311p2.resize(Nsample);
      TT11p2.resize(Nsample);
      V11p1p2.resize(Nsample);
      V31p1p2.resize(Nsample);
      TT1p1p2.resize(Nsample);
      V111p2p.resize(Nsample);
      V311p2p.resize(Nsample);
      TT11p2p.resize(Nsample);
      V11p1p2p.resize(Nsample);
      V31p1p2p.resize(Nsample);
      TT1p1p2p.resize(Nsample);
      X11p.resize(Nsample);
      Y11p.resize(Nsample);
      U1p.resize(Nsample);
      X1p1p.resize(Nsample);
      Y1p1p.resize(Nsample);

    } catch(std::exception& e) {myexception(e);}
    delta=(qmax-qmin)/(Nsample-1);
    for (int i=0; i<Nsample; ++i) {
      double q = qmin+i*delta;
      std::vector<double> efMG=calcEfuncsMG(q);
      qq[i]   = q;
      X11MG[i]  = efMG[ 0];
      X22p[i]  = efMG[ 1];
      X13p[i]  = efMG[ 2];
      Y11MG[i]  = efMG[ 3];
      Y22p[i]  = efMG[ 4];
      Y13p[i]  = efMG[ 5];
      X12p10[i]= efMG[ 6];
      Y12p10[i]= efMG[ 7];
      V1112p[i]   = efMG[ 8];
      V3112p[i]   = efMG[ 9];
      TT112p[i]   = efMG[10];
      U1MG[i]   = efMG[11];
      U3p[i]   = efMG[12];
      U220p[i] = efMG[13];
      U211p[i] = efMG[14];
      S2Dp[i]  = efMG[15];

      X2p2p[i]  = efMG[16];
      X1p3[i]  = efMG[17];
      Y2p2p[i]  = efMG[18];
      Y1p3[i]  = efMG[19];
      X1p3p[i]  = efMG[20];
      Y1p3p[i]  = efMG[21];

      X1p210[i]= efMG[22];
      Y1p210[i]= efMG[23];
      X1p2p10[i]= efMG[24];
      Y1p2p10[i]= efMG[25];

      V111p2[i]   = efMG[26];
      V311p2[i]   = efMG[27];
      TT11p2[i]   = efMG[28];

      V11p1p2[i]   = efMG[29];
      V31p1p2[i]   = efMG[30];
      TT1p1p2[i]   = efMG[31];

      V111p2p[i]   = efMG[32];
      V311p2p[i]   = efMG[33];
      TT11p2p[i]   = efMG[34];

      V11p1p2p[i]   = efMG[35];
      V31p1p2p[i]   = efMG[36];
      TT1p1p2p[i]   = efMG[37];

      X11p[i]   = efMG[38];
      Y11p[i]   = efMG[39];
      U1p[i]   = efMG[40];

      X1p1p[i]   = efMG[41];
      Y1p1p[i]   = efMG[42];
    }
    // then fit splines and retabulate it onto a finer grid.
    Spline X11MGSpline(qq,X11MG);
    Spline X22pSpline(qq,X22p);
    Spline X13pSpline(qq,X13p);
    Spline Y11MGSpline(qq,Y11MG);
    Spline Y22pSpline(qq,Y22p);
    Spline Y13pSpline(qq,Y13p);
    Spline X12p10Spline(qq,X12p10);
    Spline Y12p10Spline(qq,Y12p10);
    Spline V1112pSpline(qq,V1112p);
    Spline V3112pSpline(qq,V3112p);
    Spline TT112pSpline(qq,TT112p);
    Spline U1MGSpline(qq,U1MG);
    Spline U3pSpline(qq,U3p);
    Spline U220pSpline(qq,U220p);
    Spline U211pSpline(qq,U211p);
    Spline S2DpSpline(qq,S2Dp);

    Spline X2p2pSpline(qq,X2p2p);
    Spline X1p3Spline(qq,X1p3);
    Spline Y2p2pSpline(qq,Y2p2p);
    Spline Y1p3Spline(qq,Y1p3);
    Spline X1p3pSpline(qq,X1p3p);
    Spline Y1p3pSpline(qq,Y1p3p);
    Spline X1p210Spline(qq,X1p210);
    Spline Y1p210Spline(qq,Y1p210);
    Spline X1p2p10Spline(qq,X1p2p10);
    Spline Y1p2p10Spline(qq,Y1p2p10);
    Spline V111p2Spline(qq,V111p2);
    Spline V311p2Spline(qq,V311p2);
    Spline TT11p2Spline(qq,TT11p2);
    Spline V11p1p2Spline(qq,V11p1p2);
    Spline V31p1p2Spline(qq,V31p1p2);
    Spline TT1p1p2Spline(qq,TT1p1p2);
    Spline V111p2pSpline(qq,V111p2p);
    Spline V311p2pSpline(qq,V311p2p);
    Spline TT11p2pSpline(qq,TT11p2p);
    Spline V11p1p2pSpline(qq,V11p1p2p);
    Spline V31p1p2pSpline(qq,V31p1p2p);
    Spline TT1p1p2pSpline(qq,TT1p1p2p);
    Spline X11pSpline(qq,X11p);
    Spline Y11pSpline(qq,Y11p);
    Spline U1pSpline(qq,U1p);
    Spline X1p1pSpline(qq,X1p1p);
    Spline Y1p1pSpline(qq,Y1p1p);

    try {
      X11MG.resize(NqTable);
      X22p.resize(NqTable);
      X13p.resize(NqTable);
      Y11MG.resize(NqTable);
      Y22p.resize(NqTable);
      Y13p.resize(NqTable);
      X12p10.resize(NqTable);
      Y12p10.resize(NqTable);
      V1112p.resize(NqTable);
      V3112p.resize(NqTable);
      TT112p.resize(NqTable);
      U1MG.resize(NqTable);
      U3p.resize(NqTable);
      U220p.resize(NqTable);
      U211p.resize(NqTable);
      S2Dp.resize(NqTable);
      X2p2p.resize(NqTable);
      X1p3.resize(NqTable);
      Y2p2p.resize(NqTable);
      Y1p3.resize(NqTable);
      X1p3p.resize(NqTable);
      Y1p3p.resize(NqTable);
      X1p210.resize(NqTable);
      Y1p210.resize(NqTable);
      X1p2p10.resize(NqTable);
      Y1p2p10.resize(NqTable);
      V111p2.resize(NqTable);
      V311p2.resize(NqTable);
      TT11p2.resize(NqTable);
      V11p1p2.resize(NqTable);
      V31p1p2.resize(NqTable);
      TT1p1p2.resize(NqTable);
      V111p2p.resize(NqTable);
      V311p2p.resize(NqTable);
      TT11p2p.resize(NqTable);
      V11p1p2p.resize(NqTable);
      V31p1p2p.resize(NqTable);
      TT1p1p2p.resize(NqTable);
      X11p.resize(NqTable);
      Y11p.resize(NqTable);
      U1p.resize(NqTable);
      X1p1p.resize(NqTable);
      Y1p1p.resize(NqTable);

    } catch(std::exception& e) {myexception(e);}
    delta=(qmax-qmin)/(NqTable-1);
    for (int i=0; i<NqTable; ++i) {
      double q  = qmin+i*delta;
      X11MG[i]  = X11MGSpline(q);
      X22p[i]  = X22pSpline(q);
      X13p[i]  = X13pSpline(q);
      Y11MG[i]  = Y11MGSpline(q);
      Y22p[i]  = Y22pSpline(q);
      Y13p[i]  = Y13pSpline(q);
      X12p10[i]= X12p10Spline(q);
      Y12p10[i]= Y12p10Spline(q);
      V1112p[i]   = V1112pSpline(q);
      V3112p[i]   = V3112pSpline(q);
      TT112p[i]   = TT112pSpline(q);
      U1MG[i]   = U1MGSpline(q);
      U3p[i]   = U3pSpline(q);
      U220p[i] = U220pSpline(q);
      U211p[i] = U211pSpline(q);
      S2Dp[i]  = S2DpSpline(q);

      X2p2p[i]  = X2p2pSpline(q); 
      X1p3[i]  = X1p3Spline(q);
      Y2p2p[i]  = Y2p2pSpline(q);
      Y1p3[i]  = Y1p3Spline(q);
      X1p3p[i]  = X1p3pSpline(q);
      Y1p3p[i]  = Y1p3pSpline(q);
      X1p210[i]= X1p210Spline(q);
      Y1p210[i]= Y1p210Spline(q);
      X1p2p10[i]= X1p2p10Spline(q);
      Y1p2p10[i]= Y1p2p10Spline(q);
      V111p2[i]   = V111p2Spline(q);
      V311p2[i]   = V311p2Spline(q);
      TT11p2[i]   = TT11p2Spline(q);
      V11p1p2[i]   = V11p1p2Spline(q);
      V31p1p2[i]   = V31p1p2Spline(q);
      TT1p1p2[i]   = TT1p1p2Spline(q);
      V111p2p[i]   = V111p2pSpline(q);
      V311p2p[i]   = V311p2pSpline(q);
      TT11p2p[i]   = TT11p2pSpline(q);
      V11p1p2p[i]   = V11p1p2pSpline(q);
      V31p1p2p[i]   = V31p1p2pSpline(q);
      TT1p1p2p[i]   = TT1p1p2pSpline(q);
      X11p[i]   = X11pSpline(q);
      Y11p[i]   = Y11pSpline(q);
      U1p[i]   = U1pSpline(q);
      X1p1p[i]   = X1p1pSpline(q);
      Y1p1p[i]   = Y1p1pSpline(q);
    }
}

std::vector<double>
LSM::interpEfuncsMG(const double q) {
    // Does a linear interpolation to return the "extra" functions.
    std::vector<double> efMG(43);
    int k=(NqTable-1)*(q-qmin)/(qmax-qmin);
    if (q>qmin && q<qmax) {
      double dq = (q-(qmin+k*delta))/delta;
      //std::cout  << X22p[k] << std::endl;
      efMG[ 0]=X11MG[k]+dq*(X11MG[k+1]-X11MG[k]);
      efMG[ 1]=X22p[k]+dq*(X22p[k+1]-X22p[k]);
      efMG[ 2]=X13p[k]+dq*(X13p[k+1]-X13p[k]);
      efMG[ 3]=Y11MG[k]+dq*(Y11MG[k+1]-Y11MG[k]);
      efMG[ 4]=Y22p[k]+dq*(Y22p[k+1]-Y22p[k]);
      efMG[ 5]=Y13p[k]+dq*(Y13p[k+1]-Y13p[k]);
      efMG[ 6]=X12p10[k]+dq*(X12p10[k+1]-X12p10[k]);
      efMG[ 7]=Y12p10[k]+dq*(Y12p10[k+1]-Y12p10[k]);
      efMG[ 8]=V1112p[k]+dq*(V1112p[k+1]-V1112p[k]);
      efMG[ 9]=V3112p[k]+dq*(V3112p[k+1]-V3112p[k]);
      efMG[10]=TT112p[k]+dq*(TT112p[k+1]-TT112p[k]);
      efMG[11]=U1MG[k]+dq*(U1MG[k+1]-U1MG[k]);
      efMG[12]=U3p[k]+dq*(U3p[k+1]-U3p[k]);
      efMG[13]=U220p[k]+dq*(U220p[k+1]-U220p[k]);
      efMG[14]=U211p[k]+dq*(U211p[k+1]-U211p[k]);
      efMG[15]=S2Dp[k]+dq*(S2Dp[k+1]-S2Dp[k]);

      efMG[16]=X2p2p[k]+dq*(X2p2p[k+1]-X2p2p[k]);
      efMG[17]=X1p3[k]+dq*(X1p3[k+1]-X1p3[k]);
      efMG[18]=Y2p2p[k]+dq*(Y2p2p[k+1]-Y2p2p[k]);
      efMG[19]=Y1p3[k]+dq*(Y1p3[k+1]-Y1p3[k]);
      efMG[20]=X1p3p[k]+dq*(X1p3p[k+1]-X1p3p[k]);
      efMG[21]=Y1p3p[k]+dq*(Y1p3p[k+1]-Y1p3p[k]);
      efMG[22]=X1p210[k]+dq*(X1p210[k+1]-X1p210[k]);
      efMG[23]=Y1p210[k]+dq*(Y1p210[k+1]-Y1p210[k]);
      efMG[24]=X1p2p10[k]+dq*(X1p2p10[k+1]-X1p2p10[k]);
      efMG[25]=Y1p2p10[k]+dq*(Y1p2p10[k+1]-Y1p2p10[k]);
      efMG[26]=V111p2[k]+dq*(V111p2[k+1]-V111p2[k]);
      efMG[27]=V311p2[k]+dq*(V311p2[k+1]-V311p2[k]);
      efMG[28]=TT11p2[k]+dq*(TT11p2[k+1]-TT11p2[k]);
      efMG[29]=V11p1p2[k]+dq*(V11p1p2[k+1]-V11p1p2[k]);
      efMG[30]=V31p1p2[k]+dq*(V31p1p2[k+1]-V31p1p2[k]);
      efMG[31]=TT1p1p2[k]+dq*(TT1p1p2[k+1]-TT1p1p2[k]);
      efMG[32]=V111p2p[k]+dq*(V111p2p[k+1]-V111p2p[k]);
      efMG[33]=V311p2p[k]+dq*(V311p2p[k+1]-V311p2p[k]);
      efMG[34]=TT11p2p[k]+dq*(TT11p2p[k+1]-TT11p2p[k]);
      efMG[35]=V11p1p2p[k]+dq*(V11p1p2p[k+1]-V11p1p2p[k]);
      efMG[36]=V31p1p2p[k]+dq*(V31p1p2p[k+1]-V31p1p2p[k]);
      efMG[37]=TT1p1p2p[k]+dq*(TT1p1p2p[k+1]-TT1p1p2p[k]);
      efMG[38]=X11p[k]+dq*(X11p[k+1]-X11p[k]);
      efMG[39]=Y11p[k]+dq*(Y11p[k+1]-Y11p[k]);
      efMG[40]=U1p[k]+dq*(U1p[k+1]-U1p[k]);
      efMG[41]=X1p1p[k]+dq*(X1p1p[k+1]-X1p1p[k]);
      efMG[42]=Y1p1p[k]+dq*(Y1p1p[k+1]-Y1p1p[k]);
    }
    else {
      if (q<qmin) {
        efMG[0]=efMG[1]=efMG[2]=efMG[3]=efMG[4]=efMG[5]=efMG[6]=efMG[7]=efMG[8]=efMG[9]
             =efMG[10]=efMG[11]=efMG[12]=efMG[13]=efMG[14]=efMG[15]=efMG[16]=efMG[17]=efMG[18]=efMG[19]=efMG[14]=efMG[20]
             =efMG[21]=efMG[22]=efMG[23]=efMG[24]=efMG[25]=efMG[26]=efMG[27]=efMG[28]=efMG[29]=efMG[30]=efMG[31]=efMG[32]
             =efMG[34]=efMG[35]=efMG[36]=efMG[37]=efMG[38]=efMG[39]=efMG[40]=efMG[41]=efMG[42]=0;
      }
      if (q>qmax) {
      efMG[ 0]=X11MG[NqTable-1];
      efMG[ 1]=X22p[NqTable-1];
      efMG[ 2]=X13p[NqTable-1];
      efMG[ 3]=Y11MG[NqTable-1];
      efMG[ 4]=Y22p[NqTable-1];
      efMG[ 5]=Y13p[NqTable-1];
      efMG[ 6]=X12p10[NqTable-1];
      efMG[ 7]=Y12p10[NqTable-1];
      efMG[ 8]=V1112p[NqTable-1];
      efMG[ 9]=V3112p[NqTable-1];
      efMG[10]=TT112p[NqTable-1];
      efMG[11]=U1MG[NqTable-1];
      efMG[12]=U3p[NqTable-1];
      efMG[13]=U220p[NqTable-1];
      efMG[14]=U211p[NqTable-1];
      efMG[15]=S2Dp[NqTable-1];

      efMG[16]=X2p2p[NqTable-1];
      efMG[17]=X1p3[NqTable-1];
      efMG[18]=Y2p2p[NqTable-1];
      efMG[19]=Y1p3[NqTable-1];
      efMG[20]=X1p3p[NqTable-1];
      efMG[21]=Y1p3p[NqTable-1];
      efMG[22]=X1p210[NqTable-1];
      efMG[23]=Y1p210[NqTable-1];
      efMG[24]=X1p2p10[NqTable-1];
      efMG[25]=Y1p2p10[NqTable-1];
      efMG[26]=V111p2[NqTable-1];
      efMG[27]=V311p2[NqTable-1];
      efMG[28]=TT11p2[NqTable-1];
      efMG[29]=V11p1p2[NqTable-1];
      efMG[30]=V31p1p2[NqTable-1];
      efMG[31]=TT1p1p2[NqTable-1];
      efMG[32]=V111p2p[NqTable-1];
      efMG[33]=V311p2p[NqTable-1];
      efMG[34]=TT11p2p[NqTable-1];
      efMG[35]=V11p1p2p[NqTable-1];
      efMG[36]=V31p1p2p[NqTable-1];
      efMG[37]=TT1p1p2p[NqTable-1];
      efMG[38]=X11p[NqTable-1];
      efMG[39]=Y11p[NqTable-1];
      efMG[40]=U1p[NqTable-1];
      efMG[41]=X1p1p[NqTable-1];
      efMG[42]=Y1p1p[NqTable-1];
      }
    }
    return(efMG);
}

//
//void
//LSM::setupQR() {
void
LSM::setupQR(const char fname[]) {
    // Set up the Q and R's that we need, apodized.
    const int NkTemp=500;
    std::vector<double> ka(NkTemp),R1(NkTemp),R2(NkTemp),kaR(NkTemp);
    std::vector<double> Q1(NkTemp),Q2(NkTemp),Q5(NkTemp),Q8(NkTemp),Q3(NkTemp);
    std::vector<double> Qs(NkTemp),QI(NkTemp),R12(NkTemp),RI(NkTemp);
    //extra Q's to store
    std::vector<double> Q4(NkTemp),Q6(NkTemp),Q7(NkTemp),Q9(NkTemp);
    std::vector<double> Q10(NkTemp),Q11(NkTemp),Q12(NkTemp),Q13(NkTemp);
//#pragma omp parallel for
    //for (int i=1; i<NkTemp; ++i) {
      //std::cout << "Num threads=" << omp_get_num_threads() << std::endl;
      //double kk = exp( kLin[0]+i*(kLin[kLin.size()-1]-kLin[0])/(NkTemp-1) );
      //double ap = cos(M_PI/2.*kk/exp(kLin[kLin.size()-1]));
      //std::vector<double> Qn=lpt.Qn(kk);
      //std::vector<double> Rn=lpt.Rn(kk);
    //  //std::cout << "threads=" << omp_get_num_threads() << std::endl;
      //ka[i]=kk; R1[i]=Rn[0]*ap; R2[i]=Rn[1]*ap;
      //Q1[i]=Qn[1]*ap; Q2[i]=Qn[2]*ap; Q5[i]=Qn[5]*ap; Q8[i]=Qn[8]*ap;
      //Qs[i]=Qn[14]*ap;
      //Q3[i]=Qn[3]*ap; Q4[i]=Qn[4]*ap; Q6[i]=Qn[6]*ap; Q7[i]=Qn[7]*ap;
      //Q9[i]=Qn[9]*ap; Q10[i]=Qn[10]*ap; Q11[i]=Qn[11]*ap; Q12[i]=Qn[12]*ap; Q13[i]=Qn[13]*ap;
    //}
    // and fit splines to them.
  //  std::cout << NkTemp << std::endl;
  //  R1spl.init(ka,R1);
  //  R2spl.init(ka,R2);
  //  Q1spl.init(ka,Q1);
  //  Q2spl.init(ka,Q2);
  //  Q5spl.init(ka,Q5);
  //  Q8spl.init(ka,Q8);
  // Qsspl.init(ka,Qs);
  //  std::cout << NkTemp << std::endl;
    //extra Q's to store
    //Q3spl.init(ka,Q3);
    //Q4spl.init(ka,Q4);
    //Q6spl.init(ka,Q6);
    //Q7spl.init(ka,Q7);
    //Q9spl.init(ka,Q9);
    //Q10spl.init(ka,Q10);
    //Q11spl.init(ka,Q11);
    //Q12spl.init(ka,Q12);
    //Q13spl.init(ka,Q13);

    std::ifstream fsQ("./plin_Fr6z05wmap9_cleftQnew_z000.txt");
    //std::ifstream fsQ("./ps_python3/plin_Fr5z05wmap9_cleftQnew_z000.txt");
    //std::ifstream fsQ("./ps_python3/plin_N1z05wmap9_cleftQ0_z000.txt");
    //std::ifstream fsQ("./Q_zeros.txt");
    //std::ifstream fsQ("./Qfuncs.txt");
    //std::ifstream fsQ("./ps_python3/plin_N5z05wmap9_cleftQnew_z000.txt");
    //std::ifstream fsQ("./ps_python3/plin_z05wmap9_cleftQnew_z000.txt");
    //std::ifstream fsQ("./ps_python3/plin_Fr5z1planck_cleftQnew_z000.txt");
    for (int i=0; i<NkTemp; ++i) {
      std::string ssQ;
      getline(fsQ,ssQ);
               //std::istringstream(ssQ) >> ka[i] >> Q1[i] >> Q2[i] >> Q3[i] >> Q5[i] >> Q8[i] >> Qs[i];  
               std::istringstream(ssQ) >> ka[i] >> Q1[i] >> Q2[i] >> Q3[i] >> Q5[i] >> Q8[i] >> Qs[i] >> QI[i] >> Q7[i] >> Q9[i] >> Q11[i] >> Q12[i] >> Q13[i];                          
               //std::cout << QI[i] << std::endl;
    }
    fsQ.close();
    // std::cout << NkTemp << std::endl;
    Q1spl.init(ka,Q1);
    Q2spl.init(ka,Q2);
    Q5spl.init(ka,Q5);
    Q8spl.init(ka,Q8);
    Qsspl.init(ka,Qs);
    QIspl.init(ka,QI);
    
    std::ifstream fsR("./plin_Fr6z05wmap9_cleftRnew_z000.txt");
    //std::ifstream fsR("./ps_python3/plin_Fr5z05wmap9_cleftRnew_z000.txt");
    //std::ifstream fsR("./ps_python3/plin_N5z05wmap9_cleftRnew_z000.txt");
    //std::ifstream fsR("./ps_python3/plin_N1z05wmap9_cleftR0_z000.txt");
    //std::ifstream fsR("./R_zeros.txt");
    //std::ifstream fsR("./Rfuncs.txt");

    //std::ifstream fsR("./ps_python3/plin_z05wmap9_cleftRnew_z000.txt");
    //std::ifstream fsR("./ps_python3/plin_Fr5z1planck_cleftRnew_z000.txt");
    for (int i=0; i<NkTemp; ++i) {
             std::string ssR;
             getline(fsR,ssR);
               //std::istringstream(ssR) >> kaR[i] >> R1[i] >> R2[i];
               std::istringstream(ssR) >> kaR[i] >> R1[i] >> R2[i] >> R12[i] >> RI[i];
               //std::cout << RI[i] << std::endl;
    }
    fsR.close();
   // std::cout << NkTemp << std::endl;
    R1spl.init(kaR,R1);
    R2spl.init(kaR,R2);
    R12spl.init(kaR,R12);
    RIspl.init(kaR,RI);


    //George adding extra part here, to import and spline over time derivatives of Q_n and R_n functions needed for direct RSD
    //Start by declaring these extra functions
    std::vector<double> kaprime(NkTemp),R1prime(NkTemp),R2prime(NkTemp),kaRprime(NkTemp);
    std::vector<double> Q1prime(NkTemp),Q2prime(NkTemp),Q5prime(NkTemp),Q8prime(NkTemp),Q3prime(NkTemp);
    std::vector<double> Qsprime(NkTemp),QIprime(NkTemp),R12prime(NkTemp),RIprime(NkTemp), Q1halfprime(NkTemp) ;   
    //std::ifstream fsQder("./ps_python3/plin_N5z05wmap9_cleftQnewDer_z000.txt");
    std::ifstream fsQder("./plin_Fr6z05wmap9_cleftQnewDer_z000.txt");
    //std::ifstream fsQder("./ps_python3/plin_Fr5z05wmap9_cleftQnewDer_z000.txt");
    //std::ifstream fsQder("./ps_python3/plin_Fr5z1planck_cleftQnewDer_z000.txt");
    //std::ifstream fsQder("./ps_python3/plin_N1z05wmap9_cleftQ0_z000.txt");
    //std::ifstream fsQder("./Q_zeros.txt");
    //std::ifstream fsQder("./Qfuncs.txt");
    //std::ifstream fsQder("./Qdotfuncs.txt");
    for (int i=0; i<NkTemp; ++i) {
      std::string ssQder;
      getline(fsQder,ssQder);
               //std::istringstream(ssQ) >> ka[i] >> Q1[i] >> Q2[i] >> Q3[i] >> Q5[i] >> Q8[i] >> Qs[i];  
               std::istringstream(ssQder) >> kaprime[i] >> Q1prime[i] >> Q2prime[i] >> Q3prime[i] >> Q5prime[i] >> Q8prime[i] >> Qsprime[i] >> QIprime[i] >> Q1halfprime[i];
          // Switching to new format to import from Alejandos code. Old format above.
          //   std::istringstream(ssQder) >> kaprime[i] >> Q1prime[i] >> Q2prime[i] >> Q5prime[i] >> Q8prime[i] >> QIprime[i] >> Q1halfprime[i];
               //std::cout << QI[i] << std::endl;
    }
    fsQder.close();
    //And spline the new functions
    Q1primespl.init(kaprime,Q1prime);
    Q2primespl.init(kaprime,Q2prime);
    Q5primespl.init(kaprime,Q5prime);
    Q8primespl.init(kaprime,Q8prime);
    //Qprimesspl.init(kaprime,Qsprime); //Not calculating this for now.
    Qprimesspl.init(kaprime,Qs); //Just storing something in there, for the sake of storing, temporarily.
    QIprimespl.init(kaprime,QIprime);
    Q1halfprimespl.init(kaprime,Q1halfprime);

    //std::ifstream fsRder("./ps_python3/plin_N5z05wmap9_cleftRnewDer_z000.txt");    
    std::ifstream fsRder("./plin_Fr6z05wmap9_cleftRnewDer_z000.txt");
    //std::ifstream fsRder("./ps_python3/plin_Fr5z05wmap9_cleftRnewDer_z000.txt");
    //std::ifstream fsRder("./ps_python3/plin_Fr5z1planck_cleftRnewDer_z000.txt");
    //std::ifstream fsRder("./ps_python3/plin_N1z05wmap9_cleftR0_z000.txt");
   // std::ifstream fsRder("./R_zeros.txt");
    //std::ifstream fsRder("./Rfuncs.txt");
    //std::ifstream fsRder("./Rdotfuncs.txt"); //Importing from Alejandro's code.
    for (int i=0; i<NkTemp; ++i) {
             std::string ssRder;
             getline(fsRder,ssRder);
               //std::istringstream(ssR) >> kaR[i] >> R1[i] >> R2[i];
               std::istringstream(ssRder) >> kaRprime[i] >> R1prime[i] >> R2prime[i] >> R12prime[i] >> RIprime[i];
    }
    fsRder.close();
    //And spline the new functions
    R1primespl.init(kaRprime,R1prime);
    R2primespl.init(kaRprime,R2prime);
    R12primespl.init(kaRprime,R12prime);
    RIprimespl.init(kaRprime,RIprime);

    //std::ifstream grth("./fgrowth_N5z05wmap9.txt");
    std::ifstream grth("./fgrowth_Fr6z05wmap9.txt");
    //std::ifstream grth("./fgrowth_Fr5z05wmap9.txt");
    //std::ifstream grth("./fgrowth_Fr5z1planck.txt");

    //std::stringstream ss;
    //ss<<"fgrowth_"<< fname;
    
    //std::ifstream grth(ss.str().c_str());

    const int NkTempgrowth=2000;
    std::vector<double> kgrowth(NkTempgrowth), fgrowth(NkTempgrowth);
    for (int i=0; i<NkTempgrowth; ++i) {
             std::string ssgrth;
             getline(grth,ssgrth);
               std::istringstream(ssgrth) >> kgrowth[i] >> fgrowth[i];
    }
    grth.close();
    //And spline the new functions
    fgrowthspl.init(kgrowth,fgrowth);
    //std::cout  << "here" << std::endl;

}


// Returns the different contributions to the real-space correlation function (component 0)
// mean infall velocity (component 1), and the velocity dispersion (component 2)
// for locally biased tracers.
// For mean infall velocity (component 1) only the line-of-sight component is returned
// and the result should be multiplied by f and divided by 1+xi(real).
// For velocity dispersion (component 2) both sigma_perp and sigma_par are returned
// and the result should be multiplied by f^2 and divided by 1+xi(real). NOTE we return
// parallel then perpendicular/transverse.
// This is not tested for very large or small values of r.
// The integration is over x=q-r, in length and angle with the
// azimuthal integral being trivial.
std::vector<std::vector<double>>
LSM::dvsPair(const double rval)
{
    const double xmin=0;
    const double xmax=10*sqrt(sigma2);
    const double rr[3]={0,0,rval};
    const double r2   =rval*rval;
    const int    Nx=500;
    const double dx=(xmax-xmin)/Nx;

    std::vector<double> xi(12);
    std::vector<double> vv(10);
    std::vector<double> ss(16);

    for (int ixx=0; ixx<Nx; ++ixx) {
        double xx=xmin+(ixx+0.5)*dx;
        double x2=xx*xx;
        for (int imu=0; imu<gl.N; ++imu) {
            double mu = gl.x[imu];

            // Compute vec{q}=vec{r}+vec{x} with vec{r}=r.zhat,
            // so r_z=r, x_z=x*mu, cos_rq=(r_z+x_z)/q.
            double qlen = sqrt(r2+x2+2*rval*xx*mu);
            double qcos = (rval+xx*mu)/qlen;
            double qsin = sqrt(1-qcos*qcos);
            double qq[3] = {qlen*qsin,0,qlen*qcos};
            double qh[3] = {     qsin,0,     qcos};
            if (qlen>qmin && qlen<qmax) {
                // We keep the Zeldovich piece exponentiated and expand down
                // the 1-loop piece.
                double pref = x2 * zeldovichIntegrand(rr,qq,0) * gl.w[imu];
                // For the bias terms, compute U, xi and Ainv (even though in above).
                std::vector<double> qf  =interpQfuncs(qlen);
                std::vector<double> ef  =interpEfuncs(qlen);
                std::vector<double> Ainv=calcAinv(qq);
                std::vector<double> Aloop(9);
                std::vector<double> Alin(9);
                std::vector<double> Adot(9);
                std::vector<double> Addot(9);

                double Xdot=ef[0]+2*ef[1]+4*ef[2],Ydot=ef[3]+2*ef[4]+4*ef[5];
                double Xddot=ef[0]+4*ef[1]+6*ef[2],Yddot=ef[3]+4*ef[4]+6*ef[5];
                for (int i=0; i<3; ++i) {
                    for (int j=0; j<3; ++j) {
                        Aloop[3*i+j] = (ef[4]+2*ef[5])*qh[i]*qh[j]+(ef[1]+2*ef[2])*(i==j);
                        Adot[ 3*i+j] = Ydot*qh[i]*qh[j]+Xdot*(i==j);
                        Alin[ 3*i+j] = ef[3]*qh[i]*qh[j]+ef[0]*(i==j);
                        Addot[3*i+j] = Yddot*qh[i]*qh[j]+Xddot*(i==j);
	                //Alindotdot[ 3*i+j] = ef[3]*qh[i]*qh[j]+ef[0]*(i==j);
                    }
                }
                double xiL=qf[3];
                // Construct the auxilliary matrix/vectors g, G of CLPT Eq. (45)
                // and Gamma of Eq. (75).
                double g[3],UL[3],U[3],Udot[3],U20[3],U11[3],G[9],W[27],Wddot[27];
                for (int i=0; i<3; ++i) {
                    g[i]=0;
                    for (int j=0; j<3; ++j)
                        g[i] += Ainv[3*i+j]*(qq[j]-rr[j]);
                    UL[i] = ef[11]*qh[i];
                    U[ i] =(ef[11]+ef[12])*qh[i];
                    Udot[i]  =(ef[11]+3*ef[12])*qh[i];
                    U20[i]=ef[13]*qh[i];
                    U11[i]=ef[14]*qh[i];
                }

                for (int i=0; i<3; ++i)
                    for (int j=0; j<3; ++j)
                        G[3*i+j]=Ainv[3*i+j]-g[i]*g[j];

                double GA=0;
                for (int i=0; i<3; ++i)
                    for (int j=0; j<3; ++j)
                        GA += Aloop[3*i+j]*G[3*i+j];

                for (int i=0; i<3; ++i)
                    for (int j=0; j<3; ++j)
                        for (int k=0; k<3; ++k){

                            W[9*i+3*j+k] = ef[8]*qh[i]*(j==k) + ef[8]*qh[j]*(i==k) + ef[9]*qh[k]*(i==j) +ef[10]*qh[i]*qh[j]*qh[k];

                }
                for (int i=0; i<3; ++i)
                    for (int j=0; j<3; ++j)
                        for (int k=0; k<3; ++k)
                            Wddot[9*i+3*j+k] = 2*W[9*i+3*j+k]+2*W[9*i+3*k+j]+W[9*k+3*j+i];

                // We also need \ddot{A}^{10}:
                double A10[9];
                for (int i=0; i<3; ++i)
                    for (int j=0; j<3; ++j)
                        A10[3*i+j] = (4*ef[7])*qh[i]*qh[j] + (4*ef[6])*(i==j);

                double GW=0;
                double V1=ef[8];
                double V3=ef[9];
                double Tq=ef[10];
                for (int i=0; i<3; ++i) {
                    for (int j=0; j<3; ++j) {
                        for (int k=0; k<3; ++k) {
                            double Gam,W;
                            Gam = Ainv[3*i+j]*g[k]+Ainv[3*k+i]*g[j]+Ainv[3*j+k]*g[i]
                                - g[i]*g[j]*g[k];
                            W   = Tq*qh[i]*qh[j]*qh[k];
                            if (j==k) W += V1*qh[i];
                            if (i==k) W += V1*qh[j];
                            if (i==j) W += V3*qh[k];
                            GW += Gam*W;
                        }
                    }
                }
                GW *= 3;	// Account for permutations.
                double trG = 0, Ug = 0, ULg = 0, U2 = 0, gq = 0, qG = 0, gA = 0, UGA = 0, qGq = 0, gAL = 0, gAU = 0, AGA = 0;
                for (int i=0; i<3; ++i) {
                    gq += g[i]*qh[i];
                    gA += g[i]*Adot[3*2+i];
                    gAL+= g[i]*Alin[3*2+i];
                    Ug += U[i]*g[i];
                    ULg+=UL[i]*g[i];
                    U2 +=UL[i]*UL[i];
                    qG += qh[i]*G[3*2+i];
                    trG+= G[3*i+i];
                    for (int j=0; j<3; ++j) {
                        UGA += UL[i]*G[3*i+j]*Alin[3*2+j];
                        qGq += qh[i]*G[3*i+j]*qh[j];
                        gAU += g[i]*Alin[3*i+j]*UL[j];
                        AGA += Alin[3*2+i]*G[3*i+j]*Alin[3*2+j];
                    }
                }

                double GWv=0;
                for (int i=0; i<3; ++i)
                    for (int j=0; j<3; ++j) {
                        double Wdot_ijn=(3*V1+V3)*(qh[i]*(j==2)+qh[j]*(i==2))+
                                  2*(V1+V3)*qh[2]*(i==j)+
                                  4*Tq*qh[i]*qh[j]*qh[2];
                        GWv += G[3*i+j]*Wdot_ijn;
                    }

                double U20g=0,U11g=0;
                for (int i=0; i<3; ++i) {
                    U20g += U20[i]*g[i];
                    U11g += U11[i]*g[i];
                }
                double UUG=0,qqG=0;
                for (int i=0; i<3; ++i)
                    for (int j=0; j<3; ++j) {
                        UUG += G[3*i+j]*UL[i]*UL[j];
                        qqG += G[3*i+j]*qh[i]*qh[j];
                    }
                double A10G=2*trG*ef[6] + 2*qqG*ef[7];
                double d2xiLin=qf[5];

                double gA10=3*(ef[6]*g[2]+ef[7]*gq*qh[2]);

                // The mode-coupling term, then add the <s^2 Delta Delta> term:
                double shear=ef[15]*gq;
                for (int i=0; i<3; ++i)
                    for (int j=0; j<3; ++j) {
                        double upsilon= qh[i]*qh[j]*(3*qf[6]*qf[6]+4*qf[6]*qf[7]+
                                  2*qf[6]*qf[8]+2*qf[7]*qf[7]+4*qf[7]*qf[8]+
                                  qf[8]*qf[8]) + (i==j)*2*qf[7]*qf[7];
                        shear += G[3*i+j]*upsilon;
                    }
                shear *= 2;

                // The mode-coupling term, then add the <s^2 Delta Delta> term:
                double shear_v=2*ef[15]*g[2];
                for (int i=0; i<3; ++i) {
                    double upsilon= qh[i]*qh[2]*(3*qf[6]*qf[6]+4*qf[6]*qf[7]+
                                2*qf[6]*qf[8]+2*qf[7]*qf[7]+4*qf[7]*qf[8]+
                                qf[8]*qf[8]) + (i==2)*2*qf[7]*qf[7];
                    shear_v -= g[i]*upsilon;
                }
                shear_v *= 2;
                double V12=qf[9]*gq;

                // Now do the 1, Fp, Fpp, Fp^2, Fp.Fpp, Fpp^2, b_nabla^2,
                // bs, b1.bs2, b2.bs2, bs2^2 terms.
                xi[ 0] +=   pref *(1-GA/2.+GW/6.); // +Aeft*trG*eftNorm);
                xi[ 1] +=  -pref *(2*Ug+A10G);
                xi[ 2] +=  -pref *(UUG+U20g);
                xi[ 3] +=   pref *(xiL-UUG-U11g);
                xi[ 4] +=  -pref *(2*xiL*ULg);
                xi[ 5] +=   pref *xiL*xiL/2;
                xi[ 6] +=  -pref *0.5*trG;
                xi[ 7] +=   pref *d2xiLin;
                xi[ 8] +=  -pref *shear;
                xi[ 9] +=  -pref *2*V12;
                xi[10] +=   pref *qf[10];
                xi[11] +=   pref *qf[11];

                // Now do the 1, Fp, Fpp, Fp^2, Fp.Fpp, Fpp^2, grad_xi, g_los,
                // bs2 terms.
                vv[0] +=   -pref *(gA+0.5*GWv);
                vv[1] +=  2*pref *(Udot[2]-UGA-gA10);
                vv[2] +=    pref *(2*U20[2]-2*ULg*UL[2]);
                vv[3] +=   -pref *(ULg*UL[2]+ULg*UL[2]+xiL*gAL-2*U11[2]);
                vv[4] +=  2*pref *xiL*UL[2];
                vv[5] +=  0;
                vv[6] +=    pref *qf[4];
                vv[7] +=    pref *g[2];
                vv[8] +=    pref *shear_v;
                vv[9] +=    pref *V12;

                double Wg=0;
                for (int i=0; i<3; ++i) Wg += Wddot[9*i+3*2+2]*g[i];
                // Now the shear term.
                double upsilon= qh[2]*qh[2]*(3*qf[6]*qf[6]+4*qf[6]*qf[7]+
                          2*qf[6]*qf[8]+2*qf[7]*qf[7]+4*qf[7]*qf[8]+
                          qf[8]*qf[8]) + 2*qf[7]*qf[7];
                double shear_s = 2*upsilon;
                // Now do the 1, Fp, Fpp, Fp^2, Fp.Fpp, Fpp^2 terms for \sigma_par^2.
                ss[ 0] +=    pref *(Addot[3*2+2]-AGA-Wg);
                //ss[ 0] +=    pref *(Alin[3*2+2]-AGA-0*Wg); //Zeldovich
                ss[ 1] += -2*pref *(2*gAL*UL[2]+ULg*Alin[3*2+2]-A10[3*2+2]);
                //ss[ 1] += -2*pref *(2*gAL*UL[2]+ULg*Alin[3*2+2]-0*A10[3*2+2]); //Zeldovich
                ss[ 2] +=  2*pref *   UL[2]*UL[2];
                ss[ 3] +=    pref *(2*UL[2]*UL[2]+xiL*Alin[3*2+2]);
                ss[ 4]  =  0;
                ss[ 5]  =  0;
                ss[ 6]  =    pref *2*shear_s;
                ss[ 7]  =    pref *1*qf[3];
                // Next work out the trace components, i.e. summed over n=m.
                Wg=0;
                for (int i=0; i<3; ++i)
                    for (int j=0; j<3; ++j)
                        Wg += Wddot[9*i+3*j+j]*g[i];
                AGA=0;
                for (int m=0; m<3; ++m) {
                    for (int i=0; i<3; ++i) {
                        for (int j=0; j<3; ++j) {
                            AGA += Alin[3*m+i]*G[3*i+j]*Alin[3*j+m];
                        }
                    }
                }
                double trA=0,trAL=0,trA10=0;
                for (int m=0; m<3; ++m) {
                    trA  += Addot[3*m+m];
                    trAL +=  Alin[3*m+m];
                    trA10+=   A10[3*m+m];
                }
                // Now the shear term.
                upsilon=0;
                for (int m=0; m<3; ++m)
                    upsilon += qh[m]*qh[m]*(3*qf[6]*qf[6]+4*qf[6]*qf[7]+
                           2*qf[6]*qf[8]+2*qf[7]*qf[7]+4*qf[7]*qf[8]+
                           qf[8]*qf[8]) + 2*qf[7]*qf[7];
                shear_s = 2*upsilon;
                ss[ 8] +=    pref *(trA-AGA-Wg);
                //ss[ 8] +=    pref *(trAL-AGA-0*Wg); //Zeldovich
                ss[ 9] += -2*pref *(2*gAU+ULg*trAL-trA10);
                //ss[ 9] += -2*pref *(2*gAU+ULg*trAL-0*trA10); //Zeldovich
                ss[10] +=  2*pref * U2;
                ss[11] +=    pref *(2*U2+xiL*trAL);
                ss[12]  =  0;
                ss[13]  =  0;
                ss[14]  =    pref *2*shear_s;
                ss[15]  =    pref *3*qf[3];
            }
        }
    }
    for (int j=0; j<xi.size(); ++j) {
        xi[j] *= dx;      // Convert sum to integral.
        xi[j] *= 2*M_PI;  // The azimuthal integral.
    }
    xi[0] -= 1.0;	// Calculated 1+xi, subtract 1 for xi.

    for (int j=0; j<vv.size(); ++j) {
        vv[j] *= dx;      // Convert sum to integral.
        vv[j] *= 2*M_PI;  // The azimuthal integral.
    }

    // Now sigma_perp is related to the trace by s_perp^2=(1/2)[Tr-sig_par^2]
    for (int j=8; j<ss.size(); ++j)
        ss[j] = 0.5*(ss[j] - ss[j-8]);
    for (int j=0; j<ss.size(); ++j) {
        ss[j] *= dx;	// Convert sum to integral.
        ss[j] *= 2*M_PI;	// The azimuthal integral.
    }
    std::vector<std::vector<double>> res = {xi, vv, ss};
    return res;
}


// MG version of subroutine dvsPair, created by George.
// Returns the different contributions to the real-space correlation function (component 0)
// mean infall velocity (component 1), and the velocity dispersion (component 2)
// for locally biased tracers in MG
// Due to the scale dependence in MG theories, the result should NOT be multiplied by f^2, but should be divided by 1+xi(real). 
// NOTE we return
// parallel then perpendicular/transverse.
// This is not tested for very large or small values of r.
// The integration is over x=q-r, in length and angle with the
// azimuthal integral being trivial.
std::vector<std::vector<double>>
LSM::dvsPairMG(const double rval)
{
    const double xmin=0;
    const double xmax=10*sqrt(sigma2);
    //const double xmax=10*sqrt(sigma2);
    const double rr[3]={0,0,rval};
    const double r2   =rval*rval;
    const int    Nx=500;
    const double dx=(xmax-xmin)/Nx;
   //std::cout  << qmax << std::endl;
    std::vector<double> xi(12);
    std::vector<double> vv(10);
    std::vector<double> ss(16);

    for (int ixx=0; ixx<Nx; ++ixx) {
        double xx=xmin+(ixx+0.5)*dx;
        double x2=xx*xx;
        for (int imu=0; imu<gl.N; ++imu) {
            double mu = gl.x[imu];

            // Compute vec{q}=vec{r}+vec{x} with vec{r}=r.zhat,
            // so r_z=r, x_z=x*mu, cos_rq=(r_z+x_z)/q.
            double qlen = sqrt(r2+x2+2*rval*xx*mu);
           
           //std::cout<< std::scientific << qlen <<std::endl;
            double qcos = (rval+xx*mu)/qlen;
            double qsin = sqrt(1-qcos*qcos);
            double qq[3] = {qlen*qsin,0,qlen*qcos};
            double qh[3] = {     qsin,0,     qcos};
            if (qlen>qmin && qlen<qmax) {
                // We keep the Zeldovich piece exponentiated and expand down
                // the 1-loop piece.
                //    std::cout  << "here" << std::endl;
              // if (qlen>500) {
              // std::cout<< std::scientific << qlen <<std::endl;      }
               
                double pref = x2 * zeldovichIntegrand(rr,qq,0) * gl.w[imu];
                //if (rval>300){
                //std::cout<< std::setw(15)<< pref << std::setw(15)<< qlen <<std::endl; 
                // }
                // For the bias terms, compute U, xi and Ainv (even though in above).
                std::vector<double> qf  =interpQfuncs(qlen);
                std::vector<double> ef  =interpEfuncs(qlen);
                std::vector<double> efMG  =interpEfuncsMG(qlen);
                std::vector<double> Ainv=calcAinv(qq);
                std::vector<double> Aloop(9);
                std::vector<double> Alin(9),Alindot(9),Alindotdot(9);
                std::vector<double> Adot(9);
                std::vector<double> Addot(9);

               // for (int i=0; i<3; ++i) {
               //     for (int j=0; j<3; ++j) {
               //         Aloop22[3*i+j] = (ef[4])*qh[i]*qh[j]+(ef[1])*(i==j);
               //         Aloop13[3*i+j] = (2*ef[5])*qh[i]*qh[j]+(2*ef[2])*(i==j);                        
               //         Alin[ 3*i+j] = ef[3]*qh[i]*qh[j]+ef[0]*(i==j);
               //         Aloop22p[3*i+j] = (efMG[4])*qh[i]*qh[j]+(efMG[1])*(i==j);
               //         Aloop13p[3*i+j] = (2*efMG[5])*qh[i]*qh[j]+(2*efMG[2])*(i==j);
               //         Aloop1p3[3*i+j] = (2*efMG[19])*qh[i]*qh[j]+(2*efMG[17])*(i==j);
               //         Aloop2p2p[3*i+j] = (efMG[18])*qh[i]*qh[j]+(efMG[16])*(i==j);
               //         Aloop1p3p[3*i+j] = (2*efMG[21])*qh[i]*qh[j]+(2*efMG[20])*(i==j);
               //     }
               // }
                //    std::cout  << "here" << std::endl;
                //double Xdot=ef[0]+2*ef[1]+4*ef[2],Ydot=ef[3]+2*ef[4]+4*ef[5];
                double Xdot=efMG[38]+efMG[1]+efMG[2]+efMG[17],Ydot=efMG[39]+efMG[4]+efMG[5]+efMG[19];
                //double Xddot=ef[0]+4*ef[1]+6*ef[2],Yddot=ef[3]+4*ef[4]+6*ef[5];
                double Xddot=efMG[41]+efMG[16]+2*efMG[20],Yddot=efMG[42]+efMG[18]+2*efMG[21];
                for (int i=0; i<3; ++i) {
                    for (int j=0; j<3; ++j) {
                        Aloop[3*i+j] = (ef[4]+2*ef[5])*qh[i]*qh[j]+(ef[1]+2*ef[2])*(i==j);
                        Adot[ 3*i+j] = Ydot*qh[i]*qh[j]+Xdot*(i==j);
                        Alin[ 3*i+j] = ef[3]*qh[i]*qh[j]+ef[0]*(i==j);
                        Alindot[ 3*i+j] = efMG[39]*qh[i]*qh[j]+efMG[38]*(i==j);
                        Alindotdot[ 3*i+j] = efMG[42]*qh[i]*qh[j]+efMG[41]*(i==j);
	                //Alindotdot[ 3*i+j] = ef[3]*qh[i]*qh[j]+ef[0]*(i==j);
                        Addot[3*i+j] = Yddot*qh[i]*qh[j]+Xddot*(i==j);
                    }
                }
                //    std::cout  << "here" << std::endl;
                double xiL=qf[3];
                // Construct the auxilliary matrix/vectors g, G of CLPT Eq. (45)
                // and Gamma of Eq. (75).
                double g[3],UL[3],U[3],Udot[3],U20[3],U11[3],G[9],W[27],Wddot[27],Wdot[27],U20dot[3],U11dot[3], ULdot[3];
                for (int i=0; i<3; ++i) {
                    g[i]=0;
                    for (int j=0; j<3; ++j)
                        g[i] += Ainv[3*i+j]*(qq[j]-rr[j]);
                    UL[i] = ef[11]*qh[i];
                    //ULdot[i]  =(efMG[11])*qh[i];
   		    ULdot[i]  =(efMG[40])*qh[i];
                    U[ i] =(ef[11]+ef[12])*qh[i];
                    //Udot[i]  =(ef[11]+3*ef[12])*qh[i];
                    Udot[i]  =(efMG[40]+efMG[12])*qh[i];
                    U20[i]=ef[13]*qh[i];
                    U20dot[i]=efMG[13]*qh[i];
                    U11[i]=ef[14]*qh[i];
                    U11dot[i]=efMG[14]*qh[i];
                }
                //std::cout  << "here" << std::endl;
                for (int i=0; i<3; ++i)
                    for (int j=0; j<3; ++j)
                        G[3*i+j]=Ainv[3*i+j]-g[i]*g[j];

                double GA=0;
                for (int i=0; i<3; ++i)
                    for (int j=0; j<3; ++j)
                        GA += Aloop[3*i+j]*G[3*i+j];

                for (int i=0; i<3; ++i)
                    for (int j=0; j<3; ++j)
                        for (int k=0; k<3; ++k){

                            W[9*i+3*j+k] = ef[8]*qh[i]*(j==k) + ef[8]*qh[j]*(i==k) + ef[9]*qh[k]*(i==j) +ef[10]*qh[i]*qh[j]*qh[k];

                }
                double W112p[27], W11p2[27],W1p1p2[27], W11p2p[27], W1p1p2p[27];
                for (int i=0; i<3; ++i)
                    for (int j=0; j<3; ++j)
                        for (int k=0; k<3; ++k){

                            //W[9*i+3*j+k] = ef[8]*qh[i]*(j==k) + ef[8]*qh[j]*(i==k) + ef[9]*qh[k]*(i==j) +ef[10]*qh[i]*qh[j]*qh[k];
                            W112p[9*i+3*j+k] = efMG[8]*qh[i]*(j==k) + efMG[8]*qh[j]*(i==k) + efMG[9]*qh[k]*(i==j) +efMG[10]*qh[i]*qh[j]*qh[k];
                            W11p2[9*i+3*j+k] = efMG[26]*qh[i]*(j==k) + efMG[26]*qh[j]*(i==k) + efMG[27]*qh[k]*(i==j) +efMG[28]*qh[i]*qh[j]*qh[k];
                            W1p1p2[9*i+3*j+k] = efMG[29]*qh[i]*(j==k) + efMG[29]*qh[j]*(i==k) + efMG[30]*qh[k]*(i==j) +efMG[31]*qh[i]*qh[j]*qh[k];
                            W11p2p[9*i+3*j+k] = efMG[32]*qh[i]*(j==k) + efMG[32]*qh[j]*(i==k) + efMG[33]*qh[k]*(i==j) +efMG[34]*qh[i]*qh[j]*qh[k];
                            W1p1p2p[9*i+3*j+k] = efMG[35]*qh[i]*(j==k) + efMG[35]*qh[j]*(i==k) + efMG[36]*qh[k]*(i==j) +efMG[37]*qh[i]*qh[j]*qh[k];
                }
                for (int i=0; i<3; ++i)
                    for (int j=0; j<3; ++j)
                        for (int k=0; k<3; ++k)
                            //Wddot[9*i+3*j+k] = 0.733327*0.733327*(2*W[9*i+3*j+k]+2*W[9*i+3*k+j]+W[9*k+3*j+i]);
                            Wddot[9*i+3*j+k] = W11p2p[9*i+3*j+k]+W11p2p[9*i+3*k+j]+W1p1p2[9*k+3*j+i];
                            //Wddot[9*i+3*j+k] = W11p2p[9*i+3*j+k]+W11p2p[9*i+3*j+k]+W1p1p2[9*i+3*j+k];

                // We also need \ddot{A}^{10}:
                double A10[9];
                for (int i=0; i<3; ++i)
                    for (int j=0; j<3; ++j)
                        //A10[3*i+j] = 0.733327*0.733327*((4*ef[7])*qh[i]*qh[j] + (4*ef[6])*(i==j));
                        A10[3*i+j] = (2*efMG[25])*qh[i]*qh[j] + (2*efMG[24])*(i==j);

                double GW=0;
                double V1=ef[8];
                double V3=ef[9];
                double Tq=ef[10];
                for (int i=0; i<3; ++i) {
                    for (int j=0; j<3; ++j) {
                        for (int k=0; k<3; ++k) {
                            double Gam,W;
                            Gam = Ainv[3*i+j]*g[k]+Ainv[3*k+i]*g[j]+Ainv[3*j+k]*g[i]
                                - g[i]*g[j]*g[k];
                            W   = Tq*qh[i]*qh[j]*qh[k];
                            if (j==k) W += V1*qh[i];
                            if (i==k) W += V1*qh[j];
                            if (i==j) W += V3*qh[k];
                            GW += Gam*W;
                        }
                    }
                }
                GW *= 3;	// Account for permutations.
                double trG = 0, Ug = 0, ULg = 0, U2 = 0, gq = 0, qG = 0, gA = 0, UGA = 0, qGq = 0, gAL = 0, gAU = 0, AGA = 0;
                for (int i=0; i<3; ++i) {
                    gq += g[i]*qh[i];
                    gA += g[i]*Adot[3*2+i];
                    //gAL+= g[i]*Alin[3*2+i];
                    gAL+= g[i]*Alindot[3*2+i];
                    Ug += U[i]*g[i];
                    ULg+=UL[i]*g[i];
                    U2 +=ULdot[i]*ULdot[i];
                    qG += qh[i]*G[3*2+i];
                    trG+= G[3*i+i];
                    for (int j=0; j<3; ++j) {
                        //UGA += UL[i]*G[3*i+j]*Alin[3*2+j];
                        UGA += UL[i]*G[3*i+j]*Alindot[3*2+j];
                        qGq += qh[i]*G[3*i+j]*qh[j];
                        gAU += g[i]*Alindot[3*i+j]*ULdot[j];
                        //AGA += Alin[3*2+i]*G[3*i+j]*Alin[3*2+j];
                        AGA += Alindot[3*2+i]*G[3*i+j]*Alindot[3*2+j];
                    }
                }


                 for (int i=0; i<3; ++i)
                    for (int j=0; j<3; ++j)
                        for (int k=0; k<3; ++k)
                            //Wddot[9*i+3*j+k] = 2*W[9*i+3*j+k]+2*W[9*i+3*k+j]+W[9*k+3*j+i];
                            Wdot[9*i+3*j+k] = W112p[9*i+3*j+k]+W11p2[9*i+3*k+j]+W11p2[9*k+3*j+i];
                            //Wdot[9*i+3*j+k] = W112p[9*i+3*j+k]+W11p2[9*i+3*j+k]+W11p2[9*i+3*j+k];
                double GWv=0;
                for (int i=0; i<3; ++i)
                    for (int j=0; j<3; ++j) {
                        //double Wdot_ijn=(3*V1+V3)*(qh[i]*(j==2)+qh[j]*(i==2))+
                        //          2*(V1+V3)*qh[2]*(i==j)+
                        //          4*Tq*qh[i]*qh[j]*qh[2];
                        double Wdot_ijn;
                        Wdot_ijn = Wdot[9*2+3*i+j];
                        //Wdot_ijn = Wdot[9*i+3*j+2];
                        GWv += G[3*i+j]*Wdot_ijn;
                    }

                double U20g=0,U11g=0;
                for (int i=0; i<3; ++i) {
                    U20g += U20[i]*g[i];
                    U11g += U11[i]*g[i];
                }
                double UUG=0,qqG=0;
                for (int i=0; i<3; ++i)
                    for (int j=0; j<3; ++j) {
                        UUG += G[3*i+j]*UL[i]*UL[j];
                        qqG += G[3*i+j]*qh[i]*qh[j];
                    }
                double A10G=2*trG*ef[6] + 2*qqG*ef[7];
                double d2xiLin=qf[5];

                //double gA10=3*(ef[6]*g[2]+ef[7]*gq*qh[2]);
                double gA10=((efMG[6]+efMG[22])*g[2]+(efMG[7]+efMG[23])*gq*qh[2]);

                // The mode-coupling term, then add the <s^2 Delta Delta> term:
                double shear=ef[15]*gq;
                for (int i=0; i<3; ++i)
                    for (int j=0; j<3; ++j) {
                        double upsilon= qh[i]*qh[j]*(3*qf[6]*qf[6]+4*qf[6]*qf[7]+
                                  2*qf[6]*qf[8]+2*qf[7]*qf[7]+4*qf[7]*qf[8]+
                                  qf[8]*qf[8]) + (i==j)*2*qf[7]*qf[7];
                        shear += G[3*i+j]*upsilon;
                    }
                shear *= 2;

                // The mode-coupling term, then add the <s^2 Delta Delta> term:
                double shear_v=2*ef[15]*g[2];
                for (int i=0; i<3; ++i) {
                    double upsilon= qh[i]*qh[2]*(3*qf[6]*qf[6]+4*qf[6]*qf[7]+
                                2*qf[6]*qf[8]+2*qf[7]*qf[7]+4*qf[7]*qf[8]+
                                qf[8]*qf[8]) + (i==2)*2*qf[7]*qf[7];
                    shear_v -= g[i]*upsilon;
                }
                shear_v *= 2;
                double V12=qf[9]*gq;

                // Now do the 1, Fp, Fpp, Fp^2, Fp.Fpp, Fpp^2, b_nabla^2,
                // bs, b1.bs2, b2.bs2, bs2^2 terms.
                xi[ 0] +=   pref *(1-GA/2.+GW/6.); // +Aeft*trG*eftNorm);
                //if (rval>300 && xi[0]<-0.01){
                //std::cout<< std::setw(15)<< pref << std::setw(15)<< xi[ 0] << std::setw(15) << rval <<std::endl; 
                // }
                xi[ 1] +=  -0 *(2*Ug+A10G);
                xi[ 2] +=  -pref *(UUG+U20g);
                xi[ 3] +=   pref *(xiL-UUG-U11g);
                xi[ 4] +=  -pref *(2*xiL*ULg);
                xi[ 5] +=   pref *xiL*xiL/2;

                xi[ 6] +=  -pref *0.5*trG;
                xi[ 7] +=   pref *d2xiLin;
                xi[ 8] +=  -pref *shear;
                xi[ 9] +=  -pref *2*V12;
                xi[10] +=   pref *qf[10];
                xi[11] +=   pref *qf[11];

                // Now do the 1, Fp, Fpp, Fp^2, Fp.Fpp, Fpp^2, grad_xi, g_los,
                // bs2 terms.
                vv[0] +=   -pref *(gA+0.5*GWv);
                vv[1] +=  2*pref *(Udot[2]-UGA-gA10);
                //vv[2] +=    pref *(2*U20[2]-2*ULg*UL[2]);
                vv[2] +=    pref *(U20dot[2]-2*ULg*ULdot[2]);
                //vv[3] +=   -pref *(ULg*UL[2]+ULg*UL[2]+xiL*gAL-2*U11[2]);
                vv[3] +=   -pref *(2*ULg*ULdot[2]+xiL*gAL-U11dot[2]);
                vv[4] +=  2*pref *xiL*ULdot[2];
                vv[5] +=  0;

                vv[6] +=    pref *qf[4];
                vv[7] +=    pref *g[2];
                vv[8] +=    pref *shear_v;
                vv[9] +=    pref *V12;

                double Wg=0;
                for (int i=0; i<3; ++i) Wg += Wddot[9*i+3*2+2]*g[i];
                // Now the shear term.
                double upsilon= qh[2]*qh[2]*(3*qf[6]*qf[6]+4*qf[6]*qf[7]+
                          2*qf[6]*qf[8]+2*qf[7]*qf[7]+4*qf[7]*qf[8]+
                          qf[8]*qf[8]) + 2*qf[7]*qf[7];
                double shear_s = 2*upsilon;
                // Now do the 1, Fp, Fpp, Fp^2, Fp.Fpp, Fpp^2 terms for \sigma_par^2.
                ss[ 0] +=    pref *(Addot[3*2+2]-AGA-Wg);
                //ss[ 0] +=    pref *(Alindotdot[3*2+2]-AGA-0*Wg); //Zeldovich

                ss[ 1] += -2*pref *(2*gAL*ULdot[2]+ULg*Alindotdot[3*2+2]-A10[3*2+2]);
                //ss[ 1] += -2*pref *(2*gAL*ULdot[2]+ULg*Alindotdot[3*2+2]-0*A10[3*2+2]); //Zeldovich
                ss[ 2] +=  2*pref *   ULdot[2]*ULdot[2];
                ss[ 3] +=    pref *(2*ULdot[2]*ULdot[2]+xiL*Alindotdot[3*2+2]);
                ss[ 4]  =  0;

                ss[ 5]  =  0;
                ss[ 6]  =    pref *2*shear_s;
                ss[ 7]  =    pref *1*qf[3];
                // Next work out the trace components, i.e. summed over n=m.
                Wg=0;
                for (int i=0; i<3; ++i)
                    for (int j=0; j<3; ++j)
                        Wg += Wddot[9*i+3*j+j]*g[i];
                AGA=0;
                for (int m=0; m<3; ++m) {
                    for (int i=0; i<3; ++i) {
                        for (int j=0; j<3; ++j) {
                            AGA += Alindot[3*m+i]*G[3*i+j]*Alindot[3*j+m];
                        }
                    }
                }
                double trA=0,trAL=0,trALdotdot=0,trA10=0;
                for (int m=0; m<3; ++m) {
                    trA  += Addot[3*m+m];
                    trAL +=  Alin[3*m+m];
                    trALdotdot +=  Alindotdot[3*m+m];
                    trA10+=   A10[3*m+m];
                }
                // Now the shear term.
                upsilon=0;
                for (int m=0; m<3; ++m)
                    upsilon += qh[m]*qh[m]*(3*qf[6]*qf[6]+4*qf[6]*qf[7]+
                           2*qf[6]*qf[8]+2*qf[7]*qf[7]+4*qf[7]*qf[8]+
                           qf[8]*qf[8]) + 2*qf[7]*qf[7];
                shear_s = 2*upsilon;

                ss[ 8] +=    pref *(trA-AGA-Wg);
                //ss[ 8] +=    pref *(trALdotdot-AGA-0*Wg); //Zeldovich
                ss[ 9] += -2*pref *(2*gAU+ULg*trALdotdot-trA10);
                //ss[ 9] += -2*pref *(2*gAU+ULg*trALdotdot-0*trA10); //Zeldovich
                ss[10] +=  2*pref * U2;
                ss[11] +=    pref *(2*U2+xiL*trALdotdot);
                ss[12]  =  0;
                ss[13]  =  0;

                ss[14]  =    pref *2*shear_s;
                ss[15]  =    pref *3*qf[3];
            }
        }
    }
    //std::cout  << "here" << std::endl;
    for (int j=0; j<xi.size(); ++j) {
        xi[j] *= dx;      // Convert sum to integral.
        xi[j] *= 2*M_PI;  // The azimuthal integral.
    }
    xi[0] -= 1.0;	// Calculated 1+xi, subtract 1 for xi.

    for (int j=0; j<vv.size(); ++j) {
        vv[j] *= dx;      // Convert sum to integral.
        vv[j] *= 2*M_PI;  // The azimuthal integral.
    }

    // Now sigma_perp is related to the trace by s_perp^2=(1/2)[Tr-sig_par^2]
    for (int j=8; j<ss.size(); ++j)
        ss[j] = 0.5*(ss[j] - ss[j-8]);
    for (int j=0; j<ss.size(); ++j) {
        ss[j] *= dx;	// Convert sum to integral.
        ss[j] *= 2*M_PI;	// The azimuthal integral.
    }
    std::vector<std::vector<double>> res = {xi, vv, ss};
    //std::cout  << "here" << std::endl;
    return res;
}

//Added by George, calculate contributions to RSD 1-loop CLPT \xi, using "Direct" Lagrangian approach
std::vector<double> 
LSM::xiloopContributionsRSD(const double rval, const double mu, const double f1) {
    // Returns the different contributions to the redshift-space Zel'dovich
    // correlation function for locally biased tracers.
    // This is not tested for very large or small values of r.
    // The integration is over x=q-r, in length and angle with the
    // azimuthal integral being done explicitly.
    const double xmin=0;
    const double xmax=10*sqrt(sigma2);
    const double rr[3]={rval*sqrt(1-mu*mu),0,rval*mu};
    const double r2   =rval*rval;
    const int    Nx=256,Nphi=32;
    const double dx=(xmax-xmin)/Nx;
    const double dphi=2*M_PI/Nphi;
    std::vector<double> xi(6);
    for (int ixx=0; ixx<Nx; ++ixx) {
      double xx=xmin+(ixx+0.5)*dx;
      double x2=xx*xx;
      for (int imu=0; imu<gl.N; ++imu) {
        double mu = gl.x[imu];
        double st = sqrt(1-mu*mu);
        for (int iphi=0; iphi<Nphi; ++iphi) {
          double phi = (iphi+0.5)*dphi;
          double qq[3],qh[3],xv[3];
          xv[0]=xx*st*cos(phi);xv[1]=xx*st*sin(phi);xv[2]=xx*mu;
          double qlen=0;
          for (int i=0; i<3; ++i) {
            qq[i] = xv[i]+rr[i];
            qlen += qq[i]*qq[i];
          }
          qlen = sqrt(qlen);
          for (int i=0; i<3; ++i) qh[i] = qq[i]/qlen;
          if (qlen>qmin && qlen<qmax) {
            // For the unbiased tracers we only need this--all other terms
            // are multiplied by this anyway.
            double pref = x2 * zeldovichIntegrand(rr,qq,f1) * gl.w[imu];
                           // For the bias terms, compute U, xi and Ainv (even though in above).
                std::vector<double> qf  =interpQfuncs(qlen);
                std::vector<double> ef  =interpEfuncs(qlen);
                std::vector<double> Ainv=calcAinv(qq);
                std::vector<double> Aloop(9),Aloop22(9),Aloop13(9);
                std::vector<double> Alin(9);

                for (int i=0; i<3; ++i) {
                    for (int j=0; j<3; ++j) {
                        Aloop22[3*i+j] = (ef[4])*qh[i]*qh[j]+(ef[1])*(i==j);
                        Aloop13[3*i+j] = (2*ef[5])*qh[i]*qh[j]+(2*ef[2])*(i==j);                        
                        Alin[ 3*i+j] = ef[3]*qh[i]*qh[j]+ef[0]*(i==j);
                    }
                }
                //Adding RSD contributions to Ainv and Alin
                for (int i=0; i<3; ++i) {
                   Ainv[3*i+2] /= (1+f1);
                   Ainv[3*2+i] /= (1+f1);
                   Alin[3*i+2] *= (1+f1);
                   Alin[3*2+i] *= (1+f1);
                   Aloop22[3*i+2] *= (1+2*f1);
                   Aloop22[3*2+i] *= (1+2*f1);
                   Aloop13[3*i+2] *= (1+f1);
                   Aloop13[3*2+i] *= (1+3*f1);
                  }
                 Ainv[9] /= (1+f1)*(1+f1); //RSD change to determinant of Ainv

                //Add two contributions to Aloop, after we shifted each of them to RSD
                for (int i=0; i<3; ++i) {
                    for (int j=0; j<3; ++j) {
                         Aloop[3*i+j] = Aloop22[3*i+j]+ Aloop13[3*i+j] ;                        
                    }
                }


                double xiL=qf[3];
                // Construct the auxilliary matrix/vectors g, G of CLPT Eq. (45)
                // and Gamma of Eq. (75).
                double g[3],UL[3],U[3],Udot[3],U20[3],U11[3],G[9],W[27];
                for (int i=0; i<3; ++i) {
                    g[i]=0;
                    for (int j=0; j<3; ++j)
                        g[i] += Ainv[3*i+j]*(qq[j]-rr[j]);
                    UL[i] = ef[11]*qh[i];
                    U[ i] =(ef[11]+ef[12])*qh[i];
                    U20[i]=ef[13]*qh[i];
                    U11[i]=ef[14]*qh[i];
                }

                UL[2] *= 1+f1;	// Correct U as U->RU.
                U[2] = UL[2] + (1+3*f1)*(ef[12])*qh[2];
                U20[2] *= 1+2*f1;
                U11[2] *= 1+2*f1;

                for (int i=0; i<3; ++i)
                    for (int j=0; j<3; ++j)
                        G[3*i+j]=Ainv[3*i+j]-g[i]*g[j];

                double GA=0;
                for (int i=0; i<3; ++i)
                    for (int j=0; j<3; ++j)
                        GA += Aloop[3*i+j]*G[3*i+j];

                for (int i=0; i<3; ++i)
                    for (int j=0; j<3; ++j)
                        for (int k=0; k<3; ++k)
                            W[9*i+3*j+k] = ef[8]*qh[i]*(j==k) + ef[8]*qh[j]*(i==k)
                                + ef[9]*qh[k]*(i==j) +ef[10]*qh[i]*qh[j]*qh[k];
                //Add RSD contribution to Wijk
                for (int i=0; i<3; ++i){
                    for (int j=0; j<3; ++j){
                            W[9*i+3*j+2] *=(1+f1) ;
                            W[9*i+3*2+j] *=(1+f1) ;
                            W[9*2+3*i+j] *=(1+2*f1) ;
                  }
                }
                // We also need \ddot{A}^{10}:
                double A10[9];
                for (int i=0; i<3; ++i)
                    for (int j=0; j<3; ++j)
                        A10[3*i+j] = (4*ef[7])*qh[i]*qh[j] + (4*ef[6])*(i==j);

                //Adding RSD contributions to A10
                for (int i=0; i<3; ++i) {
                   A10[3*i+2] *= (1+f1);
                   A10[3*2+i] *= (1+2*f1);
                  }

                double GW=0;
                double V1=ef[8];
                double V3=ef[9];
                double Tq=ef[10];
                for (int i=0; i<3; ++i) {
                    for (int j=0; j<3; ++j) {
                        for (int k=0; k<3; ++k) {
                            double Gam,Wel;
                            Gam = Ainv[3*i+j]*g[k]+Ainv[3*k+i]*g[j]+Ainv[3*j+k]*g[i]
                                - g[i]*g[j]*g[k];
                            //Wel   = Tq*qh[i]*qh[j]*qh[k];
                            //if (j==k) Wel += V1*qh[i];
                            //if (i==k) Wel += V1*qh[j];
                            //if (i==j) Wel += V3*qh[k];
                            Wel = W[9*i+3*j+k];
                            GW += Gam*Wel;
                        }
                    }
                }
                GW *= 3;	// Account for permutations.
                double trG = 0, Ug = 0, ULg = 0, U2 = 0, gq = 0, qG = 0, gA = 0, UGA = 0, qGq = 0, gAL = 0, gAU = 0, AGA = 0;
                for (int i=0; i<3; ++i) {
                    gq += g[i]*qh[i];
                    gAL+= g[i]*Alin[3*2+i];
                    Ug += U[i]*g[i];
                    ULg+=UL[i]*g[i];
                    U2 +=UL[i]*UL[i];
                    qG += qh[i]*G[3*2+i];
                    trG+= G[3*i+i];
                    for (int j=0; j<3; ++j) {
                        UGA += UL[i]*G[3*i+j]*Alin[3*2+j];
                        qGq += qh[i]*G[3*i+j]*qh[j];
                        gAU += g[i]*Alin[3*i+j]*UL[j];
                        AGA += Alin[3*2+i]*G[3*i+j]*Alin[3*2+j];
                    }
                }

                double U20g=0,U11g=0;
                for (int i=0; i<3; ++i) {
                    U20g += U20[i]*g[i];
                    U11g += U11[i]*g[i];
                }
                double UUG=0,qqG=0;
                for (int i=0; i<3; ++i)
                    for (int j=0; j<3; ++j) {
                        UUG += G[3*i+j]*UL[i]*UL[j];
                        qqG += G[3*i+j]*qh[i]*qh[j];
                    }
                double A10G=2*trG*ef[6] + 2*qqG*ef[7];





                // Now do the 1, Fp, Fpp, Fp^2, Fp.Fpp, Fpp^2
                xi[ 0] +=   pref *(1-GA/2.+GW/6.); // +Aeft*trG*eftNorm);
                xi[ 1] +=  -pref *(2*Ug+A10G);
                xi[ 2] +=  -pref *(UUG+U20g);
                xi[ 3] +=   pref *(xiL-UUG-U11g);
                xi[ 4] +=  -pref *(2*xiL*ULg);
                xi[ 5] +=   pref *xiL*xiL/2;


          }
        }
      }
    }
    for (int j=0; j<xi.size(); ++j) {
      xi[j] *= dx*dphi;	// Convert sum to integral.
    }
    xi[0] -= 1.0;	// Calculated 1+xi, subtract 1 for xi.
    return(xi);
  }


std::vector<double> 
LSM::xiloopContributionsRSDmult(const double rval, const double f1) {
    // Returns the contributions to the multipoles of the redshift-space
    // correlation function for locally biased tracers.
    // This is not tested for very large or small values of r.
    const int Nmu=4;
    GaussLegendre gg = GaussLegendre(2*Nmu);	// Must be even.
    // For even lengths, can sum over half of the points.
    std::vector<double> xiell;
    try{xiell.resize(12);}catch(std::exception& e) {myexception(e);}
    for (int i=0; i<Nmu; ++i) {
      std::vector<double> ximu = xiloopContributionsRSD(rval,gg.x[i],f1);
      double p0=1.0;
      double p2=0.5*(3*gg.x[i]*gg.x[i]-1);
      for (int j=0; j<ximu.size(); ++j) {
        xiell[0*ximu.size()+j] += ximu[j]*gg.w[i] * p0 * 1;
        xiell[1*ximu.size()+j] += ximu[j]*gg.w[i] * p2 * 5;
      }
    }
    return(xiell);
  }


//Added by George, calculate contributions to RSD 1-loop CLPT \xi in MG, using "Direct" Lagrangian approach
std::vector<double> 
LSM::xiloopContributionsRSDMG(const double rval, const double mu, const double f1) {
    // Returns the different contributions to the redshift-space Zel'dovich
    // correlation function for locally biased tracers.
    // This is not tested for very large or small values of r.
    // The integration is over x=q-r, in length and angle with the
    // azimuthal integral being done explicitly.
    const double xmin=0;
    const double xmax=10*sqrt(sigma2);
    const double rr[3]={rval*sqrt(1-mu*mu),0,rval*mu};
    const double r2   =rval*rval;
    const int    Nx=256,Nphi=32;
    const double dx=(xmax-xmin)/Nx;
    const double dphi=2*M_PI/Nphi;
    std::vector<double> xi(6);
    for (int ixx=0; ixx<Nx; ++ixx) {
      double xx=xmin+(ixx+0.5)*dx;
      double x2=xx*xx;
      for (int imu=0; imu<gl.N; ++imu) {
        double mu = gl.x[imu];
        double st = sqrt(1-mu*mu);
        for (int iphi=0; iphi<Nphi; ++iphi) {
          double phi = (iphi+0.5)*dphi;
          double qq[3],qh[3],xv[3];
          xv[0]=xx*st*cos(phi);xv[1]=xx*st*sin(phi);xv[2]=xx*mu;
          double qlen=0;
          for (int i=0; i<3; ++i) {
            qq[i] = xv[i]+rr[i];
            qlen += qq[i]*qq[i];
          }
          qlen = sqrt(qlen);
          for (int i=0; i<3; ++i) qh[i] = qq[i]/qlen;
          if (qlen>qmin && qlen<qmax) {
            // For the unbiased tracers we only need this--all other terms
            // are multiplied by this anyway.
            double pref = x2 * zeldovichIntegrandMG(rr,qq,f1) * gl.w[imu];
                           // For the bias terms, compute U, xi and Ainv (even though in above).
                std::vector<double> qf  =interpQfuncs(qlen);
                std::vector<double> ef  =interpEfuncs(qlen);
                //std::cout  << "here" << std::endl;
                std::vector<double> efMG  =interpEfuncsMG(qlen);
                //std::cout  << "here" << std::endl;
                std::vector<double> Ainv=calcAinv(qq);
                std::vector<double> Ainvprime=calcAinvprime(qq);
                std::vector<double> Ainvprimeprime=calcAinvprimeprime(qq);
                std::vector<double> Aloop(9),Aloop22(9),Aloop13(9),Aloop22p(9),Aloop13p(9),Aloop2p2p(9),Aloop1p3(9),Aloop1p3p(9),Aloop31(9);
                std::vector<double> Alin(9);

                for (int i=0; i<3; ++i) {
                    for (int j=0; j<3; ++j) {
                        Aloop22[3*i+j] = (ef[4])*qh[i]*qh[j]+(ef[1])*(i==j);
                        Aloop13[3*i+j] = (1*ef[5])*qh[i]*qh[j]+(1*ef[2])*(i==j);       
                        Aloop31[3*i+j] = (1*ef[5])*qh[i]*qh[j]+(1*ef[2])*(i==j);                 
                        Alin[ 3*i+j] = ef[3]*qh[i]*qh[j]+ef[0]*(i==j);
                        Aloop22p[3*i+j] = (efMG[4])*qh[i]*qh[j]+(efMG[1])*(i==j);
                        Aloop13p[3*i+j] = (1*efMG[5])*qh[i]*qh[j]+(1*efMG[2])*(i==j);
                        Aloop1p3[3*i+j] = (1*efMG[19])*qh[i]*qh[j]+(1*efMG[17])*(i==j);
                        Aloop2p2p[3*i+j] = (efMG[18])*qh[i]*qh[j]+(efMG[16])*(i==j);
                        Aloop1p3p[3*i+j] = (1*efMG[21])*qh[i]*qh[j]+(1*efMG[20])*(i==j);
                    }
                }
                //Adding RSD contributions to Ainv and Aloop, don't care about Alin for now
                //First, correct zz components of Ainv (Zeldovich)
                 for (int i=0; i<2; ++i) {
                   Ainv[3*i+2] = Ainvprime[3*i+2];
                   Ainv[3*2+i] = Ainvprime[3*2+i];
                 }
                 Ainv[8] = Ainvprimeprime[8]; //added by George        
                 double dett= Ainv[0]*(Ainv[4]*Ainv[8]-Ainv[7]*Ainv[5])- Ainv[1]*(Ainv[3]*Ainv[8]-Ainv[6]*Ainv[5])+Ainv[2]*(Ainv[3]*Ainv[7]-Ainv[6]*Ainv[4]);
                 Ainv[9] = dett;

                //Correct Aloop
                for (int i=0; i<3; ++i) {
                   Alin[3*i+2] *= (1+f1);
                   Alin[3*2+i] *= (1+f1);
                   //Aloop22[3*i+2] *= (1+2*f1);
                   //Aloop22[3*2+i] *= (1+2*f1);
                   Aloop22[3*i+2] += Aloop22p[3*i+2];
                   Aloop22[3*2+i] += Aloop22p[3*2+i];
                   //Aloop13[3*i+2] *= (1+f1);
                   //Aloop13[3*2+i] *= (1+3*f1);
                   Aloop13[3*i+2] += Aloop1p3[3*i+2];
                   Aloop13[3*2+i] += Aloop13p[3*2+i];
                   Aloop31[3*i+2] += Aloop13p[3*i+2];
                   Aloop31[3*2+i] += Aloop1p3[3*2+i];
                  }
                 Aloop22[8] += Aloop2p2p[8];
                 Aloop13[8] += Aloop1p3p[8];
                 Aloop31[8] += Aloop1p3p[8];
                //Add two contributions to Aloop, after we shifted each of them to RSD
                for (int i=0; i<3; ++i) {
                    for (int j=0; j<3; ++j) {
                         Aloop[3*i+j] = Aloop22[3*i+j]+ Aloop13[3*i+j]+ Aloop31[3*i+j] ;                        
                    }
                }


                double xiL=qf[3];
                // Construct the auxilliary matrix/vectors g, G of CLPT Eq. (45)
                // and Gamma of Eq. (75).
                double g[3],UL[3],U[3],Udot[3],U20[3],U11[3],G[9],W[27],ULp[3],Up[3],U20p[3],U11p[3];
                for (int i=0; i<3; ++i) {
                    g[i]=0;
                    for (int j=0; j<3; ++j)
                        g[i] += Ainv[3*i+j]*(qq[j]-rr[j]);
                    UL[i] = ef[11]*qh[i];
                    U[ i] =(ef[11]+ef[12])*qh[i];
                    U20[i]=ef[13]*qh[i];
                    U11[i]=ef[14]*qh[i];

                    ULp[i] = efMG[40]*qh[i];
                    Up[ i] =(efMG[40]+efMG[12])*qh[i];
                    U20p[i]=efMG[13]*qh[i];
                    U11p[i]=efMG[14]*qh[i];
                }

                //UL[2] *= 1+f1;	// Correct U as U->RU.
                //U[2] = UL[2] + (1+3*f1)*(ef[12])*qh[2];
                //U20[2] *= 1+2*f1;
                //U11[2] *= 1+2*f1;
                UL[2] += ULp[2];	// Correct U 
                U[2] = UL[2] + (ef[12]+efMG[12])*qh[2];
                U20[2] += U20p[2];
                U11[2] += U11p[2];

                for (int i=0; i<3; ++i)
                    for (int j=0; j<3; ++j)
                        G[3*i+j]=Ainv[3*i+j]-g[i]*g[j];

                double GA=0;
                for (int i=0; i<3; ++i)
                    for (int j=0; j<3; ++j)
                        GA += Aloop[3*i+j]*G[3*i+j];
                double W112p[27], W11p2[27],W1p1p2[27], W11p2p[27], W1p1p2p[27];
                for (int i=0; i<3; ++i)
                    for (int j=0; j<3; ++j)
                        for (int k=0; k<3; ++k){

                            W[9*i+3*j+k] = ef[8]*qh[i]*(j==k) + ef[8]*qh[j]*(i==k) + ef[9]*qh[k]*(i==j) +ef[10]*qh[i]*qh[j]*qh[k];
                            W112p[9*i+3*j+k] = efMG[8]*qh[i]*(j==k) + efMG[8]*qh[j]*(i==k) + efMG[9]*qh[k]*(i==j) +efMG[10]*qh[i]*qh[j]*qh[k];
                            W11p2[9*i+3*j+k] = efMG[26]*qh[i]*(j==k) + efMG[26]*qh[j]*(i==k) + efMG[27]*qh[k]*(i==j) +efMG[28]*qh[i]*qh[j]*qh[k];
                            W1p1p2[9*i+3*j+k] = efMG[29]*qh[i]*(j==k) + efMG[29]*qh[j]*(i==k) + efMG[30]*qh[k]*(i==j) +efMG[31]*qh[i]*qh[j]*qh[k];
                            W11p2p[9*i+3*j+k] = efMG[32]*qh[i]*(j==k) + efMG[32]*qh[j]*(i==k) + efMG[33]*qh[k]*(i==j) +efMG[34]*qh[i]*qh[j]*qh[k];
                            W1p1p2p[9*i+3*j+k] = efMG[35]*qh[i]*(j==k) + efMG[35]*qh[j]*(i==k) + efMG[36]*qh[k]*(i==j) +efMG[37]*qh[i]*qh[j]*qh[k];
                }
                //Add RSD contribution to Wijk
                for (int i=0; i<3; ++i){
                    for (int j=0; j<3; ++j){
                            //W[9*i+3*j+2] *=(1+f1) ;
                            //W[9*i+3*2+j] *=(1+f1) ;
                            //W[9*2+3*i+j] *=(1+2*f1) ;
                            W[9*i+3*j+2] += W11p2[9*i+3*j+2] ;
                            W[9*i+3*2+j] += W11p2[9*i+3*2+j] ;
                            //W[9*2+3*i+j] += f1*2*W[9*2+3*i+j];
                            W[9*2+3*i+j] += W112p[9*2+3*i+j];
                  }
                }
                for (int i=0; i<3; ++i){
                            W[9*i+3*2+2] += W1p1p2[9*i+3*2+2] ;
                            W[9*2+3*2+i] += W11p2p[9*2+3*2+i] ;
                            W[9*2+3*i+2] += W11p2p[9*2+3*i+2];
                }                
                W[26] += W1p1p2p[26];

                // We also need \ddot{A}^{10}:
                double A10[9],A12p[9],A1p2[9],A1p2p[9];
                for (int i=0; i<3; ++i)
                    for (int j=0; j<3; ++j){
                        A10[3*i+j] = (4*ef[7])*qh[i]*qh[j] + (4*ef[6])*(i==j);
                        A12p[3*i+j] = (4*efMG[7])*qh[i]*qh[j] + (4*efMG[6])*(i==j);
                        A1p2[3*i+j] = (4*efMG[23])*qh[i]*qh[j] + (4*efMG[22])*(i==j);
                        A1p2p[3*i+j] = (4*efMG[25])*qh[i]*qh[j] + (4*efMG[24])*(i==j);
                }

                //Adding RSD contributions to A10
                for (int i=0; i<3; ++i) {
                   //A10[3*i+2] *= (1+f1);
                   //A10[3*2+i] *= (1+2*f1);
                   A10[3*i+2] += A1p2[3*i+2];
                   A10[3*2+i] += A12p[3*2+i];
                  }
                A10[8] +=A1p2p[8];

                double GW=0;
                double V1=ef[8];
                double V3=ef[9];
                double Tq=ef[10];
                for (int i=0; i<3; ++i) {
                    for (int j=0; j<3; ++j) {
                        for (int k=0; k<3; ++k) {
                            double Gam,Wel;
                            Gam = Ainv[3*i+j]*g[k]+Ainv[3*k+i]*g[j]+Ainv[3*j+k]*g[i]
                                - g[i]*g[j]*g[k];
                            //Wel   = Tq*qh[i]*qh[j]*qh[k];
                            //if (j==k) Wel += V1*qh[i];
                            //if (i==k) Wel += V1*qh[j];
                            //if (i==j) Wel += V3*qh[k];
                            Wel = W[9*i+3*j+k];
                            GW += Gam*Wel;
                        }
                    }
                }
                GW *= 3;	// Account for permutations.

                double trG = 0, Ug = 0, ULg = 0, U2 = 0, gq = 0, qG = 0, gA = 0, UGA = 0, qGq = 0, gAL = 0, gAU = 0, AGA = 0;
                for (int i=0; i<3; ++i) {
                    gq += g[i]*qh[i];
                    gAL+= g[i]*Alin[3*2+i];
                    Ug += U[i]*g[i];
                    ULg+=UL[i]*g[i];
                    U2 +=UL[i]*UL[i];
                    qG += qh[i]*G[3*2+i];
                    trG+= G[3*i+i];
                    for (int j=0; j<3; ++j) {
                        UGA += UL[i]*G[3*i+j]*Alin[3*2+j];
                        qGq += qh[i]*G[3*i+j]*qh[j];
                        gAU += g[i]*Alin[3*i+j]*UL[j];
                        AGA += Alin[3*2+i]*G[3*i+j]*Alin[3*2+j];
                    }
                }

                double U20g=0,U11g=0;
                for (int i=0; i<3; ++i) {
                    U20g += U20[i]*g[i];
                    U11g += U11[i]*g[i];
                }
                double UUG=0,qqG=0;
                for (int i=0; i<3; ++i)
                    for (int j=0; j<3; ++j) {
                        UUG += G[3*i+j]*UL[i]*UL[j];
                        qqG += G[3*i+j]*qh[i]*qh[j];
                    }
                double A10G=2*trG*ef[6] + 2*qqG*ef[7];





                // Now do the 1, Fp, Fpp, Fp^2, Fp.Fpp, Fpp^2
                xi[ 0] +=   pref *(1-GA/2.+GW/6.); // +Aeft*trG*eftNorm);
                xi[ 1] +=  -pref *(2*Ug+A10G);
                xi[ 2] +=  -pref *(UUG+U20g);
                xi[ 3] +=   pref *(xiL-UUG-U11g);
                xi[ 4] +=  -pref *(2*xiL*ULg);
                xi[ 5] +=   pref *xiL*xiL/2;


          }
        }
      }
    }
    for (int j=0; j<xi.size(); ++j) {
      xi[j] *= dx*dphi;	// Convert sum to integral.
    }
    xi[0] -= 1.0;	// Calculated 1+xi, subtract 1 for xi.
    return(xi);
  }

std::vector<double> 
LSM::xiloopContributionsRSDmultMG(const double rval, const double f1) {
    // Returns the contributions to the multipoles of the redshift-space
    // correlation function for locally biased tracers.
    // This is not tested for very large or small values of r.
    const int Nmu=4;
    GaussLegendre gg = GaussLegendre(2*Nmu);	// Must be even.
    // For even lengths, can sum over half of the points.
    std::vector<double> xiell;
    try{xiell.resize(12);}catch(std::exception& e) {myexception(e);}
    for (int i=0; i<Nmu; ++i) {
      std::vector<double> ximu = xiloopContributionsRSDMG(rval,gg.x[i],f1);
      double p0=1.0;
      double p2=0.5*(3*gg.x[i]*gg.x[i]-1);
      for (int j=0; j<ximu.size(); ++j) {
        xiell[0*ximu.size()+j] += ximu[j]*gg.w[i] * p0 * 1;
        xiell[1*ximu.size()+j] += ximu[j]*gg.w[i] * p2 * 5;
      }
    }
    return(xiell);
  }

//
void
LSM::init(const char fname[], const double f,
          const double b1, const double b2, const double bs,
          const double Aeft, const double Aeft1, const double Aeft2) {
    // Set up the Zel'dovich class.
    Zeldovich::init(fname);
    // Initialize the LPT class with our newly populated kLin/pLin.
    lpt.init(kLin,pLin);
    //std::cout << kLin[0] << std::endl;
    // We can now set up the "extra" functions we need.  Check
    // to see whether we have this pretabulated.
    std::stringstream ss;
    ss<<fname<<".lesmSave";
    std::ifstream fs(ss.str().c_str());
    //if (!fs) {
      // Set up the Q and R's that we need.
     //setupQR();
     std::cout << 600 << std::endl;
     setupQR(fname);
     std::cout << 600 << std::endl;
     // and tabulate/calculate the "extra" functions.
     tabulateEfuncs();
     std::cout << 600 << std::endl;
     tabulateEfuncsMG();
     std::cout << 600 << std::endl;
     writeSaveFile(ss.str().c_str());
    //}
    //else {
      fs.close();
      //readSaveFile(ss.str().c_str());
    //}
    // Now tabulate the functions we need for the streaming model.
    // The point at zero lag is known analytically.
    std::vector<double>	rrvec,xivec,vvvec,stvec,spvec;
    std::vector<double>	rrvecMG,xivecMG,vvvecMG,stvecMG,spvecMG,xivecmat;
    double rr,xi,zi,vv,st,sp,xiZel,xiZelmat;
    double rrMG,xiMG,ziMG,vvMG,stMG,spMG;
    // Step up to the maximum distance.
    rr=xi=vv=st=sp=0;
    rrMG=xiMG=vvMG=stMG=spMG=xiZel=xiZelmat=0;
    //const double dr=2,rmax=250;
    //const double dr=1,rmax=290;
    const double dr=1,rmax=780;
    do {
      try {
        rrvec.push_back(rr);
      } catch(std::exception& e) {myexception(e);}
      rr += dr;
    } while(rr<rmax);
    try {
      xivec.resize( rrvec.size() );
      vvvec.resize( rrvec.size() );
      stvec.resize( rrvec.size() );
      spvec.resize( rrvec.size() );
      xivecMG.resize( rrvec.size() );
      vvvecMG.resize( rrvec.size() );
      stvecMG.resize( rrvec.size() );
      spvecMG.resize( rrvec.size() );
      xivecmat.resize( rrvec.size() );
    } catch(std::exception& e) {myexception(e);}
#pragma omp parallel for
    for (int i=1; i<rrvec.size(); ++i) {
      std::vector<double> zC = xiContributions(rrvec[i],0);
      //std::vector<double> vC = v12(rrvec[i]);
      //std::cout  << rrvec[rrvec.size()-1] << std::endl;
      auto allPairs = dvsPair(rrvec[i]);
      std::vector<double> xC = allPairs[0];
      std::vector<double> vC = allPairs[1];
      std::vector<double> sC = allPairs[2];
      xiZel = zC[0]+b1*zC[1]+b2*zC[2]+b1*b1*zC[3]+b1*b2*zC[ 4]+b2*b2*zC[ 5];
      xiZelmat = zC[0];

      xi = xC[0]+b1*xC[1]+b2*xC[2]+b1*b1*xC[3]+b1*b2*xC[ 4]+b2*b2*xC[ 5]+
               Aeft*xC[6]+ 0*xC[7]+   bs*xC[8]+b1*bs*xC[ 9]+b2*bs*xC[10]+
              bs*bs*xC[11];
      vv = vC[0]+b1*vC[1]+b2*vC[2]+b1*b1*vC[3]+b1*b2*vC[ 4]+b2*b2*vC[ 5]+
              Aeft1*vC[6]+Aeft2*vC[7]+bs*vC[8]+b1*bs*vC[ 9];
      sp = sC[0]+b1*sC[1]+b2*sC[2]+b1*b1*sC[ 3]+b1*b2*sC[ 4]+b2*b2*sC[ 5]+
           bs*sC[6]+0*sC[7];
      st = sC[8]+b1*sC[9]+b2*sC[10]+b1*b1*sC[11]+b1*b2*sC[12]+b2*b2*sC[13]+
           bs*sC[14]+0*sC[15];
      vv*=   f/(1+xi);
      sp*= f*f/(1+xi);
      st*= f*f/(1+xi);
      //sp*= f*f/(1+xiZel);
      //st*= f*f/(1+xiZel);
      xivec[i] = xi*rrvec[i]*rrvec[i];  // Actually stores r^2.xi
      vvvec[i] = vv;
      spvec[i] = sp;
      stvec[i] = st;
    //MG counterpart starts here, added by George
      //std::cout  << "here" << std::endl;
      auto allPairsMG = dvsPairMG(rrvec[i]);
      std::vector<double> xCMG = allPairsMG[0];
      std::vector<double> vCMG = allPairsMG[1];
      std::vector<double> sCMG = allPairsMG[2];
      //std::cout  << "here" << std::endl;
      //xiZel = zC[0]+b1*zC[1]+b2*zC[2]+b1*b1*zC[3]+b1*b2*zC[ 4]+b2*b2*zC[ 5];
      std::vector<double> qff = interpQfuncs(rrvec[i]);
      //if (rrvec[i] < 204 ){
      xiMG = xCMG[0]+b1*xCMG[1]+b2*xCMG[2]+b1*b1*xCMG[3]+b1*b2*xCMG[ 4]+b2*b2*xCMG[ 5]+
               Aeft*xCMG[6]+ 0*xCMG[7]+   bs*xCMG[8]+b1*bs*xCMG[ 9]+b2*bs*xCMG[10]+
          bs*bs*xCMG[11]; //}
      //else {
       //   xiMG = (1+b1)*(1+b1)*qff[3];
      //}
      //if (rrvec[i] < 204 ){
      vvMG = vCMG[0]+b1*vCMG[1]+b2*vCMG[2]+b1*b1*vCMG[3]+b1*b2*vCMG[ 4]+b2*b2*vCMG[ 5]+
          Aeft1*vC[6]+Aeft2*vC[7]+f*bs*vCMG[8]+f*b1*bs*vCMG[ 9];//}
      //    else {
            //  vvMG = 2*f*(1+b1)*qff[2];
       //     vvMG = 2*(1+b1)*qff[14]*(1+xiMG);
       //   }
      //if (rrvec[i] < 204 ){
      spMG = sCMG[0]+b1*sCMG[1]+b2*sCMG[2]+b1*b1*sCMG[ 3]+b1*b2*sCMG[ 4]+b2*b2*sCMG[ 5]+
          f*f*bs*sCMG[6]+0*sC[7];//}
     // else{
      //  spMG = 2*(sigma2primeprime-qff[16])*(1+xiMG);
      //}
     // if (rrvec[i] < 204 ){
      stMG = sCMG[8]+b1*sCMG[9]+b2*sCMG[10]+b1*b1*sCMG[11]+b1*b2*sCMG[12]+b2*b2*sCMG[13]+
          f*f*bs*sCMG[14]+0*sC[15];//}
     // else{
     //  stMG = 2*(sigma2primeprime-qff[15])*(1+xiMG);
     // }
      vvMG*=   1.0/(1+xiMG); //Careful, we do NOT multiply by power of f in MG. It's inside the integral.
      spMG*= 1.0*1.0/(1+xiMG);
      stMG*= 1.0*1.0/(1+xiMG);
      //spMG*= 1.0*1.0/(1+xiZel);//Zeldovich
      //stMG*= 1.0*1.0/(1+xiZel);//Zeldovich
      //std::cout  << "here" << std::endl;
      xivecMG[i] = xiMG*rrvec[i]*rrvec[i];  // Actually stores r^2.xi
      vvvecMG[i] = vvMG;
      spvecMG[i] = spMG;
      stvecMG[i] = stMG;

      xivecmat[i] = xiZelmat*rrvec[i]*rrvec[i];
      //std::cout  << "here" << std::endl;
    }
    // and fit splines to them for later use.
    xispl.init(rrvec,xivec);
    vvspl.init(rrvec,vvvec);
    stspl.init(rrvec,stvec);
    spspl.init(rrvec,spvec);
    //MG splined ingredients for GSM
    xisplMG.init(rrvec,xivecMG);
    vvsplMG.init(rrvec,vvvecMG);
    stsplMG.init(rrvec,stvecMG);
    spsplMG.init(rrvec,spvecMG);

    xiZelmatspl.init(rrvec,xivecmat);
   // return(xivec);
    //Importing xi, sigma, v12 from simulations and getting splines for them 

    //std::ifstream xiN("/home/astrosun2/gvalogiannis/swot-1.1.0/xilightF5z1_bin1160.out");
    std::ifstream xiN("./xiF6sim160.txt");
    //const int Nksim=30;
    const int Nksim=160;
    std::vector<double> rxisim(Nksim), xisim(Nksim), dum(Nksim);
    for (int i=0; i<Nksim; ++i) {
             std::string sxiN;
             getline(xiN,sxiN);
               std::istringstream(sxiN) >> rxisim[i] >> xisim[i] >> dum[i] >> dum[i];
               xisim[i] *= rxisim[i]*rxisim[i];
              //std::cout << rxisim[i] << xisim[i] << std::endl; 
    }
    xiN.close();
    //And spline the new functions
    xisimspl.init(rxisim,xisim);

    //std::ifstream v12N("/home/astrosun2/gvalogiannis/swot-1.1.0/velF5bin1160_1.out");
    std::ifstream v12N("./velF6sim160.txt");
    //const int Nksim=30;
    std::vector<double> rv12sim(Nksim), v12sim(Nksim);
    for (int i=0; i<Nksim; ++i) {
             std::string sv12N;
             getline(v12N,sv12N);
               std::istringstream(sv12N) >> rv12sim[i] >> v12sim[i] >> dum[i] >> dum[i];
              //std::cout << v12sim[i] << std::endl; 
    }
    v12N.close();
    //And spline the new functions
    v12simspl.init(rv12sim,v12sim);

    //std::ifstream sigpN("/home/astrosun2/gvalogiannis/swot-1.1.0/sigparF5bin1160.out");
    std::ifstream sigpN("./sigparF6sim160.txt");
    //const int Nksim=30;
    std::vector<double> rsigparsim(Nksim), sigparsim(Nksim);
    for (int i=0; i<Nksim; ++i) {
             std::string ssigpN;
             getline(sigpN,ssigpN);
               std::istringstream(ssigpN) >> rsigparsim[i] >> sigparsim[i] >> dum[i] >> dum[i];
              //std::cout << sigparsim[i] << std::endl; 
             //sigparsim[i] -= v12sim[i]*v12sim[i];
    }
    sigpN.close();
    //And spline the new functions
    sigparsimspl.init(rsigparsim,sigparsim);

    //std::ifstream sigperN("/home/astrosun2/gvalogiannis/swot-1.1.0/sigperpF5bin1160.out");
    std::ifstream sigperN("./sigperpF6sim160.txt");
    //const int Nksim=30;
    std::vector<double> rsigperpsim(Nksim), sigperpsim(Nksim);
    for (int i=0; i<Nksim; ++i) {
             std::string ssigperN;
             getline(sigperN,ssigperN);
               std::istringstream(ssigperN) >> rsigperpsim[i] >> sigperpsim[i] >> dum[i] >> dum[i];
              //std::cout << sigperpsim[i] << std::endl; 
    }
    sigperN.close();
    //And spline the new functions
    sigperpsimspl.init(rsigperpsim,sigperpsim);
}

double
LSM::xiRZ(const double R, const double Z, const double s2fog) {
    // The 2D correlation function for the streaming model.
    // Does the integral over the "true" line-of-sight separation
    // using Simpson's rule.
    const double R2=R*R;
    const double ymax=50;
    const int    Ny=500;
    const double hh=2*ymax/Ny;
    int errcnt=0;
    double xi=0;
    // Careful throwing exceptions from threads...
#pragma omp parallel for reduction(+:xi,errcnt)
    for (int i=1; i<Ny; ++i) {
      double yy  = -ymax + i*hh;        // Actually Z-y
      double rr  = sqrt(R2+(Z-yy)*(Z-yy));
      if (errcnt>0 || rr>=xispl.xmax() || rr<=xispl.xmin()) {
        errcnt=1;
      }
      else {
        double mu  = (Z-yy)/rr;
        double xip1= 1.0 + xispl(rr)/rr/rr;
        double xiZelmat1= 1.0 + xiZelmatspl(rr)/rr/rr;
        double vr  = mu*vvspl(rr);
        double expt= yy-vr;
        //double s2  = mu*mu*spspl(rr)+(1-mu*mu)*stspl(rr)-vr*vr + s2fog*(xiZelmat1)/(xip1);
        double s2  = mu*mu*spspl(rr)+(1-mu*mu)*stspl(rr)-vr*vr + s2fog;
        int    wt  = 2+2*(i%2);
        if (s2>0)
          xi += xip1*exp(-0.5*expt*expt/s2)/sqrt(s2) * wt;
      }
    }
    if (errcnt>0) {myexit(1);}
    xi *= hh/3.0 / sqrt(2*M_PI);
    xi -= 1.0;
    return(xi);
}

std::vector<double>
LSM::xiEll(const double ss, const double s2fog,
                        const double Apar, const double Aperp) {
    // The multipoles of the correlation function for the streaming model.
    // Integrates the 2D correlation function using Gauss-Legendre integration.
    //std::vector<double> xiell(2);
    std::vector<double> xiell(3);
    const int Nmu=16;
    GaussLegendre gg = GaussLegendre(2*Nmu);	// Must be even.
    // For even lengths, can sum over half of the points.
    for (int i=0; i<Nmu; ++i) {
      double ximu = xiRZ(ss*sqrt(1-gg.x[i]*gg.x[i])*Aperp,
                         ss*gg.x[i]*Apar,s2fog);
      double p0=1.0;
      double p2=0.5*(3*gg.x[i]*gg.x[i]-1);
      double p4=(1./8.)*(35*gg.x[i]*gg.x[i]*gg.x[i]*gg.x[i]-30*gg.x[i]*gg.x[i]+3);
      xiell[0] += ximu*gg.w[i] * p0 * 1;
      xiell[1] += ximu*gg.w[i] * p2 * 5;
      xiell[2] += ximu*gg.w[i] * p4 * 9;
    }
    return(xiell);
}

double
LSM::xiRZMG(const double R, const double Z, const double s2fog) {
    // The 2D correlation function for the streaming model.
    // Does the integral over the "true" line-of-sight separation
    // using Simpson's rule.
    const double R2=R*R;
    const double ymax=50;
    const int    Ny=500;
    const double hh=2*ymax/Ny;
    int errcnt=0;
    double xi=0;
    // Careful throwing exceptions from threads...
#pragma omp parallel for reduction(+:xi,errcnt)
    for (int i=1; i<Ny; ++i) {
      double yy  = -ymax + i*hh;        // Actually Z-y
      double rr  = sqrt(R2+(Z-yy)*(Z-yy));
      if (errcnt>0 || rr>=xisplMG.xmax() || rr<=xisplMG.xmin()) {
        errcnt=1;
      }
      else {
        double mu  = (Z-yy)/rr;
        double xip1= 1.0 + xisplMG(rr)/rr/rr; 
        double xiZelmat1= 1.0 + xiZelmatspl(rr)/rr/rr;
        double vr  = mu*vvsplMG(rr);
        //double vr  = mu*vvspl(rr);
        double expt= yy-vr;
        //double s2  = mu*mu*spsplMG(rr)+(1-mu*mu)*stsplMG(rr)-vr*vr + s2fog*(xiZelmat1)/(xip1);
        double s2  = mu*mu*spsplMG(rr)+(1-mu*mu)*stsplMG(rr)-vr*vr + s2fog;
        //double s2  = mu*mu*spspl(rr)+(1-mu*mu)*stspl(rr)-vr*vr + s2fog;
        int    wt  = 2+2*(i%2);
        if (s2>0)
          xi += xip1*exp(-0.5*expt*expt/s2)/sqrt(s2) * wt;
      }
    }
    if (errcnt>0) {myexit(1);}
    xi *= hh/3.0 / sqrt(2*M_PI);
    xi -= 1.0;
    return(xi);
}

std::vector<double>
LSM::xiEllMG(const double ss, const double s2fog,
                        const double Apar, const double Aperp) {
    // The multipoles of the correlation function for the streaming model.
    // Integrates the 2D correlation function using Gauss-Legendre integration.
    //std::vector<double> xiell(2);
    std::vector<double> xiell(3);
    const int Nmu=16;
    GaussLegendre gg = GaussLegendre(2*Nmu);	// Must be even.
    // For even lengths, can sum over half of the points.
    for (int i=0; i<Nmu; ++i) {
      double ximu = xiRZMG(ss*sqrt(1-gg.x[i]*gg.x[i])*Aperp,
                         ss*gg.x[i]*Apar,s2fog);
      double p0=1.0;
      double p2=0.5*(3*gg.x[i]*gg.x[i]-1);
      double p4=(1./8.)*(35*gg.x[i]*gg.x[i]*gg.x[i]*gg.x[i]-30*gg.x[i]*gg.x[i]+3);
      xiell[0] += ximu*gg.w[i] * p0 * 1;
      xiell[1] += ximu*gg.w[i] * p2 * 5;
      xiell[2] += ximu*gg.w[i] * p4 * 9;
    }
    return(xiell);
}

double
LSM::xiRZsim(const double R, const double Z, const double s2fog) {
    // The 2D correlation function for the streaming model.
    // Does the integral over the "true" line-of-sight separation
    // using Simpson's rule.
    const double R2=R*R;
    const double ymax=50;
    const int    Ny=500;
    const double hh=2*ymax/Ny;
    int errcnt=0;
    double xi=0, rrmin=1000;
    // Careful throwing exceptions from threads...
#pragma omp parallel for reduction(+:xi,errcnt)
    for (int i=1; i<Ny; ++i) {
      double yy  = -ymax + i*hh;        // Actually Z-y
      double rr  = sqrt(R2+(Z-yy)*(Z-yy));
      if (rr < rrmin){
      rrmin=rr;     }
      if (errcnt>0 || rr>=xisimspl.xmax() || rr<=xisimspl.xmin()) {
        errcnt=1;
        //std::cout << rr << std::endl; 
      }
      else {
        double mu  = (Z-yy)/rr;
        double xip1= 1.0 + xisimspl(rr)/rr/rr;
        //double xip1= 1.0 + xisplMG(rr)/rr/rr; 
        //double xiZelmat1= 1.0 + xiZelmatspl(rr)/rr/rr;
        double off=-0.0;
        double vr  = mu*(v12simspl(rr)+off);
        //double vr  = mu*vvsplMG(rr);
        double expt= yy-vr;
        //double s2  = mu*mu*spsplMG(rr)+(1-mu*mu)*stsplMG(rr)-vr*vr + s2fog;
        double s2  = mu*mu*sigparsimspl(rr)+(1-mu*mu)*sigperpsimspl(rr)-vr*vr + s2fog;
        int    wt  = 2+2*(i%2);
        if (s2>0)
          xi += xip1*exp(-0.5*expt*expt/s2)/sqrt(s2) * wt;
      }
    }
    //std::ofstream storefs("rrmin.txt",std::ios_base::app);
    //storefs << std::scientific << std::setw(20) << std::setprecision(9) << rrmin << std::endl;
    //std::cout << rrmin << std::endl; 
    //storefs.close();

    //if (errcnt>0) {myexit(1);}
    xi *= hh/3.0 / sqrt(2*M_PI);
    xi -= 1.0;
    return(xi);
}

std::vector<double>
LSM::xiEllsim(const double ss, const double s2fog,
                        const double Apar, const double Aperp) {
    // The multipoles of the correlation function for the streaming model.
    // Integrates the 2D correlation function using Gauss-Legendre integration.
    //std::vector<double> xiell(2);
    std::vector<double> xiell(3);
    const int Nmu=16;
    GaussLegendre gg = GaussLegendre(2*Nmu);	// Must be even.
    // For even lengths, can sum over half of the points.
    for (int i=0; i<Nmu; ++i) {
      double ximu = xiRZsim(ss*sqrt(1-gg.x[i]*gg.x[i])*Aperp,
                         ss*gg.x[i]*Apar,s2fog);
      double p0=1.0;
      double p2=0.5*(3*gg.x[i]*gg.x[i]-1);
      double p4=(1./8.)*(35*gg.x[i]*gg.x[i]*gg.x[i]*gg.x[i]-30*gg.x[i]*gg.x[i]+3);
      xiell[0] += ximu*gg.w[i] * p0 * 1;
      xiell[1] += ximu*gg.w[i] * p2 * 5;
      xiell[2] += ximu*gg.w[i] * p4 * 9;
    }
    return(xiell);
}

void
LSM::printzFuncs(const char fbase[]) {
    // Print the "extra" functions.
    std::ostringstream ss;
    ss << fbase << ".zFuncs";
    std::ofstream fs(ss.str().c_str());
    if (!fs) {std::cerr<<"Unable to open file."<<std::endl;myexit(1);}
    fs<<"# q-dependent functions computed for Zeldovich."<<std::endl;
    fs<<"# Order is"<<std::endl;
    fs<<"#  1) q [Mpc/h]"<<std::endl
      <<"#  2) eta_per" <<std::endl
      <<"#  3) eta_par" <<std::endl
      <<"#  4) U^{(1)}" <<std::endl
      <<"#  5) xi_L"    <<std::endl
      <<"#  6) xi_L\'"  <<std::endl
      <<"#  7) d^2 xi_L"<<std::endl
      <<"#  8) J_2"     <<std::endl
      <<"#  9) J_3"     <<std::endl
      <<"# 10) J_4"     <<std::endl
      <<"# 11) V_i^{12}"<<std::endl
      <<"# 12) chi^{12}"<<std::endl
      <<"# 13) zeta"    <<std::endl;
    for (int i=0; i<120; ++i) {
      double q = (i+1.0);
      std::vector<double> qf = interpQfuncs(q);
      fs<<std::scientific<<std::setw(15)<<std::setprecision(5)<<q;
      for (int j=0; j<qf.size(); ++j)
        fs<<std::scientific<<std::setw(15)<<std::setprecision(5)<<qf[j];
      fs<<std::endl;
    }
    fs.close();
}

void
LSM::printqFuncs(const char fbase[]) {
    // Print the "extra" functions.
    std::ostringstream ss;
    ss << fbase << ".qFuncs";
    std::ofstream fs(ss.str().c_str());
    if (!fs) {std::cerr<<"Unable to open file."<<std::endl;myexit(1);}
    fs<<"# q-dependent functions used by CLPT."<<std::endl;
    fs<<"# Order is"<<std::endl;
    fs<<"#  1) q [Mpc/h]"<<std::endl
      <<"#  2) X^{(11)}"<<std::endl
      <<"#  3) X^{(22)}"<<std::endl
      <<"#  4) X^{(13)}"<<std::endl
      <<"#  5) Y^{(11)}"<<std::endl
      <<"#  6) Y^{(22)}"<<std::endl
      <<"#  7) Y^{(13)}"<<std::endl
      <<"#  8) X^{(12)}_{10}"<<std::endl
      <<"#  9) Y^{(12)}_{10}"<<std::endl
      <<"# 10) V^{(112)}_{1}"<<std::endl
      <<"# 11) V^{(112)}_{3}"<<std::endl
      <<"# 12) T^{(112)}"<<std::endl
      <<"# 13) U^{(1)}"<<std::endl
      <<"# 14) U^{(3)}"<<std::endl
      <<"# 15) U^{(2)}_{20}"<<std::endl
      <<"# 16) U^{(2)}_{11}"<<std::endl
      <<"# 17) V^{10}"<<std::endl;
    for (int i=0; i<60; ++i) {
      double q = 2*(i+0.5);
      std::vector<double> ef = interpEfuncs(q);
      fs<<std::scientific<<std::setw(15)<<std::setprecision(5)<<q;
      for (int j=0; j<ef.size(); ++j)
        fs<<std::scientific<<std::setw(15)<<std::setprecision(5)<<ef[j];
      fs<<std::endl;
    }
    fs.close();
}

void
LSM::printXiStuff(const char fbase[]) {
    // Print the contributions to Xi
    std::ostringstream ss;
    ss << "xireal_"<< fbase;
    std::ofstream fs(ss.str().c_str());
    //std::ofstream fs("xiGRz1planck.txt");
    //std::ofstream fs("xiGRDESIZelz05.txt");
    if (!fs) {std::cerr<<"Unable to open file."<<std::endl;myexit(1);}
    fs << "# Contributions to xi_real."<<std::endl
       << "# Order is r [Mpc/h], xi_L,"
       << " 1, b1, b2, b1^2, b1.b2, b2^2, Aeft, d2xiLin, bs, b1.bs, b2.bs, bs^2"
       << std::endl;
    //for (int ir=0; ir<140; ++ir) {
    for (int ir=0; ir<230; ++ir) {
      //double rr = 10.0 + 2*(ir+0.5);
      double rr = 0.0 + 1*(ir+0.5);
      std::vector<double> qf = interpQfuncs(rr); //xiContributions
      std::vector<double> xC = dvsPairMG(rr)[0];
      //std::vector<double> xC = xiContributions(rr,0);
      fs<<std::scientific<<std::setw(12)<<std::setprecision(4)<<rr;
      fs<<std::scientific<<std::setw(12)<<std::setprecision(4)<<xisplMG(rr)/rr/rr;
      fs<<std::scientific<<std::setw(12)<<std::setprecision(4)<<qf[3];
      for (int j=0; j<xC.size(); ++j)
        fs<<std::scientific<<std::setw(12)<<std::setprecision(4)<<xC[j];
      fs<<std::endl;
    }
    fs.close();
}
//xisplMG.init(rrvec,xivecMG);
//vvsplMG.init(rrvec,vvvecMG);
//stsplMG.init(rrvec,stvecMG);
//spsplMG

void
LSM::printXiStuffRSD(const char fbase[]) {
    // Print the contributions to Xi
    std::ostringstream ss;
    ss << fbase << ".xiGRDESIZelz05dir.txt";
    //std::ofstream fs(ss.str().c_str());
    //std::ofstream fs("xiGRz1planck.txt");
    std::ofstream fs("xiGRDESIZelz05dir.txt");
    if (!fs) {std::cerr<<"Unable to open file."<<std::endl;myexit(1);}
    fs << "# Contributions to xi_rsd."<<std::endl
       << "# Order is r [Mpc/h],"
       << " 1, b1, b2, b1^2, b1.b2, b2^2"
       << std::endl;
    for (int ir=0; ir<140; ++ir) {
      //double rr = 10.0 + 2*(ir+0.5);
      double rr = 0.0 + 1*(ir+0.5);
      std::vector<double> qf = interpQfuncs(rr);
      //std::vector<double> xC = dvsPair(rr)[0];
      std::vector<double> xC = xiContributionsRSDmult(rr, 0.720427);
      //std::vector<double> xC = xiloopContributionsRSDmultMG(rr, 0.733327);
      //std::vector<double> ximu = xiContributionsRSDMG(rr, 1, 0);
      fs<<std::scientific<<std::setw(12)<<std::setprecision(4)<<rr;
      //  <<std::scientific<<std::setw(12)<<std::setprecision(4)<<qf[3];
      for (int j=0; j<xC.size(); ++j)
        fs<<std::scientific<<std::setw(12)<<std::setprecision(4)<<xC[j];

      fs<<std::endl;
    }
    fs.close();
}

void
LSM::printVpStuff(const char fbase[]) {
    // Print the contributions to Vpair
    std::ostringstream ss;
    //ss << fbase << ".vplinGRz1.txt";
    ss << "vpreal_"<< fbase;
    //std::ofstream fs(ss.str().c_str());
    std::ofstream fs(ss.str().c_str());
    //std::ofstream fs("vplinGRz1.txt");
    if (!fs) {std::cerr<<"Unable to open file."<<std::endl;myexit(1);}
    fs << "# Contributions to (1/f)(1+xi)v_{12}."<<std::endl
       << "# Order is r [Mpc/h] "
       << "v_L, 1, b1, b2, b1^2, b1.b2, b2^2, gradXi, g_los, bs2, bs2.b1"
       << std::endl;
    //for (int ir=0; ir<140; ++ir) {
    for (int ir=0; ir<230; ++ir) {
      double rr = 0.0 + 1*(ir+0.5);
        
      std::vector<double> qf = interpQfuncs(rr);
      std::vector<double> vC = dvsPairMG(rr)[1];
      //std::vector<double> vC = dvsPair(rr)[1];
      fs<<std::scientific<<std::setw(12)<<std::setprecision(4)<<rr;
      //fs  <<std::scientific<<std::setw(12)<<std::setprecision(4)<<2*qf[2];
      //fs  <<std::scientific<<std::setw(12)<<std::setprecision(4)<<2*qf[14];
      fs<<std::scientific<<std::setw(12)<<std::setprecision(4)<<vvsplMG(rr);
      for (int j=0; j<vC.size(); ++j)
        fs<<std::scientific<<std::setw(12)<<std::setprecision(4)<<vC[j];
      fs<<std::endl;
    }
    fs.close();
}

void
LSM::printS2Stuff(const char fbase[]) {
    // Print the contributions to Sigma^2
    std::ostringstream ss;
    //ss << fbase << ".s2Stuffz05Zel.txt";
    //std::ofstream fs(ss.str().c_str());
    //std::ofstream fs("s2Stuffz05Zel.txt");
    ss << "vpreal_"<< fbase;
    //std::ofstream fs(ss.str().c_str());
    std::ofstream fs(ss.str().c_str());
    if (!fs) {std::cerr<<"Unable to open file."<<std::endl;myexit(1);}
    fs << "# Contributions to (1/f^2)(1+xi)sigma^2."<<std::endl
       << "# Order is r [Mpc/h]"<<std::endl
       << "# 1, b1, b2, b1^2, b1.b2, b2^2, bs2, beta"<<std::endl
       << "# first for sig2_par then for sig2_perp."<<std::endl;
    //for (int ir=0; ir<140; ++ir) {
    for (int ir=0; ir<230; ++ir) {
      double rr = 0.0 + 1*(ir+0.5);
      //std::vector<double> sC = dvsPair(rr)[2];
      std::vector<double> sC = dvsPairMG(rr)[2];
      fs<<std::scientific<<std::setw(12)<<std::setprecision(4)<<rr;
      fs<<std::scientific<<std::setw(12)<<std::setprecision(4)<<stsplMG(rr);
      fs<<std::scientific<<std::setw(12)<<std::setprecision(4)<<spsplMG(rr);
      for (int j=0; j<sC.size(); ++j)
        fs<<std::scientific<<std::setw(12)<<std::setprecision(4)<<sC[j];
      fs<<std::endl;
    }
    fs.close();
}

