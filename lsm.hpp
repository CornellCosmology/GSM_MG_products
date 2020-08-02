#ifndef CLEFT_LSM_HPP
#define CLEFT_LSM_HPP

#include "lpt.hpp"
#include "zeldovich.hpp"
#include "spline.hpp"
#include "gauss_legendre.hpp"

// Implements the Lagrangian streaming model.
// Currently set up to take a linear power spectrum, f, F1 and F2
// at input, but could be modified to allow f, F1 and F2 to vary.
class LSM: public Zeldovich {
public:
    LSM() {}	// Do nothing.
    LSM(const char fname[], const double f,
      const double b1, const double b2, const double bs,
      const double Aeft, const double Aeft1, const double Aeft2) {
    init(fname,f,b1,b2,bs,Aeft,Aeft1,Aeft2);
    }

    void init(const char fname[], const double f,
            const double b1, const double b2, const double bs,
            const double Aeft, const double Aeft1, const double Aeft2);

    double xiRZ(const double R, const double Z, const double s2fog);
    double xiRZMG(const double R, const double Z, const double s2fog);
    double xiRZsim(const double R, const double Z, const double s2fog);

    std::vector<double> xiEll(const double ss, const double s2fog,
                            const double Apar, const double Aperp);
    std::vector<double> xiEllMG(const double ss, const double s2fog,
                            const double Apar, const double Aperp);
    std::vector<double> xiEllsim(const double ss, const double s2fog,
                            const double Apar, const double Aperp);
    std::vector<double> xiloopContributionsRSD(const double rval, const double mu, const double f1);
    std::vector<double> xiloopContributionsRSDmult(const double rval, const double f1);
    std::vector<double> xiloopContributionsRSDMG(const double rval, const double mu, const double f1);
    std::vector<double> xiloopContributionsRSDmultMG(const double rval, const double f1);

    void printzFuncs(const char fbase[]);
    void printqFuncs(const char fbase[]);
    void printXiStuff(const char fbase[]);
    void printXiStuffRSD(const char fbase[]);
    void printVpStuff(const char fbase[]);
    void printS2Stuff(const char fbase[]);

protected:
    Spline		xispl,vvspl,stspl,spspl,xiZelmatspl;
    Spline		xisplMG,vvsplMG,stsplMG,spsplMG, xisimspl, v12simspl, sigparsimspl,sigperpsimspl;
    Spline		R1spl,R2spl,Q1spl,Q2spl,Q5spl,Q8spl,Qsspl,QIspl,RIspl,R12spl;
    Spline		R1primespl,R2primespl,Q1primespl,Q2primespl,Q5primespl,Q8primespl,Qprimesspl,QIprimespl,RIprimespl,R12primespl, fgrowthspl,Q1halfprimespl;
    std::vector<double>	X11,X22,X13,Y11,Y22,Y13,X1210,Y1210;
    std::vector<double>	X11MG, X22p, X13p, Y11MG, Y22p, Y13p, X12p10, Y12p10, V1112p, V3112p, TT112p, U1MG, U3p,U220p, U211p,    S2Dp,    X2p2p,   X1p3,   Y2p2p,  Y1p3,   X1p3p, Y1p3p, Y1p,  X1p210, Y1p210, X1p2p10, Y1p2p10, Y1p2p1, V111p2, V311p2,  TT11p2, V11p1p2, V31p1p2,  TT1p1p2, V111p2p, V311p2p,  TT11p2p, V11p1p2p, V31p1p2p, TT1p1p2p, X11p, Y11p, U1p, X1p1p, Y1p1p;
    std::vector<double>	V1,V3,TT;
    std::vector<double>	U1,U3,U220,U211,S2D;
    LPT			lpt;

    std::vector<double> calcEfuncs(const double q);
    std::vector<double> calcEfuncsMG(const double q);

    void tabulateEfuncs();
    void tabulateEfuncsMG();
    void writeSaveFile(const char fname[]);
    void readSaveFile(const char fname[]);
    std::vector<double> interpEfuncs(const double q);
    std::vector<double> interpEfuncsMG(const double q);
    //void setupQR();
    void setupQR(const char fname[]);

    // dpair, vpair, spair.
    std::vector<std::vector<double>> dvsPair(const double rval);
    std::vector<std::vector<double>> dvsPairMG(const double rval);
};

#endif






