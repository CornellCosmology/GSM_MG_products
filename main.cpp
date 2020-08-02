#include	<iostream>
#include	<iomanip>
#include <fstream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "utils.hpp"
#include "lsm.hpp"
#include "zeldovich.hpp"

//#include "omp.h"

//omp_set_num_threads(32);

int	main(int argc, char **argv)
{
  //if (argc != 9) {
  //  std::cout<<"Usage: lesm "
  //           <<"<Pk-file> <ff> <b1> <b2> <bs2> <Aeft> <Aeftv> <s2FoG>"
  //           <<std::endl;
  //  myexit(1);
  //}
    if (argc != 10) {
      std::cout<<"Usage: lesm "
               <<"<Pk-file> <ff> <b1> <b2> <bs2> <Aeft> <Aeftv> <s2FoG> <ngrav>"
               //<<"<Pk-file> <ff> <b1> <b2> <bs2> <Aeft> <Aeftv> <s2FoG>"
               <<std::endl;
      myexit(1);
    }

  const double ff   = atof(argv[2]);
  const double b1   = atof(argv[3]);
  const double b2   = atof(argv[4]);
  const double bs2  = atof(argv[5]);
  const double Apar = 1;
  const double Aperp= 1;
  const double Aeft = atof(argv[6]);
  const double Aeft1= atof(argv[7]);
  const double Aeft2= 0;
  const double s2FoG= atof(argv[8]);
  const double ngrav= atof(argv[9]); //Extra new parameter for determing GR or MG run
  try {
    LSM lsm(argv[1],ff,b1,b2,bs2,Aeft,Aeft1,Aeft2);
    //LSM lsm(argv[1],ff,b1,b2,bs2,Aeft,Aeft1,Aeft2,ngrav); //Adding new MG parameter
    Zeldovich zel(argv[1]);

#ifdef	PRINTSTUFF
    lsm.printzFuncs( argv[1]);
    lsm.printqFuncs( argv[1]);
    //lsm.printXiStuff(argv[1]);
    //lsm.printXiStuffRSD(argv[1]);
    //lsm.printVpStuff(argv[1]);
    //lsm.printS2Stuff(argv[1]);
    return(0);
#endif
    //zel.print_eta();
    lsm.printXiStuff(argv[1]);
    lsm.printVpStuff("textVp.txt");
    lsm.printS2Stuff("textS2.txt");
    //lsm.printXiStuffRSD("text.txt");
    //const int Nr=61; //Default settings
    //const double rmin=5.0,rmax=135;
    //const int Nr=140;
    const int Nr=380;
    //const double rmin=0.5,rmax=139.5;
    const double rexpmin=log(0.001),rexpmax=log(600);
    //const int Nr=35;
    //const double rmin=25.,rmax=599.;
    //const double rexpmin=log(25.),rexpmax=log(599.);

    //const int Nr=25;
    //const double rmin=5,rmax=139;
    //const int Nr=35;
    //const double rmin=25,rmax=139;
    //const int Nr=50;
    //const double rmin=13,rmax=132.5;
    std::vector<double> xiell;
    // new stuff 
    //std::vector<double> xivec;

    std::stringstream ss;
    ss<<"xi0GSM_"<< argv[1];
    //std::cout<<ss.str().c_str()<<std::endl;
    std::ofstream ffr(ss.str().c_str());
    //std::ofstream ffr("xi0GSMDESIZelz05.txt");

    //std::ofstream fs("xiquadGSMGRloopz05.txt");
    //std::ofstream ffr("xi0GSMGRloopz1bin3.txt");
    //std::ofstream fs("xiquadGSMGRloopz1bin3.txt");
    //std::ofstream ffr("xi0ZSMF5z1bin2.txt");
    //std::ofstream fs("xiquadZSMF5z1bin2.txt");
    //std::ofstream ffr("xi0GSMFr5loopz1.txt");
    //std::ofstream fs("xiquadGSMFr5loopz1.txt");
    //std::ofstream ffr("xi0GSMGRloopz1bin3EFT.txt");
    //std::ofstream fs("xiquadGSMGRloopz1bin3EFT.txt");
    //std::ofstream ffr("xi0GSMGRloopz05EFT.txt");
    //std::ofstream fs("xiquadGSMGRloopz05EFT.txt");
    //std::ofstream ffr("xi0GSMFr6loopz05fit.txt");
    //std::ofstream fs("xiquadGSMFr6loopz05fit.txt");
    //std::ofstream ffr("xi0GSMN5loopz05.txt");
    //std::ofstream fs("xiquadGSMN5loopz05.txt");
    //std::ofstream ffr("xi0ZSMN5z05.txt");
    //std::ofstream fs("xiquadZSMN5z05.txt");
    //std::ofstream ffr("xi0ZSMGRloopz1bin2.txt");
    //std::ofstream fs("xiquadZSMGRloopz1bin2.txt");
    //std::ofstream ffr("xi0GSMF6z05sim.txt");
    //std::ofstream fs("xiquadGSMF6z05sim.txt");
    //std::ofstream fs("xihexaGSMF5z1simbin1.txt");
    //std::ofstream fs("xihexaGSMN5z05EFT.txt");
    for (int i=0; i<Nr; ++i) {
      //double rr = rmin + i*(rmax-rmin)/(Nr-1);
      //logarithmic binning
      double rr = rexpmin + i*(rexpmax-rexpmin)/(Nr-1);
      rr = exp(rr);
        if (ngrav==0){
           xiell=lsm.xiEll(rr,s2FoG,Apar,Aperp);
      }
      else if (ngrav==1) {
           xiell=lsm.xiEllMG(rr,s2FoG,Apar,Aperp);
      }
      //xiell=lsm.xiEllMG(rr,s2FoG,Apar,Aperp);
      //xiell=lsm.xiEllsim(rr,s2FoG,Apar,Aperp);
      //xivec=lsm.init("test",ff,b1,b2,bs2,Aeft,Aeft1,Aeft2);
      std::cout<<std::fixed<<std::setw(10)<<std::setprecision(2)<<rr
               <<std::fixed<<std::setw(15)<<std::setprecision(5)<<xiell[0]*rr*rr
               <<std::fixed<<std::setw(15)<<std::setprecision(5)<<xiell[1]*rr*rr
               <<std::fixed<<std::setw(15)<<std::setprecision(5)<<xiell[2]*rr*rr
               <<std::endl;
     // fs<<std::scientific<<std::setw(12)<<std::setprecision(4)<<rr
     //   <<std::scientific<<std::setw(12)<<std::setprecision(4)<<xiell[1]*rr*rr
     //   <<std::endl;
      ffr<<std::scientific<<std::setw(10)<<std::setprecision(4)<<rr
         <<std::scientific<<std::setw(18)<<std::setprecision(10)<<xiell[0]*rr*rr
         <<std::scientific<<std::setw(18)<<std::setprecision(10)<<xiell[1]*rr*rr
         <<std::scientific<<std::setw(18)<<std::setprecision(10)<<xiell[2]*rr*rr
         <<std::endl;
      //fs<<std::scientific<<std::setw(12)<<std::setprecision(4)<<rr
      //  <<std::scientific<<std::setw(12)<<std::setprecision(4)<<xiell[2]*rr*rr
      //  <<std::endl;
    }
  } catch(std::exception& e) {myexception(e);}
  return(0);
}
