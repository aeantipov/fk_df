#ifndef	__INCLUDE_FK_OPTIONPARSER_DF_H
#define	__INCLUDE_FK_OPTIONPARSER_DF_H

#include"FKCommon.h"
#include"Logger.h"
#include"OptionParser.h"

/**
 * A class to store parameters specified by command-line arguments
 */
class FKOptionParserDF : public optparse {
public:
    enum class SC : size_t { DFCubic1d, DFCubic2d, DFCubic3d, DFCubic4d, Fail };
	FK::RealType beta ;
	FK::RealType U    ;
	FK::RealType t    ;
	FK::RealType mu   ;
	FK::RealType e_d  ;
	FK::RealType mix  ;
	size_t n_freq;
	size_t kpts;
	size_t n_dual_freq;
	size_t NDMFTRuns;
	size_t NDFRuns;
	size_t DFNumberOfSelfConsistentIterations;
    FK::RealType DFSCMixing;
	size_t DFNumberOfBSIterations;
    bool DFEvaluateBSSelfConsistent;
    FK::RealType DFBSMixing;
	std::string sc_type;
    SC sc_index;
    size_t extraops;
	std::string help;

	FKOptionParserDF() : 
          beta(10), 
          U(4.0), 
          t(1.0), 
          mu(2.0), 
          e_d(0.0),
          mix(1.0), 
          n_freq(1024), 
          kpts(32),
          n_dual_freq(128), 
          NDMFTRuns(1000), 
          NDFRuns(100), 
          DFNumberOfSelfConsistentIterations(10), 
          DFSCMixing(0.5),
          DFNumberOfBSIterations(1), 
          DFEvaluateBSSelfConsistent(false),
          DFBSMixing(1.0),
          sc_type(""), 
          sc_index(SC::Fail), 
          extraops(0),
          help("") 
          {}

	BEGIN_OPTION_MAP_INLINE()
		ON_OPTION_WITH_ARG(SHORTOPT('b') || LONGOPT("beta"))
			beta = std::atof(arg);
            //n_freq = (int)std::fabs(beta+0.5);
			used_args = 1;	// Notify the parser of a consumption of argument.

        ON_OPTION_WITH_ARG(SHORTOPT('U') || LONGOPT("U"))
			U = std::atof(arg);
            mu = U/2;
            //n_freq = (int)std::fabs(beta+0.5);
			used_args = 1;	// Notify the parser of a consumption of argument.

        ON_OPTION_WITH_ARG(SHORTOPT('t') || LONGOPT("t"))
			t = std::atof(arg);
            //n_freq = (int)std::fabs(beta+0.5);
			used_args = 1;	// Notify the parser of a consumption of argument.


        ON_OPTION_WITH_ARG(LONGOPT("mu"))
			mu = std::atof(arg);
			used_args = 1;	// Notify the parser of a consumption of argument.
        
        ON_OPTION_WITH_ARG(LONGOPT("ed"))
			e_d = std::atof(arg);
			used_args = 1;	// Notify the parser of a consumption of argument.

        ON_OPTION_WITH_ARG(LONGOPT("mix"))
			mix = std::atof(arg);
			used_args = 1;	// Notify the parser of a consumption of argument.

        ON_OPTION_WITH_ARG(SHORTOPT('k') || LONGOPT("kpoints"))
			kpts = std::atoi(arg);
			used_args = 1;	// Notify the parser of a consumption of argument.
			// no need of the notification: used_args variable will be set to 1.

        ON_OPTION_WITH_ARG(LONGOPT("ndmftiter"))
			NDMFTRuns = std::atoi(arg);
			used_args = 1;	// Notify the parser of a consumption of argument.

        ON_OPTION_WITH_ARG(LONGOPT("ndfiter"))
			NDFRuns = std::atoi(arg);
			used_args = 1;	// Notify the parser of a consumption of argument.

        ON_OPTION_WITH_ARG(LONGOPT("ndfsciter"))
			DFNumberOfSelfConsistentIterations = std::atoi(arg);
			used_args = 1;	// Notify the parser of a consumption of argument.

        ON_OPTION_WITH_ARG(LONGOPT("dfscmix"))
		    DFSCMixing = std::atof(arg);
			used_args = 1;	// Notify the parser of a consumption of argument.
        
        ON_OPTION_WITH_ARG(LONGOPT("ndfbsiter"))
            DFNumberOfBSIterations = std::atoi(arg);
            used_args = 1;

        ON_OPTION_WITH_ARG(LONGOPT("dfbsmix"))
            DFBSMixing = std::atof(arg);
            used_args = 1;

		ON_OPTION_WITH_ARG(SHORTOPT('m') || LONGOPT("nfreq"))
			n_freq = std::atoi(arg);
            n_dual_freq = n_freq;
			used_args = 1;	// Notify the parser of a consumption of argument.
			// no need of the notification: used_args variable will be set to 1.

        ON_OPTION_WITH_ARG(LONGOPT("ndfreq"))
			n_dual_freq = std::atoi(arg);
			used_args = 1;	// Notify the parser of a consumption of argument.
			// no need of the notification: used_args variable will be set to 1.

        ON_OPTION_WITH_ARG(LONGOPT("extraops"))
            extraops = std::atoi(arg);
			used_args = 1;	// Notify the parser of a consumption of argument.

        ON_OPTION(LONGOPT("dfevalbssc"))
            DFEvaluateBSSelfConsistent = true;
            used_args = 1;

        ON_OPTION(SHORTOPT('s') || LONGOPT("sc"))
            sc_type = arg;
            if (sc_type == "dfcubic1d") sc_index = SC::DFCubic1d;
            else if (sc_type == "dfcubic2d") sc_index = SC::DFCubic2d;
            else if (sc_type == "dfcubic3d") sc_index = SC::DFCubic3d;
            else if (sc_type == "dfcubic4d") sc_index = SC::DFCubic4d;
            used_args = 1;

        ON_OPTION(SHORTOPT('h') || LONGOPT("help"))
            std::cout << "Usage: fk_DF [options]" << std::endl;
            std::cout << "Options: " << std::endl;
            std::cout << "-h     --help        : Show this help message" << std::endl;
            std::cout << "-b     --beta        : The value of inverse temperature. Default: " << beta << std::endl;
            std::cout << "-U     --U           : The value of U. Default: " << U << std::endl;
            std::cout << "-t     --t           : The value of t. Default: " << t << std::endl;
            std::cout << "--ed                 : The value of e_d. Default: " << e_d << std::endl;
            std::cout << "--sc                 : The type of self-consistency. Default: " << sc_type << std::endl;
            std::cout << "Possible values: dfcubic1d; dfcubic2d; dfcubic3d; dfcubic4d " << sc_type << std::endl;
            std::cout << "-m     --matsubaras  : Amount of Matsubara frequencies. Default: " << n_freq<< std::endl;
            std::cout << "--ndfreq             : Amount of bosonic Matsubara frequencies for DF calc. Default: " << n_dual_freq<< std::endl;
            std::cout << "--ndmftiter          : Total amount of DMFT iterations. Default: " << NDMFTRuns<< std::endl;
            std::cout << "--ndfiter            : Total amount of DF iterations. Default: " << NDFRuns<< std::endl;
            std::cout << "DF related:" << std::endl;
            std::cout << "--dfscmix            : Mixing for DF updates of dual Green's function. Default: " << DFSCMixing << std::endl;
            std::cout << "--ndfsciter          : Amount of DF self-consistency iterations. Default: " << DFNumberOfSelfConsistentIterations<< std::endl;
            std::cout << "--dfevalbssc         : Evaluate BS equation self-consistently. Default: " << std::boolalpha << DFEvaluateBSSelfConsistent << std::endl;
            std::cout << "--ndfbsiter          : Amount of DF Bethe-Salpeter iterations (if needed). Default: " << DFNumberOfBSIterations<< std::endl;
            std::cout << "--dfbsmix            : Mixing for Bethe-Salpeter iterations (if needed). Default: " << DFBSMixing << std::endl;
            //std::cout << "--calc_vertex        : Defines whether the program will calculate a vertex or not. Default: false." << std::endl;
            exit(0);

	END_OPTION_MAP()
};

#endif // endif :: #ifndef __INCLUDE_FK_OPTIONPARSER_DF_H
