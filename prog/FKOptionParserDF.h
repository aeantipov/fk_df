#ifndef	__INCLUDE_FK_OPTIONPARSER_DF_H
#define	__INCLUDE_FK_OPTIONPARSER_DF_H

#include"OptionParser.h"

/**
 * A class to store parameters specified by command-line arguments
 */
class FKOptionParserDF : public optparse {
public:
	FK::real_type beta = 10.0;
	FK::real_type U    = 4.0 ;
	FK::real_type t    = 1.0 ;
	FK::real_type tp   = 0.0 ;
	FK::real_type mu   = 2.0 ;
	FK::real_type e_d  = 0.0 ;
	FK::real_type w_0  = 0.5 ;
	FK::real_type w_1  = 0.5 ;
    bool update_weights = true;
    bool read_weights   = false;
	FK::real_type mix  = 1.0;
	size_t n_freq = 96;
	size_t kpts = 16;
	size_t n_dual_freq = 12;
	size_t NDMFTRuns = 1000;
	size_t NDFRuns = 1;
    FK::real_type DFCutoff = 1e-8;
	size_t DFNumberOfSelfConsistentIterations = 200;
    FK::real_type DFSCMixing = 0.8;
    FK::real_type DFSCCutoff = 1e-7;
	size_t DFNumberOfBSIterations = 1;
    bool DFEvaluateBSSelfConsistent = false;
    bool DFEvaluateStaticDiagrams = true;
    bool DFEvaluateDynamicDiagrams = false;
    FK::real_type DFBSMixing = 1.0;
    size_t extraops = 0;
    bool update_mixing = true;
    bool read_dual_sigma = false;
    std::string dual_sigma_file = "SigmaDwk.dat";
	std::string help = "";

	FKOptionParserDF() = default; 

	BEGIN_OPTION_MAP_INLINE()
		ON_OPTION_WITH_ARG(SHORTOPT('b') || LONGOPT("beta"))
			beta = std::atof(arg);
            //n_freq = (int)std::fabs(beta+0.5);
			used_args = 1;	// Notify the parser of a consumption of argument.

		ON_OPTION_WITH_ARG(SHORTOPT('T') || LONGOPT("T"))
			beta = 1.0/std::atof(arg);
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

        ON_OPTION_WITH_ARG(LONGOPT("tp"))
			tp = std::atof(arg);
			used_args = 1;	// Notify the parser of a consumption of argument.

        ON_OPTION_WITH_ARG(LONGOPT("mu"))
			mu = std::atof(arg);
			used_args = 1;	// Notify the parser of a consumption of argument.
        
        ON_OPTION_WITH_ARG(LONGOPT("ed"))
			e_d = std::atof(arg);
			used_args = 1;	// Notify the parser of a consumption of argument.

        ON_OPTION_WITH_ARG(LONGOPT("w_0"))
			w_0 = std::atof(arg);
			w_1 = 1.0-w_0;
			used_args = 1;	// Notify the parser of a consumption of argument.

        ON_OPTION_WITH_ARG(LONGOPT("update_weights"))
			update_weights = std::atoi(arg);
			used_args = 1;	// Notify the parser of a consumption of argument.

        ON_OPTION_WITH_ARG(LONGOPT("read_weights"))
			read_weights = std::atoi(arg);
			used_args = 1;	// Notify the parser of a consumption of argument.

        ON_OPTION_WITH_ARG(LONGOPT("mix"))
			mix = std::atof(arg);
			used_args = 1;	// Notify the parser of a consumption of argument.

        ON_OPTION_WITH_ARG(SHORTOPT('k') || LONGOPT("kpoints"))
			kpts = std::atoi(arg);
			used_args = 1;	// Notify the parser of a consumption of argument.
			// no need of the notification: used_args variable will be set to 1.

        ON_OPTION_WITH_ARG(LONGOPT("dfcutoff"))
		    DFCutoff = std::atof(arg);
			used_args = 1;	// Notify the parser of a consumption of argument.

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
        
        ON_OPTION_WITH_ARG(LONGOPT("dfsccutoff"))
		    DFSCCutoff = std::atof(arg);
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

        ON_OPTION_WITH_ARG(LONGOPT("dfevalbssc"))
            DFEvaluateBSSelfConsistent = std::atoi(arg);
            used_args = 1;

        ON_OPTION_WITH_ARG(LONGOPT("dfevalstatic"))
            DFEvaluateStaticDiagrams = std::atoi(arg);
            used_args = 1;
        ON_OPTION_WITH_ARG(LONGOPT("dfevaldynamic"))
            DFEvaluateDynamicDiagrams = std::atoi(arg);
            used_args = 1;

        ON_OPTION_WITH_ARG(LONGOPT("sigmad"))
            read_dual_sigma = true;
            dual_sigma_file = arg;
            used_args = 1;

        ON_OPTION_WITH_ARG(LONGOPT("readsigmad"))
            read_dual_sigma = std::atoi(arg);
            used_args = 1;

        ON_OPTION(SHORTOPT('h') || LONGOPT("help"))
            std::cout << "Usage: fk_DF [options]" << std::endl;
            std::cout << "Options: " << std::endl;
            std::cout << "-h     --help        : Show this help message" << std::endl;
            std::cout << "-b     --beta        : The value of inverse temperature. Default: " << beta << std::endl;
            std::cout << "-U     --U           : The value of U. Default: " << U << std::endl;
            std::cout << "-t     --t           : The value of t. Default: " << t << std::endl;
            std::cout << "--ed                 : The value of e_d. Default: " << e_d << std::endl;
            std::cout << "-m     --matsubaras  : Amount of Matsubara frequencies. Default: " << n_freq<< std::endl;
            std::cout << "--ndfreq             : Amount of bosonic Matsubara frequencies for DF calc. Default: " << n_dual_freq<< std::endl;
            std::cout << "--ndmftiter          : Total amount of DMFT iterations. Default: " << NDMFTRuns<< std::endl;
            std::cout << "--ndfiter            : Total amount of DF iterations. Default: " << NDFRuns<< std::endl;
            std::cout << "DF related:" << std::endl;
            std::cout << "--dfscmix            : Mixing for DF updates of dual Green's function. Default: " << DFSCMixing << std::endl;
            std::cout << "--dfsccutoff         : Threshold for DF updating. Default: " << DFSCCutoff << std::endl;
            std::cout << "--ndfsciter          : Amount of DF self-consistency iterations. Default: " << DFNumberOfSelfConsistentIterations<< std::endl;
            std::cout << "--dfevalbssc         : Evaluate BS equation self-consistently. Default: " << std::boolalpha << DFEvaluateBSSelfConsistent << std::endl;
            std::cout << "--ndfbsiter          : Amount of DF Bethe-Salpeter iterations (if needed). Default: " << DFNumberOfBSIterations<< std::endl;
            std::cout << "--dfbsmix            : Mixing for Bethe-Salpeter iterations (if needed). Default: " << DFBSMixing << std::endl;
            std::cout << "--dfevalstatic       : Evaluate static DF diagram series. Default: " << std::boolalpha << DFEvaluateStaticDiagrams << std::endl;
            std::cout << "--dfevaldynamic      : Evaluate dynamic DF diagram series. Default: " << std::boolalpha << DFEvaluateDynamicDiagrams << std::endl;
            std::cout << "--sigmad             : Read dual self-energy from file. If not specified - start from zero. Default: " << dual_sigma_file << std::endl;
            std::cout << "--readsigmad         : Flag to read dual self-energy. Default: " << std::boolalpha << read_dual_sigma << std::endl;
            exit(0);

	END_OPTION_MAP()
};

#endif // endif :: #ifndef __INCLUDE_FK_OPTIONPARSER_DF_H
