#ifndef	__INCLUDE_FK_OPTIONPARSER_DMFT_H
#define	__INCLUDE_FK_OPTIONPARSER_DMFT_H

#include"OptionParser.h"

/**
 * A class to store parameters specified by command-line arguments
 */
class FKOptionParserDMFT : public optparse {
public:
	FK::RealType beta = 10.  ;
	FK::RealType T    = 0.1  ;
	FK::RealType U    = 4.0  ;
	FK::RealType t    = 1.0  ;
	FK::RealType tp   = 0.0  ;
	FK::RealType mu   = 2.0  ;
	FK::RealType e_d  = 0.0  ;
	FK::RealType w_0  = 0.5  ;
	FK::RealType w_1  = 0.5  ;
    bool update_weights = true;
	FK::RealType mix = 1.0;
    FK::RealType cutoff = 100*std::numeric_limits<FK::RealType>::epsilon();
	size_t n_freq = 256;
	size_t kpts = 32;
	size_t n_iter = 1000;
    size_t extra_ops = 0;
	std::string help = "";

	FKOptionParserDMFT(){}

	BEGIN_OPTION_MAP_INLINE()
        ON_OPTION_WITH_ARG(SHORTOPT('b') || LONGOPT("beta"))
			beta = std::atof(arg);
            T = 1.0/beta;
			used_args = 1;	// Notify the parser of a consumption of argument.

		ON_OPTION_WITH_ARG(SHORTOPT('T') || LONGOPT("T"))
			T = std::atof(arg);
            beta = 1.0/T;
			used_args = 1;	// Notify the parser of a consumption of argument.


        ON_OPTION_WITH_ARG(SHORTOPT('U') || LONGOPT("U"))
			U = std::atof(arg);
            mu = U/2;
            //n_freq = (int)std::fabs(beta+0.5);
			used_args = 1;	// Notify the parser of a consumption of argument.

        ON_OPTION_WITH_ARG(SHORTOPT('t') || LONGOPT("t"))
			t = std::atof(arg);
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

        ON_OPTION_WITH_ARG(LONGOPT("mix"))
			mix = std::atof(arg);
			used_args = 1;	// Notify the parser of a consumption of argument.

        ON_OPTION_WITH_ARG(LONGOPT("cutoff"))
			cutoff = std::atof(arg);
			used_args = 1;	// Notify the parser of a consumption of argument.

        ON_OPTION_WITH_ARG(SHORTOPT('n') || LONGOPT("niter"))
			n_iter = std::atoi(arg);
			used_args = 1;	// Notify the parser of a consumption of argument.
			// no need of the notification: used_args variable will be set to 1.

        ON_OPTION_WITH_ARG(SHORTOPT('k') || LONGOPT("kpoints"))
			kpts = std::atoi(arg);
			used_args = 1;	// Notify the parser of a consumption of argument.
			// no need of the notification: used_args variable will be set to 1.


       	ON_OPTION_WITH_ARG(SHORTOPT('m') || LONGOPT("nfreq"))
			n_freq = std::atoi(arg);
			used_args = 1;	// Notify the parser of a consumption of argument.
			// no need of the notification: used_args variable will be set to 1.

        ON_OPTION_WITH_ARG(LONGOPT("extraops"))
			extra_ops = std::atoi(arg);
			used_args = 1;	// Notify the parser of a consumption of argument.
			// no need of the notification: used_args variable will be set to 1.

        ON_OPTION(SHORTOPT('h') || LONGOPT("help"))
            std::cout << "Usage: fk_DF [options]" << std::endl;
            std::cout << "Options: " << std::endl;
            std::cout << "-b     --beta        : The value of inverse temperature. Default: " << beta << std::endl;
            std::cout << "-U     --U           : The value of U. Default: " << U << std::endl;
            std::cout << "--mu                 : The value of mu. Default: " << U << std::endl;
            std::cout << "-t     --t           : The value of t. Default: " << t << std::endl;
            std::cout << "--ed                 : The value of e_d. Default: " << e_d << std::endl;
            std::cout << "-m     --matsubaras  : Amount of Matsubara frequencies. Default: " << n_freq<< std::endl;
            std::cout << "-n     --niter       : Amount of iterations. Default: " << n_iter<< std::endl;
            std::cout << "--extraopt           : Generate additional statistics" << std::endl;
            std::cout << "-h     --help        : Show this help message" << std::endl;
            //std::cout << "--calc_vertex        : Defines whether the program will calculate a vertex or not. Default: false." << std::endl;
            exit(0);

	END_OPTION_MAP()
};

#endif // endif :: #ifndef __INCLUDE_FK_OPTIONPARSER_DMFT_H
