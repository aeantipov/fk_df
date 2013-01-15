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
    enum class SC : size_t { DFCubic1d, DFCubic2d, DFCubic3d, DFCubic4d };
	FK::RealType beta ;
	FK::RealType U    ;
	FK::RealType t    ;
	FK::RealType mu   ;
	FK::RealType e_d  ;
	FK::RealType mix  ;
	size_t n_freq;
	size_t n_dual_freq;
	size_t n_dmft_iter;
	size_t n_df_iter;
	size_t n_df_sc_iter;
    FK::RealType df_sc_mix;
	std::string sc_type;
    SC sc_index;
    bool extra_ops;
	std::string help;

	FKOptionParserDF() : 
          beta(10), 
          U(4.0), 
          t(1.0), 
          mu(2.0), 
          e_d(0.0),
          mix(1.0), 
          n_freq(1024), 
          n_dual_freq(128), 
          n_dmft_iter(1000), 
          n_df_iter(100), 
          n_df_sc_iter(10), 
          df_sc_mix(1.0),
          sc_type(""), 
          sc_index(SC::DFCubic1d), 
          extra_ops(false),
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

        ON_OPTION_WITH_ARG(LONGOPT("ndmftiter"))
			n_dmft_iter = std::atoi(arg);
			used_args = 1;	// Notify the parser of a consumption of argument.

        ON_OPTION_WITH_ARG(LONGOPT("ndfiter"))
			n_df_iter = std::atoi(arg);
			used_args = 1;	// Notify the parser of a consumption of argument.

        ON_OPTION_WITH_ARG(LONGOPT("ndfsciter"))
			n_df_sc_iter = std::atoi(arg);
			used_args = 1;	// Notify the parser of a consumption of argument.

        ON_OPTION_WITH_ARG(LONGOPT("dfscmix"))
		    df_sc_mix = std::atof(arg);
			used_args = 1;	// Notify the parser of a consumption of argument.

		ON_OPTION_WITH_ARG(SHORTOPT('m') || LONGOPT("nfreq"))
			n_freq = std::atoi(arg);
            n_dual_freq = n_freq;
			used_args = 1;	// Notify the parser of a consumption of argument.
			// no need of the notification: used_args variable will be set to 1.

        ON_OPTION_WITH_ARG(LONGOPT("ndfreq"))
			n_dual_freq = std::atoi(arg);
			used_args = 1;	// Notify the parser of a consumption of argument.
			// no need of the notification: used_args variable will be set to 1.

        ON_OPTION(LONGOPT("extraops"))
            extra_ops = true;
			used_args = 1;	// Notify the parser of a consumption of argument.
			// no need of the notification: used_args variable will be set to 1.



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
            std::cout << "-b     --beta        : The value of inverse temperature. Default: " << beta << std::endl;
            std::cout << "-U     --U           : The value of U. Default: " << U << std::endl;
            std::cout << "-t     --t           : The value of t. Default: " << t << std::endl;
            
            std::cout << "--ed                 : The value of e_d. Default: " << e_d << std::endl;
            std::cout << "--sc                 : The type of self-consistency. Default: " << sc_type << std::endl;
            std::cout << "Possible values: bethe; dmftcubic1d; dmftcubic2d; dmftcubic3d; dmftcubic4d; dmftcubicinfd; dfcubic1d; dfcubic2d; dfcubic3d; dfcubic4d " << sc_type << std::endl;
            std::cout << "-m     --matsubaras  : Amount of Matsubara frequencies. Default: " << n_freq<< std::endl;
            std::cout << "--ndfreq             : Amount of Matsubara frequencies for DF calc. Default: " << n_dual_freq<< std::endl;
            std::cout << "--ndmftiter          : Amount of DMFT iterations. Default: " << n_dmft_iter<< std::endl;
            std::cout << "--ndfiter            : Amount of DF iterations. Default: " << n_df_iter<< std::endl;
            std::cout << "--ndfsciter          : Amount of DF self-consistency iterations. Default: " << n_df_sc_iter<< std::endl;
            std::cout << "-h     --help        : Show this help message" << std::endl;
            //std::cout << "--calc_vertex        : Defines whether the program will calculate a vertex or not. Default: false." << std::endl;
            exit(0);

	END_OPTION_MAP()
};

#endif // endif :: #ifndef __INCLUDE_FK_OPTIONPARSER_DF_H
