/*
 *		An event-driven parser for command-line arguments.
 *  
 *		Copyright (c) 2004-2005 by N.Okazaki
 *
 * This software is provided 'as-is', without any express or implied
 * warranty.  In no event will the authors be held liable for any damages
 * arising from the use of this software.
 *
 * Permission is granted to anyone to use this software for any purpose,
 * including commercial applications, and to alter it and redistribute it
 * freely, subject to the following restrictions (known as zlib license):
 *
 * 1. The origin of this software must not be misrepresented; you must not
 *    claim that you wrote the original software. If you use this software
 *    in a product, an acknowledgment in the product documentation would be
 *    appreciated but is not required.
 * 2. Altered source versions must be plainly marked as such, and must not be
 *    misrepresented as being the original software.
 * 3. This notice may not be removed or altered from any source distribution.
 *
 * Naoaki Okazaki <okazaki at chokkan.org>
 *
 */

/* $Id: optparse.h 2 2006-10-31 00:57:57Z naoaki $ */

/*
 * Class 'optparse' implements a parser for GNU-style command-line arguments.
 * Inherit this class to define your own option variables and to implement an
 * option handler with macros, BEGIN_OPTION_MAP, ON_OPTION(_WITH_ARG), and
 * END_OPTION_MAP. Consult the sample program attached at the bottom of this
 * source code.
 *
 * This code was comfirmed to be compiled with MCVC++ 2003 and gcc 3.3.
 * Define _BUILD_NCL_SAMPLE if you want to build a sample program.
 *	$ g++ -D_BUILD_NCL_SAMPLE -xc++ optparse.h
 */

#ifndef	__INCLUDE_OPTIONPARSER_H
#define	__INCLUDE_OPTIONPARSER_H

#include <cstring>
#include <sstream>
#include <stdexcept>
#include <string>

#ifdef	USE_NCL_NAMESPACE
namespace ncl {
#endif/*USE_NCL_NAMESPACE*/


/**
 * An event-driven parser for command-line arguments.
 *	@author	Naoaki Okazaki
 */
class optparse {
public:
	/**
	 * Exception class for unrecognized options.
	 */
	class unrecognized_option : public std::invalid_argument {
	public:
		unrecognized_option(char shortopt)
			: std::invalid_argument(std::string("-") + shortopt) {}
		unrecognized_option(const std::string& longopt)
			: std::invalid_argument(std::string("--") + longopt) {}
	};

	/**
	 * Exception class for invalid values.
	 */
	class invalid_value : public std::invalid_argument {
	public:
		std::string optionstr;

		invalid_value(const std::string& message)
			: std::invalid_argument(message) {}
		invalid_value(char shortopt, const char *longopt, const std::string& message) :
			std::invalid_argument(message),
			optionstr(
				shortopt ?
					(std::string("-") + shortopt) : 
					(longopt ? (std::string("--") + longopt) : std::string(""))
				) {}

        ~invalid_value() throw (){};
		const std::string& option() const {return optionstr; }
	};

public:
	/** Construct. */
	optparse() {}
	/** Destruct. */
	virtual ~optparse() {}

	/**
	 * Parse options.
	 *	@param	argv		array of null-terminated strings to be parsed
	 *	@param	num_argv	specifies the number, in strings, of the array
	 *	@return				the number of used arguments
	 *	@throws				optparse_exception
	 */
	int parse(char * const argv[], int num_argv)
	{
		int i;
		for (i = 0;i < num_argv;++i) {
			const char *token = argv[i];
			if (*token++ == '-') {
				const char *next_token = (i+1 < num_argv) ? argv[i+1] : "";
				if (!*token) {
					break;	// only '-' was found.
				} else if (*token == '-') {
					const char *arg = std::strchr(++token, '=');
					if (arg) {
						arg++;
					} else {
						arg = next_token;
					}
					int ret = handle_option(0, token, arg);
					if (ret < 0) {
						throw unrecognized_option(token);
					}
					if (arg == next_token) {
						i += ret;
					}
				} else {
					char c;
					while ((c = *token++) != '\0') {
						const char *arg = *token ? token : next_token;
						int ret = handle_option(c, token, arg);
						if (ret < 0) {
							throw unrecognized_option(c);
						}
						if (ret > 0) {
							if (arg == token) {
								token = "";
							} else {
								i++;
							}
						}
					} // while
				} // else (*token == '-') 
			} else {
				break;	// a non-option argument was fonud.
			} 
		} // for (i)

		return i;
	}

protected:
	/**
	 * Option handler
	 *	This function should be overridden by inheritance class.
	 *	@param	c			short option character, 0 for long option
	 *	@param	longname	long option name
	 *	@param	arg			an argument for the option
	 *	@return				0 (success);
							1 (success with use of an argument);
							-1 (failed, unrecognized option)
	 *	@throws				option_parser_exception
	 */
	virtual int handle_option(char c, const char *longname, const char *arg)
	{
		return 0;
	}

	int __optstrcmp(const char *option, const char *longname)
	{
		const char *p = std::strchr(option, '=');
		return p ?
			std::strncmp(option, longname, p-option) :
			std::strcmp(option, longname);
	}
};


/** The begin of inline option map. */
#define	BEGIN_OPTION_MAP_INLINE() \
	virtual int handle_option(char __c, const char *__longname, const char *arg) \
	{ \
		int used_args = 0; \
		if (0) { \

/** Define of option map. */
#define	DEFINE_OPTION_MAP() \
	virtual int handle_option(char __c, const char *__longname, const char *arg);

/** Begin of option map implimentation. */
#define	BEGIN_OPTION_MAP(_Class) \
	int _Class::handle_option(char __c, const char *__longname, const char *arg) \
	{ \
		int used_args = 0; \
		if (0) { \

/** An entry of option map */
#define	ON_OPTION(test) \
			return used_args; \
		} else if (test) { \
			used_args = 0; \

#define	ON_OPTION_WITH_ARG(test) \
			return used_args; \
		} else if (test) { \
			used_args = 1; \

/** The end of option map implementation */
#define	END_OPTION_MAP() \
			return used_args; \
		} \
		return -1; \
	} \

/** A predicator for short options */
#define	SHORTOPT(x)		(__c == x)
/** A predicator for long options */
#define	LONGOPT(x)		(!__c && __optstrcmp(__longname, x) == 0)


#ifdef	USE_NCL_NAMESPACE
};
#endif/*USE_NCL_NAMESPACE*/

#include"FKCommon.h"
#include"Logger.h"

/**
 * A class to store parameters specified by command-line arguments
 */
class FKOptionParser : public optparse {
public:
    enum class SC : size_t { Bethe, DMFTCubic1d, DMFTCubic2d, DMFTCubic3d, DMFTCubic4d, DMFTCubicInfd, DFCubic1d, DFCubic2d, DFCubic3d, DFCubic4d };
	FK::RealType beta ;
	FK::RealType U    ;
	FK::RealType t    ;
	FK::RealType mu   ;
	FK::RealType e_d  ;
	FK::RealType mix  ;
	size_t n_freq;
	size_t n_dual_freq;
	size_t n_iter;
	std::string sc_type;
    SC sc_index;
	std::string help;
    //bool calc_vertex;

	FKOptionParser() : beta(10), U(4.0), t(1.0), mu(2.0), e_d(0.0), mix(1.0), n_freq(1024), n_dual_freq(1024), n_iter(100), sc_type("bethe"), sc_index(SC::Bethe), help("") {}

	BEGIN_OPTION_MAP_INLINE()
		ON_OPTION(SHORTOPT('b') || LONGOPT("beta"))
			beta = std::atof(arg);
            //n_freq = (int)std::fabs(beta+0.5);
			used_args = 1;	// Notify the parser of a consumption of argument.

        ON_OPTION(SHORTOPT('U') || LONGOPT("U"))
			U = std::atof(arg);
            mu = U/2;
            //n_freq = (int)std::fabs(beta+0.5);
			used_args = 1;	// Notify the parser of a consumption of argument.

        ON_OPTION(SHORTOPT('t') || LONGOPT("t"))
			t = std::atof(arg);
            //n_freq = (int)std::fabs(beta+0.5);
			used_args = 1;	// Notify the parser of a consumption of argument.


        ON_OPTION(LONGOPT("mu"))
			mu = std::atof(arg);
			used_args = 1;	// Notify the parser of a consumption of argument.
        
        ON_OPTION(LONGOPT("ed"))
			e_d = std::atof(arg);
			used_args = 1;	// Notify the parser of a consumption of argument.

        ON_OPTION(LONGOPT("mix"))
			mix = std::atof(arg);
			used_args = 1;	// Notify the parser of a consumption of argument.

        ON_OPTION_WITH_ARG(SHORTOPT('n') || LONGOPT("niter"))
			n_iter = std::atoi(arg);
			used_args = 1;	// Notify the parser of a consumption of argument.
			// no need of the notification: used_args variable will be set to 1.


		ON_OPTION_WITH_ARG(SHORTOPT('m') || LONGOPT("nfreq"))
			n_freq = std::atoi(arg);
            n_dual_freq = n_freq;
			used_args = 1;	// Notify the parser of a consumption of argument.
			// no need of the notification: used_args variable will be set to 1.

        ON_OPTION_WITH_ARG(LONGOPT("ndfreq"))
			n_dual_freq = std::atoi(arg);
			used_args = 1;	// Notify the parser of a consumption of argument.
			// no need of the notification: used_args variable will be set to 1.

        ON_OPTION(SHORTOPT('s') || LONGOPT("sc"))
            sc_type = arg;
            if (sc_type == "bethe") sc_index = SC::Bethe;
            else if (sc_type == "dmftcubic1d") sc_index = SC::DMFTCubic1d;
            else if (sc_type == "dmftcubic2d") sc_index = SC::DMFTCubic2d;
            else if (sc_type == "dmftcubic3d") sc_index = SC::DMFTCubic3d;
            else if (sc_type == "dmftcubic4d") sc_index = SC::DMFTCubic4d;
            else if (sc_type == "dmftcubicinfd") sc_index = SC::DMFTCubicInfd;
            else if (sc_type == "dfcubic1d") sc_index = SC::DFCubic1d;
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
            std::cout << "-n     --niter       : Amount of Matsubara frequencies. Default: " << n_iter<< std::endl;
            std::cout << "-h     --help        : Show this help message" << std::endl;
            //std::cout << "--calc_vertex        : Defines whether the program will calculate a vertex or not. Default: false." << std::endl;
            exit(0);

	END_OPTION_MAP()
};

#endif // endif :: #ifndef __INCLUDE_OPTIONPARSER_H
