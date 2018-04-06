/*
 *  Arguments.cpp
 *  Trogdor6
 *
 *  Created by Paul Hansen on 8/11/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#include "Arguments.h"
#include "Log.h"
#include "Precision.h"
#include "Version.h"

using namespace std;

namespace po = boost::program_options;

po::variables_map
handleArguments(int argc, char* const argv[])
{
	// Declare the supported options.
	// Options allowed only on the command line
	po::options_description generic("Generic");
	generic.add_options()
		("help", "produce help message")
		("version,v", "print complete version information")
		("numerics", "print numerical environment information")
	;
	
	// Options allowed on the command line or in a config file
	po::options_description config("Configuration");
	config.add_options()
        ("adjoint,a", "run in adjoint mode")
        ("noAveraging", "force each cell to take an unmixed permittivity value")
        ("materials", po::value<string>()->default_value("mixed"),
            "harmonic, arithmetic or mixed material averaging")
		("timesteps,t", po::value<int>(), "override number of timesteps")
		("geometry,g", "write 3D geometry file for each grid to an m-file")
        ("oSmoothing", po::value<int>(), "smoothing cells for orientation")
        ("fSmoothing", po::value<int>(), "smoothing cells for fill factors")
        ("derivDelta", po::value<double>(), "step size for internal finite differences")
        ("noNormalizedOrientation", "don't normalize orientation matrices to unit trace")
        ("fullTensor", "use full permittivity tensor")
        ("flipOffDiagonalSign", "multiply off-diagonal permittivity elements by -1")
        ("excludeCorners", "exclude corners and diagonal faces from sensitivity calculation")
        ("offDiagonalOctant", po::value<string>()->default_value("i"),
            "octant location of off-diagonal permittivities (i, j, 0, or 7)")
        ("sensitivity", "write permittivity sensitivities to file")
        ("savePermittivity", "write all update coefficients to file")
        ("saveOrientation", "write all orientation matrices to file")
        ("saveOrientationSensitivity", "write all orientation sensitivities to file")
        ("saveFillFactors", "write all fill factors to file")
        ("saveFillFactorSensitivity", "write all fill factor sensitivities to file")
        ("noPMLE", "turn off the PML correction for E fields")
        ("noPMLH", "turn off the PML correction for H fields")
        ("printRunlines", "duh")
        ("outputDirectory", po::value<string>()->default_value("output"),
            "place to put all the simulation outputs")
        ("outputAuxFields", "save E and H in all auto-generated grids")
        ("outputAuxCurrents", "save J and M in all auto-generated grids")
		("nosim", "do not run simulation")
		("nodragon", "don't draw the dragon")
        ("fastaxis,f", po::value<char>()->default_value('x'),
            "axis along which memory is allocated")
        ("performance", "save detailed performance logfile")
        ("booleans", "use Boolean operations")
	;
	
	// Invisible options
	po::options_description hidden("Allowed options");
	hidden.add_options()
		("input-file", 
			po::value<string>()->default_value("params.xml"),
			"simulation description file (XML)")
	;
	
	// Group the options into pertinent sets
	po::options_description cmdlineOptions;
	cmdlineOptions.add(generic).add(config).add(hidden);
	
	po::options_description configFileOptions;
	configFileOptions.add(config).add(hidden);
	
	po::options_description visibleOptions("Allowed options");
	visibleOptions.add(generic).add(config);
	
	po::positional_options_description positional;
	positional.add("input-file", -1);
	
	po::variables_map variablesMap;
	
	try {
		po::store(po::command_line_parser(argc, const_cast<char**>(argv)).
			options(cmdlineOptions).positional(positional).run(), variablesMap);
	}
	catch (exception & e)
	{
		LOGMORE << "** There was an error parsing the command line.\n";
		LOGMORE << e.what() << endl;
		LOGMORE << "Use --help to see allowed options." << endl;
		exit(1);
	}
	
	try {
		ifstream ifs("trogdor.cfg");
		po::store(po::parse_config_file(ifs, configFileOptions), variablesMap);
	}
	catch (exception & e)
	{
		LOGMORE << "** There was an error parsing the configuration file.\n";
		LOGMORE << e.what() << endl;
		LOGMORE << "Use trogdor --help to see allowed options." << endl;
		exit(1);
	}
	
	po::notify(variablesMap);
	
	LOGF << "Command line invocation: " << endl;
	for (int nn = 0; nn < argc; nn++)
		LOGFMORE << argv[nn] << " ";
	LOGFMORE << "\n";
	LOGF << "(no more command line options)" << endl;
	// to do: dump the config file here
	LOGF << "Not dumping the config file yet.\n";
	
	if (variablesMap.count("help"))
	{
		cout << visibleOptions << "\n";
		return variablesMap;
	}
	
	if (variablesMap.count("version"))
	{
        cout << "Trogdor version: " << TROGDOR_VERSION_TEXT << endl;
		cout << "Compile date: " << __DATE__ << endl;
		//cout << "OS type: " << TROGDOR_OS << endl;
		cout << "Boost version: " << BOOST_LIB_VERSION << endl;
//		cout << "TinyXML version: " << TIXML_MAJOR_VERSION << "."
//			<< TIXML_MINOR_VERSION << "." << TIXML_PATCH_VERSION << endl;
		//cout << "vmlib version: 5.1" << endl;
		cout << "calc version: 4.7" << endl;
		
		return variablesMap;
	}
	
	if (variablesMap.count("numerics"))
	{
        vector<int*> v;
        cout << "sizeof(int) = " << sizeof(int) << endl;
        cout << "sizeof(long) = " << sizeof(long) << endl;
        cout << "sizeof(long long) = " << sizeof(long long) << endl;
        cout << "sizeof(int*) = " << sizeof(int*) << endl;
        cout << "sizeof(size_t) = " << sizeof(size_t) << endl;
        cout << "vector<int*>::max_size() = " << v.max_size() << endl;
		typedef numeric_limits<Precision::Float> f;
		cout << "Values from std::limits:\n";
		cout << "digits = " << f::digits << endl;
		cout << "digits10 = " << f::digits10 << endl;
		cout << "epsilon = " << f::epsilon() << endl;
		cout << "min = " << f::min() << endl;
		cout << "min_exponent = " << f::min_exponent << endl;
		cout << "min_exponent10 = " << f::min_exponent10 << endl;
		cout << "max = " << f::max() << endl;
		cout << "max_exponent = " << f::max_exponent << endl;
		cout << "max_exponent10 = " << f::max_exponent10 << endl;
		if (f::has_denorm == denorm_present)
		{
			cout << "has_denorm = denorm_present\n";
			cout << " denorm_min = " << f::denorm_min() << endl;
			cout << " has_denorm_loss = " << f::has_denorm_loss << endl;
		}
		else if (f::has_denorm == denorm_absent)
			cout << "has_denorm = denorm_absent\n";
		else if (f::has_denorm == denorm_indeterminate)
			cout << "has_denorm = denorm_indeterminate\n";
		cout << "has_infinity = " << f::has_infinity << endl;
		cout << "has_quiet_NaN = " << f::has_quiet_NaN << endl;
		cout << "has_signaling_NaN = " << f::has_signaling_NaN << endl;
		cout << "round_error = " << f::round_error() << endl;
		cout << "round_style = ";
		switch (f::round_style)
		{
			case round_indeterminate:
				cout << "round_indeterminate\n";
				break;
			case round_toward_zero:
				cout << "round_toward_zero\n";
				break;
			case round_to_nearest:
				cout << "round_to_nearest\n";
				break;
			case round_toward_infinity:
				cout << "round_toward_infinity\n";
				break;
			case round_toward_neg_infinity:
				cout << "round_toward_neg_infinity\n";
				break;
			default:
				cout << "unknown\n";
				break;
		}
		cout << "tinyness_before = " << f::tinyness_before << endl;
		cout << "traps = " << f::traps << endl;
	}
	
	return variablesMap;
}


