#include "xhmmDriver.hpp"
#include "xhmmCmdline.h"

#include <Exception.hpp>

#include <sys/times.h>
#include <unistd.h>

#include <cstdlib>
#include <string>
#include <iostream>
using namespace std;

#define HEADER_SEPARATOR "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"

int main(int argc, char** argv) {
	struct gengetopt_args_info args_info;
	if (cmdline_parser(argc, argv, &args_info) != 0) {
		exit(1);
	}

	if (argc <= 1) {
		cmdline_parser_print_help();
		exit(1);
	}

	try {
		cerr << "command executed:";
		for (int i = 0; i < argc; ++i)
			cerr << " " << argv[i];
		cerr << endl << endl;

		cerr << HEADER_SEPARATOR << endl;
		cerr << "All parameters and default values:" << endl;
		cerr << HEADER_SEPARATOR << endl;
		int i = 0;
		while (gengetopt_args_info_help[i]) {
			string optStr = gengetopt_args_info_help[i++];
			if (optStr.find("--help") == string::npos && optStr.find("--detailed-help") == string::npos && optStr.find("--full-help") == string::npos && optStr.find("--version") == string::npos)
				cerr << optStr << endl;
		}
		cerr << endl;

		cerr << HEADER_SEPARATOR << endl;
		cerr << "Command-line parameter values:" << endl;
		cerr << HEADER_SEPARATOR << endl;
		cmdline_parser_dump(stderr, &args_info);
		cerr << HEADER_SEPARATOR << endl;
		cerr << endl;

		//Start times:
		time_t startTime = time(NULL);

		XHMM::xhmmDriver* xhmm = new XHMM::xhmmDriver(args_info);
		bool printRunTime = xhmm->run();
		delete xhmm;

		if (printRunTime) {
			//End times:
			time_t endTime = time(NULL);
			struct tms finishCpuTimeStats;
			times(&finishCpuTimeStats);

			clock_t cpuTime =
					(finishCpuTimeStats.tms_utime + finishCpuTimeStats.tms_stime +
							finishCpuTimeStats.tms_cutime + finishCpuTimeStats.tms_cstime);

			clock_t ownUserTime = finishCpuTimeStats.tms_utime;
			clock_t ownSystemTime = finishCpuTimeStats.tms_stime;

			clock_t childrenUserTime = finishCpuTimeStats.tms_cutime;
			clock_t childrenSystemTime = finishCpuTimeStats.tms_cstime;

			const long CLK_TCKS_PER_SECOND = sysconf(_SC_CLK_TCK);
			cpuTime /= CLK_TCKS_PER_SECOND;

			ownUserTime /= CLK_TCKS_PER_SECOND;
			ownSystemTime /= CLK_TCKS_PER_SECOND;
			childrenUserTime /= CLK_TCKS_PER_SECOND;
			childrenSystemTime /= CLK_TCKS_PER_SECOND;

			//Print times:
			cerr << endl
					<< "Total CPU time for processing this job: " << cpuTime << " seconds" << endl;
			cerr << "[user time: " << ownUserTime << endl
					<< "system time: " << ownSystemTime << endl
					<< "children user time: " << childrenUserTime << endl
					<< "children system time: " << childrenSystemTime;
			cerr << "]" << endl;
			cerr << "Total time for processing this job: " << endTime-startTime << " seconds" << endl;
		}

		exit(0);
	}
	catch (Exception* e) {
		e->printErrorMessage();
		delete e;
		exit(1);
	}
}
