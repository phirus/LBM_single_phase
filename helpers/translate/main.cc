#include <iostream>
#include<ctime>

#include"../../src/BasicIO.h"
#include<boost/program_options.hpp>

using namespace std;
namespace po = boost::program_options;

int main(int argc, char** argv){

    boost::program_options::options_description desc("Allowed options");
	desc.add_options()
        ("help,h", "produce help message")
        ("input,i", boost::program_options::value<string> (), "specify preprocess parameter file")
        ;
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::parse_command_line(argc,argv,desc),vm);
    boost::program_options::notify(vm);

    if(vm.count("help")){
        cout << desc << endl;
        return 1;
    }

    Preprocess prepro = read_preprocess_file("preprocessFile");
    Timetrack timetrack = read_timetrack_file("preprocessFile");

    if (vm.count("input")) {
        cout << "preprocess file is: " << vm["input"].as<string>() << ".\n" << endl ;
        prepro = read_preprocess_file(vm["input"].as<string>());
        timetrack = read_timetrack_file(vm["input"].as<string>());
    }
   
     const ParamSet params = prepro.getParamSet();
     write_param_log(params);
     write_param_log_csv(params);       

    return 0;
}
