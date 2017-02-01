#include<iostream>
#include<ctime>
#include<vector>
#include<chrono>
#include<cmath>

#include"../src/Lattice.h"
#include"../src/BinaryIO2D.h"
#include<boost/program_options.hpp>

using namespace std;
namespace po = boost::program_options;

void initialSetUp(Lattice& meins, Preprocess& prepro, Boundaries& bound, int xmax, int ymax, ParamSet params, int numOfCPUs);

const double normFun(double x){return 0.5 * (sin( (x * 0.001 / PI) - PI/2 ) + 1);}; 
const BoundaryInformation setBoundVelo(const BoundaryInformation& bound, int i);

int main(int argc, char** argv){

    boost::program_options::options_description desc("Allowed options");
	desc.add_options()
        ("help,h", "produce help message")
        ("boundary,b", boost::program_options::value<string> (), "specify boundary input file")
        ("cpu,c", boost::program_options::value<int> (), "takes the number of CPUs")
        ("preprocess,p", boost::program_options::value<string> (), "specify preprocess parameter file")
        ("override,o", boost::program_options::value<string> (), "specify parameter file to bypass preprocess routine")
        ("restart,r", boost::program_options::value<string> (), "specify restart file")
        ;
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::parse_command_line(argc,argv,desc),vm);
    boost::program_options::notify(vm);

    if(vm.count("help"))
    {
        cout << desc << endl;
        return 1;
    }

    int numOfCPUs = 1;
    if(vm.count("cpu"))
    {
        numOfCPUs = vm["cpu"].as<int>();
        cout << "number of CPUs set to "<< numOfCPUs << endl;
    }

    Preprocess prepro = read_preprocess_file("preprocessFile");
    Timetrack timetrack = read_timetrack_file("preprocessFile");
    ParamSet params = prepro.getParamSet();
    Boundaries boundaries = read_boundaries_file("BoundaryInput");

    Preprocess prepro_input;
    Timetrack timetrack_input;    
    ParamSet params_input;

    if (vm.count("preprocess")) 
    {
        cout << "preprocess file is: " << vm["preprocess"].as<string>() << ".\n" << endl ;
        prepro_input = read_preprocess_file(vm["preprocess"].as<string>());
        timetrack_input = read_timetrack_file(vm["preprocess"].as<string>());
        params_input = prepro_input.getParamSet();
        
        prepro = prepro_input;
        timetrack = timetrack_input;
        params = params_input;
    }

        if (vm.count("boundary")) 
    {
        cout << "boundary file is: " << vm["boundary"].as<string>() << ".\n" << endl ;
        boundaries = read_boundaries_file(vm["boundary"].as<string>());
    }

    if (vm.count("override")) 
    {
        cout << "override preprocess file with: " << vm["override"].as<string>() << ".\n" << endl ;
        params_input = read_paramset_file(vm["override"].as<string>());
        params = params_input;
    }

    // create a Lattice   
    int ymax = prepro.getYCells();
    int xmax = prepro.getXCells();
    Lattice meins(xmax,ymax);

    if (vm.count("restart")) 
    {
        cout << "Restart file is: " << vm["restart"].as<string>() << ".\n" << endl ;
        Lattice tmpL;
        Preprocess tmpP;
        Timetrack tmpT;
        bool tmpB = read_restart_file2D(tmpL, tmpP, tmpT, vm["restart"].as<string>());
        if (tmpB == true)
        {
            meins = tmpL;
            prepro = tmpP;
            timetrack = tmpT;
        }
        else 
        {
            cout << "failed to read input file" << endl;
            return 1;
        }
        if (vm.count("preprocess"))
        {
            timetrack.setMaxCount( timetrack_input.getMaxCount() );
            timetrack.setOutputInt( timetrack_input.getOutputInt() );
            timetrack.setRestartInt( timetrack_input.getRestartInt() );
        }
    }
    else {      // vm.count("restart"), if no restart file -> initialize
       
       initialSetUp(meins, prepro, boundaries, xmax, ymax, params,numOfCPUs);
       write_vtk_output2D(meins, 0);

    }

    //time_t start,end;
    std::chrono::high_resolution_clock::time_point start,end;
    auto parallel = std::chrono::duration_cast<std::chrono::microseconds>( end - end).count();
    auto sequential = std::chrono::duration_cast<std::chrono::microseconds>( end - end).count();

//    time(&start);

    const int outputInterval = timetrack.getOutputInt();
    const int restartInterval = timetrack.getRestartInt();

    int i = 0;

    const Boundaries final_bound = meins.getBoundaries();
    Boundaries tmp_bound = meins.getBoundaries();

    tmp_bound.south = setBoundVelo(final_bound.south, i);
    meins.setBoundaries(tmp_bound);



    while (timetrack.proceed() == true)
    {
        start = std::chrono::high_resolution_clock::now();

        meins.collideAll(numOfCPUs);
        meins.evaluateBoundaries(numOfCPUs);
        meins.streamAll(numOfCPUs);

        end = std::chrono::high_resolution_clock::now();
        parallel +=  std::chrono::duration_cast<std::chrono::microseconds>( end - start).count();

        start = std::chrono::high_resolution_clock::now();
        timetrack.timestep();

        // Output if necessary
        i = timetrack.getCount();

        if(i <= 10000)
        {
            tmp_bound.south = setBoundVelo(final_bound.south, i);
            meins.setBoundaries(tmp_bound);
        }
        
        if(i%100 == 0) 
        {
            cout << i <<endl;
        }

        if(i%outputInterval == 0) 
        {
            write_vtk_output2D(meins, i);
        } 
        
        if(i%restartInterval == 0)
        {
            const string restart_file_name =  createFilename("restart", i, ".bin");
            write_restart_file2D(meins, prepro, timetrack, restart_file_name);
        }       
    }

//    time(&end);
//    cout<<"\nBerechnung beendet nach "<< difftime(end,start) <<" Sekunden"<<endl;
      cout<<"\nBerechnung beendet nach "<< sequential+parallel <<" Micro-Sekunden"<<endl;
      cout<<"\nDavon parallel:  "<< parallel <<" Micro-Sekunden"<<endl;
      cout<<"\nDavon sequentiell "<< sequential <<" Micro-Sekunden"<<endl;

    return 0;
}


void initialSetUp(Lattice& meins, Preprocess& prepro, Boundaries& bound, int xmax, int ymax, ParamSet params, int numOfCPUs)
{
    // set the parameters        
    meins.setParams(params);

    // get densities
    const Cell fluid(1,false);
    const Cell wall(0,true);

    // setup geometry (bubble at the bottom, x-centered)
    const int R1 = prepro.getResolution()/2;
    const int xm1 = xmax / 2;
    const int ym1 = R1 + xm1;

   for(int j=0; j< ymax; j++)
   {
       for(int i=0; i< xmax; i++){
           if( (i-xm1)*(i-xm1) + (j-ym1)*(j-ym1) < R1*R1 ) meins.setCell(i,j,wall);
           else meins.setCell(i,j,fluid);
       }
   }

    meins.setBoundaries(bound);
    meins.equilibriumIni();

    cout<<"Initialisierung beendet\n\n"<<endl;
}

const BoundaryInformation setBoundVelo(const BoundaryInformation& bound, int i)
{
    Vector2D v = bound.getVelocity();
    v = v * normFun(i);
    BoundaryInformation tmp_bound(bound.getType(),bound.getRho(),v);

    return tmp_bound;
}