#include"BasicIO.h"


//=========================== WRITE OUTPUT ===========================

void write_file_header(const string& filename, const string& header)
{
    ofstream Plot;
    Plot.open( filename.c_str() );
    Plot << header << "\n";
    Plot.close();
}

void write_data_plot(const std::vector<double> y, double del_x, const string& filename)
{
    ofstream RePlot;

    RePlot.open( filename.c_str() );
    RePlot << "x" << "\t" << "y" << "\n";
    for(unsigned int i = 0; i< y.size(); i++){
        RePlot << i*del_x << "\t" << y.at(i) << "\n";
    } 
    RePlot.close();
}

void write_data_plot_linewise(int time ,double y1, double y2, const string& filename)
{
    ofstream Plot;
    Plot.open( filename.c_str(), ios::app);
    Plot << time << "\t" << y1 << "\t" << y2 << "\n"; 
    Plot.close();
}

void write_csv(const nested_vector& data, const string& filename, const string& header)
{
    const auto numOfSets = data.size();
    const auto numOfLines = data.at(0).size();

    ofstream Plot;
    Plot.open( filename.c_str() );
    Plot << header << "\n";

    for (unsigned int line = 0; line  < numOfLines ; line++)
    {
        for (unsigned int item = 0; item < numOfSets; item++)
        {
            const string end_char = (item == numOfSets -1) ? "\n" : ";" ;
            Plot << data[item][line] << end_char;
        }
    } 
    Plot.close();
}

void write_csv_linewise(int i, double Posx, double PosY, double v_x, double v_y, double Re, const string& filename)
{
    ofstream Plot;
    Plot.open(filename.c_str(),ios::app);
    Plot << i << ";" << Posx << ";" << PosY << ";" << v_y << ";" << v_x << ";" << Re << "\n";
    Plot.close();
}

void write_param_log(const ParamSet& p)
{
    ofstream paramLog;

    stringstream name;
    name <<"paramLog";

    paramLog.open( name.str().c_str() );
    paramLog << "# used setup parameters" << endl;
    paramLog << "\n# omega [-]" << endl;
    paramLog << "omega = "  << p.getOmega()           << endl;
    paramLog << "\n# density, normalized [-]" << endl;
    paramLog << "rho = "    << p.getRho()               << endl;
    paramLog << "\n# discretisization [-]" << endl;
    paramLog << "dt = "         << p.getDeltaT()             << endl;
    paramLog << "dx = "         << p.getDeltaX()             << endl;
    
    RelaxationPar2D relax = p.getRelaxation2D();
    paramLog << "\n# MRT parameters [-]" << endl;
    paramLog << "s_2 = " << relax.s_2 << endl;
    paramLog << "s_3 = " << relax.s_3 << endl;
    paramLog << "s_5 = " << relax.s_5 << endl;

    paramLog.close();
}

void write_param_log_csv(const ParamSet& p)
{
    ofstream paramLog;

    stringstream name;
    name <<"paramLog.csv";
    //paramLog.precision(std::numeric_limits< double >::max_digits10);
    paramLog.open( name.str().c_str() );

    paramLog << "omega;rho;dt;dx;s_2;s_3;s_5" << endl;
    RelaxationPar2D relax = p.getRelaxation2D();
    paramLog <<  p.getOmega() << ";" << p.getRho() << ";" ;
    paramLog << p.getDeltaT() << ";" << p.getDeltaX() << ";" ;
    paramLog << relax.s_2 << ";" << relax.s_3 << ";" << relax.s_5 << endl;

    paramLog.close();
}

void write_preprocess_csv(const Preprocess& p)
{
    ofstream preproCSV;

    stringstream name;
    name <<"prepro.csv";
    preproCSV.open( name.str().c_str() );

    preproCSV << "ReynoldsMax;Morton;Eotvos;resolution;muRatio;isShearFlow;shearRate;xCells;yCells;zCells"<<endl;
    preproCSV << p.getReynoldsMax() << ";" << p.getResolution() << ";" ;
    preproCSV << p.getMuRatio() << ";" << p.getIsShearFlow() << ";" << p.getShearRate() << ";" ;
    preproCSV << p.getXCells() << ";" << p.getYCells() << ";" << p.getZCells() << endl;
    preproCSV.close();
}

//=========================== READ INPUT ===========================

const bool input_query(const string& filename, const string& query, double& value)
{
    bool success = false;
    value = 0;
    string lineString;
    string word;
    fstream file(filename.c_str(),ios::in);

    while(file.good()){
    getline(file,lineString);
    if ( lineString.substr(0,1) != "#" ){
        istringstream lineStream(lineString);
        lineStream >> word;
        if (word == query) {
            lineStream >> word; // gets the "="
            lineStream >> value;
            success = true;
            break;
            }         
        }
    }
    file.close();
    return success;
}

const map<string,double> assign_map_via_file(map<string,double> mm, const string& filename)
{
    double tmp;
    for(map<string,double>::iterator it = mm.begin(); it != mm.end(); it++){
        if( input_query(filename,it->first,tmp) == true ) it->second = tmp;
    }
    return mm;
}

const Preprocess read_preprocess_file(const string& filename)
{

    // initialzing strings and fallback values
    map<string,double> mm;
    mm.insert(pair<string,double>("Reynolds",10));
    mm.insert(pair<string,double>("resolution",30));
    mm.insert(pair<string,double>("rho",1));
    mm.insert(pair<string,double>("mu_ratio",2));
    mm.insert(pair<string,double>("s_3",1));
    mm.insert(pair<string,double>("s_5",1));
    mm.insert(pair<string,double>("isShearFlow",0)); // implicit type cast: 0 -> false, 1 -> true
    mm.insert(pair<string,double>("shearRate",0));
    mm.insert(pair<string,double>("xCells",50));
    mm.insert(pair<string,double>("yCells",50));
    mm.insert(pair<string,double>("zCells",50));

    mm = assign_map_via_file(mm, filename);

    Preprocess prepro(mm.at("Reynolds"),mm.at("resolution"),mm.at("rho"), mm.at("mu_ratio"), mm.at("s_3"), mm.at("s_5"), mm.at("isShearFlow"), mm.at("shearRate"), mm.at("xCells"), mm.at("yCells"), mm.at("zCells"));
    return prepro;
}

const Timetrack read_timetrack_file(const string& filename)
{
    // initialzing strings and fallback values
    map<string,double> mm;
    mm.insert(pair<string,double>("max_steps",1e5));
    mm.insert(pair<string,double>("output_interval",2e3));
    mm.insert(pair<string,double>("restart_interval",1e4));
    
    mm = assign_map_via_file(mm, filename);

    Timetrack time(mm.at("max_steps"), mm.at("output_interval"), mm.at("restart_interval"));
 
    return time;
}

const ParamSet read_paramset_file(const string& filename)
{
    // initialzing strings and fallback values
    map<string,double> mm;
    mm.insert(pair<string,double>("omega",1));
    mm.insert(pair<string,double>("rho",1));
    mm.insert(pair<string,double>("dt",1e-3));
    mm.insert(pair<string,double>("dx",1e-3));
    mm.insert(pair<string,double>("s_2",1));
    mm.insert(pair<string,double>("s_3",1));
    mm.insert(pair<string,double>("s_5",1));
         
    mm = assign_map_via_file(mm, filename);

    const RelaxationPar2D rel(mm.at("s_2"),mm.at("s_3"),mm.at("s_5"));

    ParamSet params(mm.at("omega"), mm.at("rho"), mm.at("dt"), mm.at("dx"), rel); /// < consructor
    return params;
}

const Boundaries read_boundaries_file(const string& filename)
{   
    // initialzing strings and fallback values
    map<string,double> mm;
    mm.insert(pair<string,double>("n_type",0));
    mm.insert(pair<string,double>("n_rho",0));
    mm.insert(pair<string,double>("nvx",0));
    mm.insert(pair<string,double>("nvy",0));

    mm.insert(pair<string,double>("s_type",0));
    mm.insert(pair<string,double>("s_rho",0));
    mm.insert(pair<string,double>("svx",0));
    mm.insert(pair<string,double>("svy",0));

    mm.insert(pair<string,double>("e_type",0));
    mm.insert(pair<string,double>("e_rho",0));
    mm.insert(pair<string,double>("evx",0));
    mm.insert(pair<string,double>("evy",0));

    mm.insert(pair<string,double>("w_type",0));
    mm.insert(pair<string,double>("w_rho",0));
    mm.insert(pair<string,double>("wvx",0));
    mm.insert(pair<string,double>("wvy",0));

    mm = assign_map_via_file(mm, filename);

    double tmp_rho;
    Vector2D tmp_u;

    Boundaries b;

    tmp_rho = mm.at("n_rho");
    tmp_u = Vector2D(mm.at("nvx"),mm.at("nvy"));
    b.north = BoundaryInformation(mm.at("n_type"),tmp_rho,tmp_u);

    tmp_rho = mm.at("s_rho");
    tmp_u = Vector2D(mm.at("svx"),mm.at("svy"));
    b.south = BoundaryInformation(mm.at("s_type"),tmp_rho,tmp_u);

    tmp_rho = mm.at("e_rho");
    tmp_u = Vector2D(mm.at("evx"),mm.at("evy"));
    b.east = BoundaryInformation(mm.at("e_type"),tmp_rho,tmp_u);

    tmp_rho = mm.at("w_rho");
    tmp_u = Vector2D(mm.at("wvx"),mm.at("wvy"));
    b.west = BoundaryInformation(mm.at("w_type"),tmp_rho,tmp_u);

    return b;
}

//=========================== AUXILIARY ===========================

const string createFilename(const string& name, int iteration, const string& type)
{
    stringstream filename;
    filename << name << iteration << type;
    return filename.str();
}