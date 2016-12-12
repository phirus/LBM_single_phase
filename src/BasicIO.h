#ifndef BASICIO_H
#define BASICIO_H

#include<fstream>
#include<iostream>
#include<sstream>
#include<map>
#include<string.h>
#include<vector>
#include<limits>

#include"Preprocess.h"
#include"Timetrack.h"
#include"Boundaries.h"

using namespace std;
typedef std::vector<std::vector<double>> nested_vector; 
/// write output
void write_file_header(const string& filename = "BubblePlot.csv", const string& header= "time , Posx , PosY , v_x , v_y , Re");
void write_data_plot(const std::vector<double> y, double del_x, const string& filename = "ReynoldsPlot.dat");
void write_data_plot_linewise(int time ,double y1, double y2, const string& filename = "BubbleVeloPlot.dat");
void write_csv(const nested_vector& data, const string& filename = "BubbleVeloPlot.csv", const string& header= "time , Posx , PosY , v_x , v_y , Re");
void write_csv_linewise(int i, double Posx, double PosY, double v_x, double v_y, double Re, const string& filename = "BubbleVeloPlot.csv");
void write_param_log(const ParamSet& p);
void write_param_log_csv(const ParamSet& p);
void write_preprocess_csv(const Preprocess& p);

/// read input
const bool input_query(const string& filename, const string& query, double& value);
const map<string,double> assign_map_via_file(map<string,double> mm, const string& filename);
const Preprocess read_preprocess_file(const string& filename);
const Timetrack read_timetrack_file(const string& filename);
const ParamSet read_paramset_file(const string& filename = "paramLog");
const Boundaries read_boundaries_file(const string& filename);

///auxiliary
const string createFilename(const string& name, int iteration, const string& type = ".bin");

#endif