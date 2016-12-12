#include"BinaryIO2D.h"

//=========================== BINARY DUMP ===========================

void write_binary2D(const Lattice& l, const string& filename){

    // setting up the file name
    stringstream name;
    name << filename;

    // setting up file
    fstream file(name.str().c_str(),ios::out | ios::binary);
    file.seekp(0);

    // start to write
    DimSet2D extent = l.getSize();
    file.write(reinterpret_cast<char*> (&extent), sizeof extent);

    ParamSet param = l.getParams();
    file.write(reinterpret_cast<char*> (&param), sizeof param);

    Boundaries bound = l.getBoundaries();
    file.write(reinterpret_cast<char*> (&bound), sizeof bound);

    field2D data = l.getData();

    for (int x = 0; x<extent[0];x++){
        for (int y = 0; y<extent[1];y++){
            file.write(reinterpret_cast<char*> (&data[x][y]), sizeof(Cell));
        }
    }
    file.close();
}

const bool read_binary2D(Lattice& outL, const string& filename){
    bool success;

    // setting up file
    fstream file(filename.c_str(),ios::in | ios::binary);
    if(file.is_open()){
        success = true;
        file.seekg(0);

        // start to read
        DimSet2D extent;
        file.read((char*) &extent, sizeof extent);

        ParamSet param;
        file.read((char*) &param, sizeof param);

        Boundaries bound;
        file.read((char*) &bound, sizeof bound);

        Cell tmpCell;
        field2D data(boost::extents[extent[0]][extent[1]]);
        for(int x = 0; x<extent[0];x++){
            for(int y=0;y<extent[1];y++){
                file.read((char*) &tmpCell, sizeof(Cell));
                data[x][y] = tmpCell;
            }
        }
        file.close();
        outL.setParams(param);
        outL.setBoundaries(bound);
        outL.setData(data, extent[0], extent[1]);
    }
    else success = false;

    return success;
}

//=========================== RESTART FILES ===========================

void write_restart_file2D(const Lattice& l, const Preprocess& p, const Timetrack time, const string& filename){

    // setting up the file name
    stringstream name;
    name << filename;

    // setting up file
    fstream file(name.str().c_str(),ios::out | ios::binary);
    file.seekp(0);

    // start to write
    DimSet2D extent = l.getSize();
    file.write(reinterpret_cast<char*> (&extent), sizeof extent);

    ParamSet param = l.getParams();
    file.write(reinterpret_cast<char*> (&param), sizeof param);

    Boundaries bound = l.getBoundaries();
    file.write(reinterpret_cast<char*> (&bound), sizeof bound);


    // write the velocity distributions
    field2D data = l.getData();
    for (int x = 0; x<extent[0];x++){
        for (int y = 0; y<extent[1];y++){
            file.write(reinterpret_cast<char*> (&data[x][y]), sizeof(Cell));
        }
    }

    // write Timetrack
    int count = time.getCount();
    int maxCount = time.getMaxCount();
    int output_interval = time.getOutputInt();
    int restart_interval = time.getRestartInt();

    file.write(reinterpret_cast<char*> (&count), sizeof count);
    file.write(reinterpret_cast<char*> (&maxCount), sizeof maxCount);
    file.write(reinterpret_cast<char*> (&output_interval), sizeof output_interval);
    file.write(reinterpret_cast<char*> (&restart_interval), sizeof restart_interval);

    // write Preprocess

    double ReynoldsMax = p.getReynoldsMax();    
    double resolution = p.getResolution();
    double rho = p.getRho();
    double mu_ratio = p.getMuRatio();
    double s_3 = p.getS_3();
    double s_5 = p.getS_5();
    bool isShearFlow = p.getIsShearFlow();
    double shearRate = p.getShearRate();
    int xCells = p.getXCells();
    int yCells = p.getYCells();
    int zCells = p.getZCells();

    file.write(reinterpret_cast<char*> (&ReynoldsMax), sizeof(double));
    file.write(reinterpret_cast<char*> (&resolution), sizeof(double));
    file.write(reinterpret_cast<char*> (&rho), sizeof(double));
    file.write(reinterpret_cast<char*> (&mu_ratio), sizeof(double));    
    file.write(reinterpret_cast<char*> (&s_3), sizeof(double));
    file.write(reinterpret_cast<char*> (&s_5), sizeof(double));
    file.write(reinterpret_cast<char*> (&isShearFlow), sizeof(bool));
    file.write(reinterpret_cast<char*> (&shearRate), sizeof(double));
    file.write(reinterpret_cast<char*> (&xCells), sizeof(int));
    file.write(reinterpret_cast<char*> (&yCells), sizeof(int));
    file.write(reinterpret_cast<char*> (&zCells), sizeof(int));

    file.close();
}

const bool read_restart_file2D(Lattice& outL, Preprocess& p, Timetrack& t, const string& filename)
{
    bool success;

    // setting up file
    fstream file(filename.c_str(),ios::in | ios::binary);
    if(file.is_open()){
        success = true;
        file.seekg(0);

        // start to read
        DimSet2D extent;
        file.read((char*) &extent, sizeof extent);

        ParamSet param;
        file.read((char*) &param, sizeof param);

        Boundaries bound;
        file.read((char*) &bound, sizeof bound);

        Cell tmpCell;
        field2D data(boost::extents[extent[0]][extent[1]]);
        for(int x = 0; x<extent[0];x++){
            for(int y=0;y<extent[1];y++){
                file.read((char*) &tmpCell, sizeof(Cell));
                data[x][y] = tmpCell;
            }
        }

        int count;
        int maxCount;
        int output_interval;
        int restart_interval;

        file.read((char*) &count, sizeof count);
        file.read((char*) &maxCount, sizeof maxCount);
        file.read((char*) &output_interval, sizeof output_interval);
        file.read((char*) &restart_interval, sizeof restart_interval);

        Timetrack time(maxCount, output_interval, restart_interval);
        time.setCount(count);

        double ReynoldsMax;
        double resolution, rho;
        double mu_ratio, s_3, s_5; 
        bool isShearFlow;
        double shearRate;
        int xCells, yCells, zCells;

        file.read((char*) &ReynoldsMax, sizeof(double));
        file.read((char*) &resolution, sizeof(double));
        file.read((char*) &rho, sizeof(double));
        file.read((char*) &mu_ratio, sizeof(double));
        file.read((char*) &s_3, sizeof(double));
        file.read((char*) &s_5, sizeof(double));
        file.read((char*) &isShearFlow, sizeof(bool));
        file.read((char*) &shearRate, sizeof(double));
        file.read((char*) &xCells, sizeof(int));
        file.read((char*) &yCells, sizeof(int));
        file.read((char*) &zCells, sizeof(int));

        Preprocess prepro(ReynoldsMax, resolution, rho, mu_ratio, s_3, s_5, isShearFlow, shearRate, xCells, yCells, zCells);
        
        file.close();
        outL.setParams(param);
        outL.setBoundaries(bound);
        outL.setData(data, extent[0], extent[1]);
        t = time;
        p = prepro;
    }
    else success = false;

    return success;
}

//=========================== WRITE OUTPUT ===========================

void write_vtk_output2D(const Lattice& l, const string& filename)
{
    ofstream VTKFile;
    Cell tmp;
    int e;
    DimSet2D extent = l.getSize();
    int xsize = extent[0];
    int ysize = extent[1];

    VTKFile.open(filename.c_str());

    VTKFile << "# vtk DataFile Version 3.1" << endl;
    VTKFile << "Lattice Boltzmann data" << endl;
    VTKFile << "ASCII" << endl;
    VTKFile << "DATASET UNSTRUCTURED_GRID" << endl;

    VTKFile << "POINTS "<< (xsize+1) * (ysize+1)  <<" INT \n";

    for (int j = 0; j <= ysize; j++)
    {
        for (int i = 0; i <= xsize; i++)
        {
            VTKFile << i << " " << j << " 0 " ;
        }
        VTKFile<<endl;
    }

    VTKFile << "\nCELLS " << (xsize) * (ysize) << " " << (xsize) * (ysize) * 5 << "\n";
    
    for (int j = 0; j < ysize; j++)
    {
        for (int i = 0; i < xsize; i++)
        {
            e = i+(xsize+1)*j;            
            VTKFile <<"4 "<< e << " " << e+1 << " "<< e + xsize +1 << " " << e + xsize + 2 << " ";
        }
        VTKFile<< endl;
    }
    VTKFile << "\nCELL_TYPES "<< (xsize) * (ysize) << "\n";
    for (int q = 0; q < (xsize * ysize); q++)
    {
        VTKFile <<"8 ";
    }

    VTKFile << "\nCELL_DATA "<< (xsize) * (ysize) << endl;

    VTKFile << "\nSCALARS Rho DOUBLE\nLOOKUP_TABLE default"<<endl;
    for (int j = 0; j < ysize; j++)
    {
        for (int i = 0; i < xsize; i++)
        {
            tmp = l.getCell(i,j);
            tmp.calcRho();
            VTKFile << tmp.getRho() << " ";
        }
    }

    VTKFile << "\nSCALARS U_abs DOUBLE\nLOOKUP_TABLE default"<<endl;
    for (int j = 0; j < ysize; j++)
    {
        for (int i = 0; i < xsize; i++)
        {
            tmp = l.getCell(i,j);
            tmp.calcRho();
            Vector2D u = tmp.getU();
            VTKFile << u.Abs() << " ";
        }
    }



    VTKFile << "\nVECTORS u DOUBLE"<<endl;
    for (int j = 0; j < ysize; j++)
    {        
        for (int i = 0; i < xsize; i++)
        {
            tmp = l.getCell(i,j);
            tmp.calcRho();

            Vector2D u = tmp.getU();
            VTKFile << u.x << " " << u.y << " 0 ";
        }
    }




    VTKFile.close();
}