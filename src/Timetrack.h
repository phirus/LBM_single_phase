/// tracks the physical time and termination conditions throughout the simulation 

#ifndef TIMETRACK_H
#define TIMETRACK_H

#include<vector>
#include<cmath>

using namespace std;

class Timetrack
{
public: 
	/// Lifecycle
	Timetrack(int t_c = 1e5, int out = 2e3, int restart = 1e4); //, double resi = 1e-3);
	
	/// operations	
	inline void timestep(){count++;};
	
	inline const bool proceed()const{return (count <= terminalCount);};

	/// accessors
	inline const int getCount()const{return count;};
	inline const int getMaxCount()const{return terminalCount;};
	inline const int getOutputInt()const{return outputInterval;};
	inline const int getRestartInt()const{return restartInterval;};
	
	inline void setCount(int c){count = c;};
	inline void setMaxCount(int t_c){terminalCount = t_c;};
	inline void setOutputInt(int tmp){outputInterval = tmp;};
	inline void setRestartInt(int tmp){restartInterval = tmp;};
	
	/// operators
	const bool operator==(const Timetrack& other)const;

private:
	int count;
	
	// termination conditions
	int terminalCount;			// maximum number of time steps
	
	// output intervals
	int outputInterval;
	int restartInterval;
};
#endif