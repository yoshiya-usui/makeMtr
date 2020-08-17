//--------------------------------------------------------------------------
// This file is part of makeMtr.
//
// makeMtr is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// makeMtr is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with makeMtr. If not, see <http://www.gnu.org/licenses/>.
//--------------------------------------------------------------------------
#ifndef DBLDEF_OBSERVATION_LINE
#define DBLDEF_OBSERVATION_LINE

#include <fstream>
#include <iostream>
#include <map>
#include "CommonParameters.h"
#include "ObservationPoint.h"

// Class of observation line
class ObservationLine{

public:
	// Default Constructer
	ObservationLine();

	// Destructer
	~ObservationLine();

	// Calculate maximum length of specified coordinate on X-Y plane ( z = 0 ) 
	double calcMaximumLengthOfPoint( const CommonParameters::XYZ& coord ) const;

	// Read data of observation points consisting line from input file
	void readObservationLineData( std::ifstream& ifs );

private:
	// Copy constructer
	ObservationLine(const ObservationLine& rhs);

	// Assignment operator
	ObservationLine& operator=(ObservationLine& rhs);

	// Pair of the points
	CommonParameters::XYZ m_endPoints[2];

	// Total number of layers
	int m_numLayers;

	// Radius for local mesh control
	double* m_radius;

	// Maximum edge length within layer around the line
	double* m_maxEdgeLengthWithinLayer;

};

#endif
