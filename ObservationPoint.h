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
#ifndef DBLDEF_OBSERVATION_POINT
#define DBLDEF_OBSERVATION_POINT

#include <fstream>
#include <iostream>
#include "CommonParameters.h"

// Class of observation point
class ObservationPoint{

public:
	// Default Constructer
	ObservationPoint();

	// Destructer
	~ObservationPoint();

	// Read data of observation point from input file
	void readObservationPointData( std::ifstream& ifs );

	// Calculate maximum length of specified coordinate
	double calcMaximumLengthOfPoint( const CommonParameters::XYZ& coord ) const;

private:
	// Copy constructer
	ObservationPoint(const ObservationPoint& rhs);

	// Assignment operator
	ObservationPoint& operator=(ObservationPoint& rhs);

	// Coordinate of the point
	CommonParameters::XYZ m_pointCoord;

	// Total number of spheres
	int m_numSpheres;

	// Radius for local mesh control
	double* m_radius;

	// Maximum edge length within sphere around the point
	double* m_maxEdgeLengthWithinSphere;

	// Oblateness of spheres
	double* m_oblateness;

	bool locateInSphere( const CommonParameters::XYZ& coord, const int iSphere ) const;

};

#endif
