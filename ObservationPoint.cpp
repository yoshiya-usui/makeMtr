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
#include "ObservationPoint.h"
#include "math.h"
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <assert.h> 

// Default constructer
ObservationPoint::ObservationPoint():
	m_numSpheres(0),
	m_radius(NULL),
	m_maxEdgeLengthWithinSphere(NULL),
	m_oblateness(NULL)
{

	m_pointCoord.X = 0.0;
	m_pointCoord.Y = 0.0;

}

// Destructer
ObservationPoint::~ObservationPoint(){

	if( m_radius != NULL ){
		delete [] m_radius;
		m_radius = NULL;
	}
	
	if( m_maxEdgeLengthWithinSphere != NULL ){
		delete [] m_maxEdgeLengthWithinSphere;
		m_maxEdgeLengthWithinSphere = NULL;
	}

	if( m_oblateness != NULL ){
		delete [] m_oblateness;
		m_oblateness = NULL;
	}

}

// Read data of observation point from input file
void ObservationPoint::readObservationPointData( std::ifstream& ifs ){

	ifs >> m_pointCoord.X >> m_pointCoord.Y >> m_pointCoord.Z;

	std::cout << "Coordinate : " << m_pointCoord.X << " " << m_pointCoord.Y << " " << m_pointCoord.Z << std::endl;

	ifs >> m_numSpheres;
	if( m_numSpheres <= 0 ){
		std::cerr << "Error : Total number of spheres less than 1. : " << m_numSpheres << std::endl;
		exit(1);
	}
	std::cout << "Number of spheres : " << m_numSpheres << std::endl;

	m_radius = new double[m_numSpheres];
	m_maxEdgeLengthWithinSphere = new double[m_numSpheres];
	m_oblateness = new double[m_numSpheres];

	for( int i = 0; i < m_numSpheres; ++i ){
		ifs >> m_radius[i] >> m_maxEdgeLengthWithinSphere[i] >> m_oblateness[i];
		if( i > 0 && ( m_radius[i] < m_radius[i-1] ) ){
			std::cerr << "Error : Inner radius ( " << m_radius[i-1] << " ) is greater than outer radius ( " << m_radius[i] << " ) !! " << std::endl;
			exit(1);
		}
		if( i > 0 && ( m_maxEdgeLengthWithinSphere[i] < m_maxEdgeLengthWithinSphere[i-1] ) ){
			std::cerr << "Error : Inner edge length ( " << m_maxEdgeLengthWithinSphere[i-1] << " ) is greater than outer edge length ( " << m_maxEdgeLengthWithinSphere[i] << " ) !! " << std::endl;
			exit(1);
		}
		if( m_oblateness[i] < 0.0 || m_oblateness[i]> 1.0 ){
			std::cerr << "Error : Oblateness must be smaller than 1 and larger than 0." << std::endl;
			exit(1);
		}

	}

	for( int i = 0; i < m_numSpheres; ++i ){
		std::cout << " Radius[km] : " << m_radius[i] << ", edge length[km] : " << m_maxEdgeLengthWithinSphere[i] << ", oblateness : " << m_oblateness[i] << std::endl;
	}
	
}


// Calculate maximum length of specified coordinate on X-Y plane ( z = 0 ) 
double ObservationPoint::calcMaximumLengthOfPoint( const CommonParameters::XYZ& coord ) const{

	if( locateInSphere( coord, 0 ) ){
		return m_maxEdgeLengthWithinSphere[0];
	}
	
	assert( m_numSpheres > 0 );

	if( !locateInSphere( coord, m_numSpheres - 1 ) ){
		return 1.0e+20;
	}

	assert( m_numSpheres > 1 );

	int iSphere = m_numSpheres - 2;
	for( ; iSphere >= 1; --iSphere ){
		if( !locateInSphere( coord, iSphere ) ){
			break;
		}
	}

	const double vecX = coord.X - m_pointCoord.X;
	const double vecY = coord.Y - m_pointCoord.Y;
	const double vecZ = coord.Z - m_pointCoord.Z;

	double val = pow( vecX / m_radius[iSphere], 2 ) + pow( vecY / m_radius[iSphere], 2 ) + pow( vecZ / ( m_radius[iSphere] * ( 1.0 - m_oblateness[iSphere] ) ), 2 );
	val = sqrt(val);

	double val1 = sqrt( 3.0 * pow( m_radius[iSphere+1] / m_radius[iSphere], 2 ) );

	return m_maxEdgeLengthWithinSphere[iSphere] + ( m_maxEdgeLengthWithinSphere[iSphere+1] - m_maxEdgeLengthWithinSphere[iSphere] ) * ( val - 1.0 ) / ( val1 - 1.0 );
}

bool ObservationPoint::locateInSphere( const CommonParameters::XYZ& coord, const int iSphere ) const{

	assert( iSphere >= 0 || iSphere < m_numSpheres );

	const double vecX = coord.X - m_pointCoord.X;
	const double vecY = coord.Y - m_pointCoord.Y;
	const double vecZ = coord.Z - m_pointCoord.Z;

	const double depth = m_radius[iSphere] * ( 1.0 - m_oblateness[iSphere] );

	double val = pow( vecX / m_radius[iSphere], 2 )
			   + pow( vecY / m_radius[iSphere], 2 )
			   + pow( vecZ / depth, 2 );

	if( val <= 1.0 ){
		return true;
	}

	return false;

}
