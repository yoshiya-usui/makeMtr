//--------------------------------------------------------------------------
// MIT License
//
// Copyright (c) 2021 Yoshiya Usui
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//--------------------------------------------------------------------------
#include "ObservationLine.h"
#include "math.h"
#include <stdlib.h>
#include <fstream>
#include <iostream>

// Default constructer
ObservationLine::ObservationLine():
	m_numLayers(0),
	m_radius(NULL),
	m_maxEdgeLengthWithinLayer(NULL)
{
}

// Destructer
ObservationLine::~ObservationLine(){
	if( m_radius != NULL ){
		delete [] m_radius;
		m_radius = NULL;
	}
	
	if( m_maxEdgeLengthWithinLayer != NULL ){
		delete [] m_maxEdgeLengthWithinLayer;
		m_maxEdgeLengthWithinLayer = NULL;
	}
}

// Calculate maximum length of specified coordinate on X-Y plane ( z = 0 ) 
double ObservationLine::calcMaximumLengthOfPoint( const CommonParameters::XYZ& coord ) const{

	double edgeLeng = 1.0e+20;
	for( int iPoint = 0; iPoint < 2; ++iPoint ){// Loop of the end points
		const double vecX = coord.X - m_endPoints[iPoint].X;
		const double vecY = coord.Y - m_endPoints[iPoint].Y;
		const double vecZ = coord.Z - m_endPoints[iPoint].Z;

		const double radius = sqrt( pow( vecX, 2.0 ) + pow( vecY, 2.0 ) + pow( vecZ, 2.0 ) );

		if( radius <= m_radius[0] ){
			edgeLeng = std::min( m_maxEdgeLengthWithinLayer[0], edgeLeng );
			break;
		}

		for( int i = m_numLayers - 1; i >= 1; --i ){
			if( radius <= m_radius[i] ){
				 edgeLeng = m_maxEdgeLengthWithinLayer[i-1] + ( radius - m_radius[i-1] ) / ( m_radius[i] - m_radius[i-1] ) * ( m_maxEdgeLengthWithinLayer[i] - m_maxEdgeLengthWithinLayer[i-1] );
			}else{
				break;
			}
		}
	}

	const double lineLength = sqrt( 
		pow(m_endPoints[1].X - m_endPoints[0].X, 2) + 
		pow(m_endPoints[1].Y - m_endPoints[0].Y, 2) + 
		pow(m_endPoints[1].Z - m_endPoints[0].Z, 2) );
	const double innerProduct =
		(coord.X - m_endPoints[0].X) * (m_endPoints[1].X - m_endPoints[0].X) + 
		(coord.Y - m_endPoints[0].Y) * (m_endPoints[1].Y - m_endPoints[0].Y) + 
		(coord.Z - m_endPoints[0].Z) * (m_endPoints[1].Z - m_endPoints[0].Z);
	const double innerProductDiv = innerProduct / lineLength;
	if( innerProductDiv >= 0.0 && innerProductDiv <= lineLength ){// Within the two end points
		const CommonParameters::XYZ vectorAlongLine = {
			innerProductDiv / lineLength * (m_endPoints[1].X - m_endPoints[0].X),
			innerProductDiv / lineLength * (m_endPoints[1].Y - m_endPoints[0].Y),
			innerProductDiv / lineLength * (m_endPoints[1].Z - m_endPoints[0].Z) };
		const CommonParameters::XYZ vectorPerpendicularToLine = {
			coord.X - m_endPoints[0].X - vectorAlongLine.X,
			coord.Y - m_endPoints[0].Y - vectorAlongLine.Y,
			coord.Z - m_endPoints[0].Z - vectorAlongLine.Z };
		const double distanceFromLine = sqrt(pow(vectorPerpendicularToLine.X, 2) + pow(vectorPerpendicularToLine.Y, 2) + pow(vectorPerpendicularToLine.Z, 2) );
		if( distanceFromLine <= m_radius[0] ){
			return std::min( m_maxEdgeLengthWithinLayer[0], edgeLeng );
		}

		for( int i = m_numLayers - 1; i >= 1; --i ){
			if( distanceFromLine <= m_radius[i] ){
				edgeLeng = m_maxEdgeLengthWithinLayer[i-1] + ( distanceFromLine - m_radius[i-1] ) / ( m_radius[i] - m_radius[i-1] ) * ( m_maxEdgeLengthWithinLayer[i] - m_maxEdgeLengthWithinLayer[i-1] );
			}else{
				break;
			}
		}
	}
	return edgeLeng;

}

// Read data of observation points consisting line from input file
void ObservationLine::readObservationLineData( std::ifstream& ifs ){

	for( int i = 0; i < 2; ++i ){
		ifs >> m_endPoints[i].X >> m_endPoints[i].Y >> m_endPoints[i].Z;
	}

	ifs >> m_numLayers;
	if( m_numLayers <= 0 ){
		std::cerr << "Error : Total number of layers is less than 1 : " << m_numLayers << std::endl;
		exit(1);
	}

	m_radius = new double[m_numLayers];
	m_maxEdgeLengthWithinLayer = new double[m_numLayers];

	for( int i = 0; i < m_numLayers; ++i ){
		ifs >> m_radius[i] >> m_maxEdgeLengthWithinLayer[i];
		if( i > 0 && ( m_radius[i] < m_radius[i-1] ) ){
			std::cerr << "Error : Inner radius ( " << m_radius[i-1] << " ) is greater than outer radius ( " << m_radius[i] << " ) !! " << std::endl;
			exit(1);
		}
		if( i > 0 && ( m_maxEdgeLengthWithinLayer[i] < m_maxEdgeLengthWithinLayer[i-1] ) ){
			std::cerr << "Error : Inner edge length ( " << m_maxEdgeLengthWithinLayer[i-1] << " ) is greater than outer edge length ( " << m_maxEdgeLengthWithinLayer[i] << " ) !! " << std::endl;
			exit(1);
		}
	}

	std::cout << " First coordinate  : " << m_endPoints[0].X << " " << m_endPoints[0].Y << " " << m_endPoints[0].Z << std::endl;
	std::cout << " Second coordinate : " << m_endPoints[1].X << " " << m_endPoints[1].Y << " " << m_endPoints[1].Z << std::endl;
	std::cout << " Number of layers  : " << m_numLayers << std::endl;
	for( int i = 0; i < m_numLayers; ++i ){
		std::cout << " Radius [km] : " << m_radius[i] << ", Edge length [km] : " << m_maxEdgeLengthWithinLayer[i] << std::endl;
	}

}
