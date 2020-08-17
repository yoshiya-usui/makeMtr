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
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <iomanip>
#include <math.h>
#include <algorithm>
#include <stdio.h>
#include <string.h>
#include "ObservationPoint.h"
#include "ObservationLine.h"

enum TypeID{
	EARTH = 0,
	AIR = 1
};

struct ObservedSite{
	double X;
	double Y;
	double Z;
	std::vector<double> radius;
	std::vector<double> edgeLength;
	double radiusInner;
	double radiusOuter;
	double edgeLengthInner;
};

struct Element{
	int nodeID[4];
	int attribute;
	int neigh[4];
};

struct FineRegionForSea{
	double X;
	double Y;
	double rotationAngle;
	std::vector<double> radius;
	std::vector<double> oblateness;
	std::vector<double> edgeLength;
};

const double PI = 3.14159265359;
const double EPS = 1.0e-12;

const double DEG2RAD = PI / 180.0;
const double RAD2DEG = 180.0 / PI;

CommonParameters::XYZ m_centerCoordOfFineSurface = { 0.0, 0.0, 0.0 };
CommonParameters::XYZ m_centerCoordOfFineRegion = { 0.0, 0.0, 0.0 };
//double m_oblateness[2] = { 0.0, 0.0 };

//double m_radiusInner = 10.0;// [km]
//double m_radiusOuter = 100.0;// [km]
//double m_edgeLengthInner = 0.1;// [km]
//double m_edgeLengthOuter = 100.0;// [km]

double m_rotationAngleFineSurface(0.0);
double m_rotationAngle(0.0);

std::set<int> m_regionAttrAirOrSea;

//double m_edgeLengthOnFineSurface = 0.0;
double m_radiusFineSurface = 0.0;
double m_edgeLengthFineSurface = 0.0;
double m_oblatenessHorizontalFineSurface = 0.0;
double m_oblatenessFineSurface[2] = { 0.0, 0.0 };

int m_numSphere = 0;
double* m_radius = NULL;
double* m_edgeLength = NULL;
double* m_oblatenessHorizontal = NULL;
double* m_oblateness[2] = { NULL, NULL };

int m_numObsSite = 0;
int m_numObsLine = 0;

ObservationPoint* m_obsSite = NULL;
ObservationLine* m_obsLine = NULL;

int m_numNodes = 0;
CommonParameters::XYZ* m_nodeCoords = NULL;

int m_numElements = 0;
Element* m_elements;
std::vector<int>* m_node2elem;

std::vector<int> m_attributesOfSea;
std::vector<FineRegionForSea> m_fineRegionForSea;

std::string extracFileNameWithoutExtension( const std::string fileNameFull );
void readParameters();
void readObsSites();
void readNodeData( const std::string& fileName );
void readElemData( const std::string& fileName );

void debugWriteNodeToVTK( const std::string& fileName );

void writeMtrFile( const std::string& fileNameRoot );
double calcEdgeLength( const int iNode );
bool locateInFineSurface( const CommonParameters::XYZ& coord );
bool locateInSphere( const CommonParameters::XYZ& coord, const int iSphere );
bool locateOnEarthSurface( const int iNode );

bool locateInTheSea( const int iNode );
double calcEdgeLengthInTheSea( const CommonParameters::XYZ& coord );

double calcEdgeLengthTransitionRegion( const CommonParameters::XYZ& coord, const int iSphere );
double calculateLengthOnEllipsoid( const double angleH, const double angleV, const int iSphere, const int iType );
double calcEdgeLengthNearSite( const CommonParameters::XYZ& coord );

int main( int argc, char* argv[] ){

	if( argc < 2 ){
		std::cerr << "You must specify node file name !!" << std::endl;
		exit(1);
	}

	std::string fileNameRoot = argv[1];

	readParameters();
	readObsSites();
	readNodeData( fileNameRoot + ".node" );
	readElemData( fileNameRoot + ".ele" );

	//debugWriteNodeToVTK( fileNameRoot + "_checkNode.vtk" );

	writeMtrFile( fileNameRoot );

	return 0;

}

std::string extracFileNameWithoutExtension( const std::string fileNameFull ){

	std::string::size_type pos;

    if( ( pos = fileNameFull.find_last_of(".") ) == std::string::npos){
        return fileNameFull;
    }
 
    return fileNameFull.substr(0, pos);
}

void readParameters(){

	std::ifstream ifs("makeMtr.param");

	if( !ifs.is_open() ){
		std::cerr << "Cannot open file makeMtr.param" << std::endl;
		exit(1);
	}

	std::cout << "Read parameters from makeMtr.param" << std::endl;	

#ifdef _FINE_SURF
	std::cout << "--- Fine surface ---"<< std::endl;

	ifs >> m_centerCoordOfFineSurface.X >> m_centerCoordOfFineSurface.Y >> m_centerCoordOfFineSurface.Z;
	std::cout << "Center coorinate of the fine region [km] : (X,Y,Z) = : " << m_centerCoordOfFineSurface.X << " " << m_centerCoordOfFineSurface.Y << " " << m_centerCoordOfFineSurface.Z << std::endl;

	ifs >> m_rotationAngleFineSurface;
	std::cout << "Rotation angle [deg] : " << m_rotationAngleFineSurface << std::endl;
	m_rotationAngleFineSurface *= DEG2RAD;

	int numAttrSeaOrAir(0);
	ifs >> numAttrSeaOrAir;
	std::cout << "Number of attributes of the sea or the air : " << numAttrSeaOrAir << std::endl;
	for( int i = 0; i < numAttrSeaOrAir; ++i ){
		int iAttr(0);
		ifs >> iAttr;
		m_regionAttrAirOrSea.insert(iAttr);
	}
	if( numAttrSeaOrAir > 0 ){
		std::cout << "List of attributes of the sea or the air : " << std::endl;
		for( std::set<int>::iterator itr = m_regionAttrAirOrSea.begin(); itr != m_regionAttrAirOrSea.end(); ++itr ){
			std::cout << std::setw(10) << *itr << std::endl;
		}
	}

	//ifs >> m_edgeLengthOnFineSurface;
	std::cout << "<Radius[km]> <Edge Length[km]> <OblatenessHorizontal> <OblatenessVerticalEarth> <OblatenessVerticalAir>" << std::endl;
	ifs >> m_radiusFineSurface >> m_edgeLengthFineSurface >> m_oblatenessHorizontalFineSurface >> m_oblatenessFineSurface[EARTH] >> m_oblatenessFineSurface[AIR];
	if( m_oblatenessFineSurface[EARTH] < 0 || m_oblatenessFineSurface[EARTH] > 1 ){
		std::cerr << "Oblateness must be smaller than 1 and larger than 0." << std::endl;
		exit(1);
	}
	if( m_oblatenessFineSurface[AIR] < 0 || m_oblatenessFineSurface[AIR] > 1 ){
		std::cerr << "Oblateness must be smaller than 1 and larger than 0." << std::endl;
		exit(1);
	}
	std::cout << std::setw(15) << std::scientific << m_radiusFineSurface;
	std::cout << std::setw(15) << std::scientific << m_edgeLengthFineSurface;
	std::cout << std::setw(15) << std::scientific << m_oblatenessHorizontalFineSurface;
	std::cout << std::setw(15) << std::scientific << m_oblatenessFineSurface[EARTH];
	std::cout << std::setw(15) << std::scientific << m_oblatenessFineSurface[AIR] << std::endl;
#endif

	//ifs >> m_oblatenessHorizontal;
	//if( m_oblatenessHorizontal < 0 || m_oblatenessHorizontal > 1 ){
	//	std::cerr << "Oblateness must be smaller than 1 and larger than 0." << std::endl;
	//	exit(1);
	//}
	//std::cout << "Oblateness of horizontal ellipsoid : " << m_oblatenessHorizontal << std::endl;

	std::cout << "--- Spheres ---"<< std::endl;

	ifs >> m_centerCoordOfFineRegion.X >> m_centerCoordOfFineRegion.Y >> m_centerCoordOfFineRegion.Z;
	std::cout << "Center coorinate of the fine region [km] : (X,Y,Z) = : " << m_centerCoordOfFineRegion.X << " " << m_centerCoordOfFineRegion.Y << " " << m_centerCoordOfFineRegion.Z << std::endl;

	ifs >> m_rotationAngle;
	m_rotationAngle *= DEG2RAD;
	std::cout << "Rotation angle of horizontal ellipsoid [deg] : " << m_rotationAngle * RAD2DEG << std::endl;

	ifs >> m_numSphere;
	std::cout << "Total number of spheres : " << m_numSphere << std::endl;

	if( m_numSphere < 1 ){
		std::cerr << "Total number of spheres must be greater than zero." << std::endl;
		exit(1);
	}

	m_radius = new double[m_numSphere];
	m_edgeLength = new double[m_numSphere];
	m_oblatenessHorizontal = new double[m_numSphere];
	m_oblateness[EARTH] = new double[m_numSphere];
	m_oblateness[AIR] = new double[m_numSphere];
	std::cout << "<Radius[km]> <Edge Length[km]> <OblatenessHorizontal> <OblatenessVerticalEarth> <OblatenessVerticalAir>" << std::endl;
	std::cout.precision(6);
	for( int i = 0; i < m_numSphere; ++i ){
		ifs >> m_radius[i] >> m_edgeLength[i] >> m_oblatenessHorizontal[i] >> m_oblateness[EARTH][i] >> m_oblateness[AIR][i];
		std::cout << std::setw(15) << std::scientific << m_radius[i];
		std::cout << std::setw(15) << std::scientific << m_edgeLength[i];
		std::cout << std::setw(15) << std::scientific << m_oblatenessHorizontal[i];
		std::cout << std::setw(15) << std::scientific << m_oblateness[EARTH][i];
		std::cout << std::setw(15) << std::scientific << m_oblateness[AIR][i] << std::endl;
	}
	for( int i = 1; i < m_numSphere; ++i ){
		if( m_radius[i] < m_radius[i-1] ){
			std::cerr << "Radius of the region " << i << " is smaller than that of the previous region." << std::endl;
			exit(1);
		}
		if( m_edgeLength[i] < m_edgeLength[i-1] ){
			std::cerr << "Edge length of the region " << i << " is smaller than that of the previous region." << std::endl;
			exit(1);
		}
		if( m_oblatenessHorizontal[i] < 0 || m_oblatenessHorizontal[i] > 1 ){
			std::cerr << "Oblateness of horizontal ellipsoid must be smaller than 1 and larger than 0." << std::endl;
			exit(1);
		}
		if( m_oblateness[EARTH][i] < 0 || m_oblateness[EARTH][i] > 1 ){
			std::cerr << "Oblateness must be smaller than 1 and larger than 0." << std::endl;
			exit(1);
		}
		if( m_oblateness[AIR][i] < 0 || m_oblateness[AIR][i] > 1 ){
			std::cerr << "Oblateness must be smaller than 1 and larger than 0." << std::endl;
			exit(1);
		}
		if( m_radius[i]*(1.0-m_oblatenessHorizontal[i]) < m_radius[i-1]*(1.0-m_oblatenessHorizontal[i-1] ) ){
			std::cerr << "Length of shorter axis of horizontal ellipsoid " << i << " is less than that of the previous ellipsoid." << std::endl;
			exit(1);
		}
		if( m_radius[i]*(1.0-m_oblateness[EARTH][i]) < m_radius[i-1]*(1.0-m_oblateness[EARTH][i-1] ) ){
			std::cerr << "Depth of sphere " << i << " is shallower than that of the previous sphere in the earth." << std::endl;
			exit(1);
		}
		if( m_radius[i]*(1.0-m_oblateness[AIR][i]) < m_radius[i-1]*(1.0-m_oblateness[AIR][i-1] ) ){
			std::cerr << "Depth of sphere " << i << " is shallower than that of the previous sphere in the air." << std::endl;
			exit(1);
		}
	}

	int numAttributesOfSea = 0;
	ifs >> numAttributesOfSea;
	for( int iAttr = 0; iAttr < numAttributesOfSea; ++iAttr ){
		int attrID(0);
		ifs >> attrID;
		m_attributesOfSea.push_back(attrID);
	}
	if( numAttributesOfSea > 0 ){
		std::cout << "Number of attributes of the sea : " << numAttributesOfSea << std::endl;
		for( std::vector<int>::const_iterator itr = m_attributesOfSea.begin(); itr != m_attributesOfSea.end(); ++itr ){
			std::cout << *itr << std::endl;
		}

		int numRegions = 0;
		ifs >> numRegions;
		std::cout << "Number of regions for definig the edge lengths in the sea : " << numRegions << std::endl;
		m_fineRegionForSea.reserve(numRegions);

		for( int iRegion = 0; iRegion < numRegions; ++iRegion ){
			FineRegionForSea tmpFineRegionForSea;
			ifs >> tmpFineRegionForSea.X >> tmpFineRegionForSea.Y;
			std::cout << "Center coorinate of the region [km] : (X,Y) = : " << tmpFineRegionForSea.X << " " << tmpFineRegionForSea.Y << std::endl;

			ifs >> tmpFineRegionForSea.rotationAngle;
			std::cout << "Rotation of the region [deg.] : " << tmpFineRegionForSea.rotationAngle << std::endl;
			tmpFineRegionForSea.rotationAngle *= DEG2RAD;

			int numRadius(0);
			ifs >> numRadius;
			std::cout << "Number of spheres definig the edge lengths : " << numRadius << std::endl;
			tmpFineRegionForSea.radius.reserve(numRadius);
			tmpFineRegionForSea.edgeLength.reserve(numRadius);
			tmpFineRegionForSea.oblateness.reserve(numRadius);

			std::cout << "<Radius[km]> <Edge Length[km]> <Oblateness>" << std::endl;
			for( int i = 0; i < numRadius; ++i ){
				double radius(0.0);
				double edgeLength(0.0);
				double oblateness(0.0);
				ifs >> radius >> edgeLength >> oblateness;
				std::cout << radius << " " << edgeLength << " " << oblateness << std::endl;
				tmpFineRegionForSea.radius.push_back(radius);
				tmpFineRegionForSea.edgeLength.push_back(edgeLength);
				tmpFineRegionForSea.oblateness.push_back(oblateness);
			}
			m_fineRegionForSea.push_back(tmpFineRegionForSea);
		}

		//for( std::vector<FineRegionForSea>::const_iterator itr = m_fineRegionForSea.begin(); itr != m_fineRegionForSea.end(); ++itr ){
		//	std::cout << "Center coorinate of the region [km] : (X,Y) = : " << itr->X << " " << itr->Y << std::endl;
		//	std::cout << "Rotation of the region [deg.] : " << itr->rotationAngle * RAD2DEG << std::endl;
		//	const int numRadius = static_cast<int>( itr->radius.size() );
		//	std::cout << "Number of spheres definig the edge lengths : " << numRadius << std::endl;
		//	std::cout << "<Radius[km]> <Edge Length[km]> <Oblateness>" << std::endl;
		//	for( int i = 0; i < numRadius; ++i ){
		//		std::cout << itr->radius[i]<< " " << itr->edgeLength[i] << " " << itr->oblateness[i] << std::endl;
		//	}
		//}		
	}

	ifs.close();
}

void readObsSites(){

	std::ifstream ifs("obs_site.dat");

	if( !ifs.is_open() ){
		std::cerr << "Cannot open file obs_site.dat" << std::endl;
		exit(1);
	}

	std::cout << "Read data of observed sites from obs_site.dat" << std::endl;	

	ifs >> m_numObsSite;

	std::cout << "Total number of observed site : " << m_numObsSite << std::endl;	

	m_obsSite = new ObservationPoint[m_numObsSite];

	for( int iObs = 0; iObs < m_numObsSite; ++iObs ){
		m_obsSite[iObs].readObservationPointData( ifs );
	}

	ifs >> m_numObsLine;

	std::cout << "Total number of lines : " << m_numObsLine << std::endl;	

	m_obsLine = new ObservationLine[m_numObsLine];

	for( int iObs = 0; iObs < m_numObsLine; ++iObs ){
		m_obsLine[iObs].readObservationLineData( ifs );
	}

	ifs.close();

}

void readNodeData( const std::string& fileName ){

	std::ifstream ifs( fileName.c_str() );

	if( !ifs.is_open() ){
		std::cerr << "Cannot open file " << fileName << std::endl;
		exit(1);
	}

	std::cout << "Read node data from " << fileName << std::endl;

	ifs >> m_numNodes;
	std::cout << "Total number of nodes : " << m_numNodes << std::endl;

	int idum(0);
	ifs >> idum >> idum >> idum;

	m_nodeCoords = new CommonParameters::XYZ[m_numNodes];
	for( int iNode = 0; iNode < m_numNodes; ++iNode ){
		ifs >> idum >> m_nodeCoords[iNode].X >> m_nodeCoords[iNode].Y >> m_nodeCoords[iNode].Z;
	}
	
	ifs.close();

}

// Read element data
void readElemData( const std::string& fileName ){

	std::ifstream ifs( fileName.c_str() );

	if( !ifs.is_open() ){
		std::cerr << "Cannot open file " << fileName.c_str() << std::endl;
		exit(1);
	}

	std::cout << "Read element data from " << fileName.c_str() << std::endl;

	ifs >> m_numElements;
	std::cout << "Total number of elements : " << m_numElements << std::endl;
	
	m_elements = new Element[m_numElements];

	int nodesPerTet(0);
	ifs >> nodesPerTet;

	if( nodesPerTet != 4 ){
		std::cerr << "Nodes number of a tetrahedron must be 4 !!" << std::endl;
		exit(1);
	}

	int iAttribute(0);
	ifs >> iAttribute;
	if( iAttribute != 1 ){
		std::cerr << "You must specify region attribute to elements." << std::endl;
		exit(1);
	}

	m_node2elem = new std::vector<int>[m_numNodes];

	for( int iElem = 0; iElem < m_numElements; ++iElem ){

		int elemID(0);
		ifs >> elemID;

		if( iElem + 1 != elemID ){
			std::cerr << "Element ID must be sequence number from 1!!" << std::endl;
			exit(1);
		}

		for( int i = 0; i < 4; ++i ){
			int nodeID(0);
			ifs >> nodeID;
			m_elements[iElem].nodeID[i] = nodeID;
			m_node2elem[nodeID-1].push_back( elemID );
		}

		ifs >> m_elements[iElem].attribute;

	}

	ifs.close();

}
void writeMtrFile( const std::string& fileNameRoot ){
	
	const std::string fileNameMtr = fileNameRoot + ".mtr";

	std::ofstream ofs( fileNameMtr.c_str() );

	if( !ofs.is_open() ){
		std::cerr << "Cannot open file " << fileNameMtr.c_str() << std::endl;
		exit(1);
	}

	std::cout << "Write mesh size to " << fileNameMtr.c_str() << std::endl;

	ofs << std::setw(10) << m_numNodes << std::setw(10) << 1 << std::endl;

	ofs.precision(6);
	std::cout.precision(6);
	for( int iNode = 0; iNode < m_numNodes; ++iNode ){
		//std::cout << std::setw(15) << std::scientific << m_nodeCoords[iNode].X;
		//std::cout << std::setw(15) << std::scientific << m_nodeCoords[iNode].Y;
		//std::cout << std::setw(15) << std::scientific << m_nodeCoords[iNode].Z << std::endl;

		const double edgeLength = calcEdgeLength( iNode );
		if( edgeLength < EPS ){
			std::cout << "Edge length is too short " << edgeLength << std::endl;
			exit(1);
		}
		ofs << std::setw(15) << std::scientific << edgeLength << std::endl;

		//std::cout << std::setw(15) << std::scientific << calcEdgeLength( m_nodeCoords[iNode] ) << std::endl;
	}

	ofs.close();

}

double calcEdgeLength( const int iNode ){

	const CommonParameters::XYZ coord = m_nodeCoords[iNode];

	double minVal = calcEdgeLengthNearSite(coord);

#ifdef _FINE_SURF
	if( locateInFineSurface( coord ) && locateOnEarthSurface( iNode ) ){
		if( m_edgeLengthFineSurface < minVal ){
			minVal = m_edgeLengthFineSurface;
		}
	}
#endif

	if( locateInTheSea( iNode ) ){
		minVal = std::min( calcEdgeLengthInTheSea( coord ), minVal );
	}

	if( locateInSphere( coord, 0 ) ){
		return std::min( m_edgeLength[0], minVal );
	}
	else if( !locateInSphere( coord, m_numSphere-1 ) ){
		return std::min( m_edgeLength[m_numSphere-1], minVal );
	}
	else{
		for( int iSphere = 1; iSphere < m_numSphere; ++iSphere ){
			if( locateInSphere( coord, iSphere ) ){
				return std::min( calcEdgeLengthTransitionRegion( coord, iSphere ), minVal );
				break;
			}
		}
	}

	std::cerr << "Wrong!!" << std::endl;
	exit(1);

	return -1.0;

}

bool locateInFineSurface( const CommonParameters::XYZ& coord ){

	const double vecXOrg = coord.X - m_centerCoordOfFineSurface.X;
	const double vecYOrg = coord.Y - m_centerCoordOfFineSurface.Y;
	// Coordinate transform
	const double vecX = vecXOrg * cos( - m_rotationAngleFineSurface ) - vecYOrg * sin( - m_rotationAngleFineSurface );
	const double vecY = vecXOrg * sin( - m_rotationAngleFineSurface ) + vecYOrg * cos( - m_rotationAngleFineSurface );
	const double vecZ = coord.Z - m_centerCoordOfFineSurface.Z;

	const int iType = vecZ < 0 ? AIR : EARTH;

	const double longAxisLength = m_radiusFineSurface;
	const double shortAxisLength = longAxisLength * ( 1.0 - m_oblatenessHorizontalFineSurface );
	const double depth = longAxisLength * ( 1.0 - m_oblatenessFineSurface[iType] );

	double val = pow( vecX / longAxisLength, 2 )
			   + pow( vecY / shortAxisLength, 2 )
			   + pow( vecZ / depth, 2 );
	
	//std::cout << "coord : " << coord.X << " " << coord.Y << " " << coord.Z << std::endl;
	//std::cout << "m_radius[iSphere] depth : " << m_radius[iSphere] << " " << depth << std::endl;
	//std::cout << "iSphere val : " << iSphere << " " << val << std::endl;

	if( val <= 1.0 ){
		return true;
	}

	return false;

}

bool locateInSphere( const CommonParameters::XYZ& coord, const int iSphere ){
	
	if( iSphere < 0 || iSphere >= m_numSphere ){
		std::cerr << "Wrong shepre ID:  " << iSphere << std::endl;
		exit(1);
	}

	const double vecXOrg = coord.X - m_centerCoordOfFineRegion.X;
	const double vecYOrg = coord.Y - m_centerCoordOfFineRegion.Y;
	// Coordinate transform
	const double vecX = vecXOrg * cos( - m_rotationAngle ) - vecYOrg * sin( - m_rotationAngle );
	const double vecY = vecXOrg * sin( - m_rotationAngle ) + vecYOrg * cos( - m_rotationAngle );
	const double vecZ = coord.Z - m_centerCoordOfFineRegion.Z;

	const int iType = vecZ < 0.0 ? AIR : EARTH;

	const double longAxisLength = m_radius[iSphere];
	const double shortAxisLength = longAxisLength * ( 1.0 - m_oblatenessHorizontal[iSphere] );
	const double depth = longAxisLength * ( 1.0 - m_oblateness[iType][iSphere] );

	double val = pow( vecX / longAxisLength, 2 )
			   + pow( vecY / shortAxisLength, 2 )
			   + pow( vecZ / depth, 2 );
	
	//std::cout << "coord : " << coord.X << " " << coord.Y << " " << coord.Z << std::endl;
	//std::cout << "m_radius[iSphere] depth : " << m_radius[iSphere] << " " << depth << std::endl;
	//std::cout << "iSphere val : " << iSphere << " " << val << std::endl;

	if( val <= 1.0 ){
		return true;
	}

	return false;

}

bool locateOnEarthSurface( const int iNode ){

	int icountLand(0);
	int icountAirOrSea(0);
	for( std::vector<int>::iterator itr = m_node2elem[iNode].begin(); itr != m_node2elem[iNode].end(); ++itr ){

		//std::cout << m_elements[*itr].attribute << std::endl;
		if( m_regionAttrAirOrSea.find( m_elements[ (*itr) - 1 ].attribute ) == m_regionAttrAirOrSea.end() ){
			++icountLand;
		}
		else{
			++icountAirOrSea;
		}

	}

	if( icountLand > 0 && icountAirOrSea > 0 ){
		//std::cout << "TRUE " << m_nodeCoords[iNode].X << " " << m_nodeCoords[iNode].Y << " " << m_nodeCoords[iNode].Z << std::endl;
		return true;
	}

	//std::cout << "FALSE " << m_nodeCoords[iNode].X << " " << m_nodeCoords[iNode].Y << " " << m_nodeCoords[iNode].Z << std::endl;
	return false;

}

bool locateInTheSea( const int iNode ){

	for( std::vector<int>::const_iterator itr = m_node2elem[iNode].begin(); itr != m_node2elem[iNode].end(); ++itr ){
		const int iElem = *itr - 1;
		for( std::vector<int>::const_iterator itr = m_attributesOfSea.begin(); itr != m_attributesOfSea.end(); ++itr ){
			if( m_elements[iElem].attribute == *itr ){
				return true;
			}
		}
	}

	return false;
}

double calcEdgeLengthInTheSea( const CommonParameters::XYZ& coord ){

	double edgeLength(1.0e+20);

	for( std::vector<FineRegionForSea>::const_iterator itr = m_fineRegionForSea.begin(); itr != m_fineRegionForSea.end(); ++itr ){

		const int numRadius = static_cast<int>( itr->radius.size() );
		for( int i = 0; i < numRadius; ++i ){
			const double vecXOrg = coord.X - itr->X;
			const double vecYOrg = coord.Y - itr->Y;

			// Coordinate transform
			const double angle = - itr->rotationAngle ;
			const double vecX = vecXOrg * cos(angle) - vecYOrg * sin(angle);
			const double vecY = vecXOrg * sin(angle) + vecYOrg * cos(angle);

			const double longAxisLength = itr->radius[i];
			const double shortAxisLength = longAxisLength * ( 1.0 - itr->oblateness[i] );

			double val = pow( vecX / longAxisLength, 2 ) + pow( vecY / shortAxisLength, 2 );
			if( val <= 1.0 ){
				edgeLength = std::min( edgeLength, itr->edgeLength[i] );
				break;
			}
		}

	}

	return edgeLength;
}

double calcEdgeLengthTransitionRegion( const CommonParameters::XYZ& coord, const int iSphere ){

	if( iSphere < 1 || iSphere >= m_numSphere ){
		std::cerr << "Wrong shepre ID:  " << iSphere << std::endl;
		exit(1);
	}

	////const double radiusInner = m_radius[iSphere-1];
	////const double radiusOuter = m_radius[iSphere];
	////const double edgeLengthInner = m_edgeLength[iSphere-1];
	////const double edgeLengthOuter = m_edgeLength[iSphere];
	////const double oblatenessInner = m_oblateness[iType][iSphere-1];
	////const double oblatenessOuter = m_oblateness[iType][iSphere];

	////const double radiusHorizontal = hypot( coord.X - m_centerCoordOfFineRegion.X, coord.Y - m_centerCoordOfFineRegion.Y );
	////const double radius = hypot( radiusHorizontal, coord.Z - m_centerCoordOfFineRegion.Z );
	////const double slope = fabs( coord.Z - m_centerCoordOfFineRegion.Z ) / radiusHorizontal;

	////const double radius0 = radiusInner * ( 1.0 - oblatenessInner ) * sqrt( ( slope*slope + 1.0 ) / ( slope*slope + (1.0-oblatenessInner)*(1.0-oblatenessInner) ) );
	////const double radius1 = radiusOuter * ( 1.0 - oblatenessOuter ) * sqrt( ( slope*slope + 1.0 ) / ( slope*slope + (1.0-oblatenessOuter)*(1.0-oblatenessOuter) ) );

	////if( radius1 < radius0 ){
	////	std::cerr << "Specified coordinate ( " << coord.X << ", " << coord.Y << ", " << coord.Z << " ) does't locate in the trasition region" << std::endl;
	////	exit(1);
	////}

	////if( radius < radius0 ){
	////	std::cerr << "Specified coordinate ( " << coord.X << ", " << coord.Y << ", " << coord.Z << " ) does't locate in the trasition region" << std::endl;
	////	std::cerr << "iSphere : " << iSphere << std::endl;
	////	std::cerr << "radiusHorizontal : " << radiusHorizontal << std::endl;
	////	std::cerr << "radius : " << radius << std::endl;
	////	std::cerr << "radius0 : " << radius0 << std::endl;
	////	std::cerr << "radius1 : " << radius1 << std::endl;
	////	exit(1);
	////}

	////return ( radius - radius0 ) / ( radius1 - radius0 ) * ( edgeLengthOuter - edgeLengthInner ) + edgeLengthInner;

	//const double vecXOrg = coord.X - m_centerCoordOfFineRegion.X;
	//const double vecYOrg = coord.Y - m_centerCoordOfFineRegion.Y;
	//// Coordinate transform
	//const double vecX = vecXOrg * cos( - m_rotationAngle ) - vecYOrg * sin( - m_rotationAngle );
	//const double vecY = vecXOrg * sin( - m_rotationAngle ) + vecYOrg * cos( - m_rotationAngle );
	//const double vecZ = coord.Z - m_centerCoordOfFineRegion.Z;

	//const double longAxisLength = m_radius[iSphere-1];
	//const double shortAxisLength = longAxisLength * ( 1.0 - m_oblatenessHorizontal );
	//const double depth = m_radius[iSphere-1] * ( 1.0 - m_oblateness[iType][iSphere-1] );

	//double val = pow( vecX / longAxisLength, 2 )
	//		   + pow( vecY / shortAxisLength, 2 )
	//		   + pow( vecZ / depth, 2 );
	//val = sqrt(val);

	//const double longAxisLength1 = m_radius[iSphere];
	//const double shortAxisLength1 = longAxisLength * ( 1.0 - m_oblatenessHorizontal );
	//const double depth1 = m_radius[iSphere] * ( 1.0 - m_oblateness[iType][iSphere] );

	//double val1 = pow( longAxisLength1 / longAxisLength, 2 )
	//		    + pow( shortAxisLength1 / shortAxisLength, 2 )
	//		    + pow( depth1 / depth, 2 );
	//val1 = sqrt(val1);

	//return m_radius[iSphere-1] + ( m_radius[iSphere] - m_radius[iSphere-1] ) * ( val - 1.0 ) / ( val1 - 1.0 );

	const double vecXOrg = coord.X - m_centerCoordOfFineRegion.X;
	const double vecYOrg = coord.Y - m_centerCoordOfFineRegion.Y;
	// Coordinate transform
	const double vecX = vecXOrg * cos( - m_rotationAngle ) - vecYOrg * sin( - m_rotationAngle );
	const double vecY = vecXOrg * sin( - m_rotationAngle ) + vecYOrg * cos( - m_rotationAngle );
	const double vecZ = coord.Z - m_centerCoordOfFineRegion.Z;

	const int iType = vecZ < 0.0 ? AIR : EARTH;

	const double angleHorizontal = atan2( vecY, vecX );
	const double lengthHorizontal = hypot( vecY, vecX );
	const double angleVertical = atan2( vecZ, lengthHorizontal );
	const double length = hypot( lengthHorizontal, vecZ );

	const double length0 = calculateLengthOnEllipsoid( angleHorizontal, angleVertical, iSphere-1, iType );
	const double length1 = calculateLengthOnEllipsoid( angleHorizontal, angleVertical, iSphere, iType );

	const double eps = 1.0e-10;
	if( length <= length0 - eps || length >= length1 + eps ){
		std::cerr << "Specified coordinate ( " << coord.X << ", " << coord.Y << ", " << coord.Z << " ) does't locate in the trasition region" << std::endl;
		std::cerr << "iSphere : " << iSphere << std::endl;

		std::cerr << "vecX : " << vecX << std::endl;
		std::cerr << "vecY : " << vecY << std::endl;
		std::cerr << "vecZ : " << vecZ << std::endl;

		std::cerr << "lengthHorizontal      : " << lengthHorizontal << std::endl;
		std::cerr << "angleH : " << angleHorizontal << std::endl;
		std::cerr << "angleV : " << angleVertical << std::endl;

		std::cerr << "length  : " << length << std::endl;
		std::cerr << "length0 : " << length0 << std::endl;
		std::cerr << "length1 : " << length1 << std::endl;
		exit(1);
	}

	return ( length - length0 ) / ( length1 - length0 ) * ( m_edgeLength[iSphere] - m_edgeLength[iSphere-1] ) + m_edgeLength[iSphere-1];

}

// Calculate length on ellipsoid
double calculateLengthOnEllipsoid( const double angleH, const double angleV, const int iSphere, const int iType ){

	if( iSphere < 0 || iSphere >= m_numSphere ){
		std::cerr << "Wrong shepre ID:  " << iSphere << std::endl;
		exit(1);
	}

	if( angleH < -PI || angleH > PI ){
		std::cerr << "Horizontal angle is improper : " << angleH << std::endl;
		exit(1);
	}

	if( angleV < -PI || angleV > PI ){
		std::cerr << "Vertical angle is improper : " << angleV << std::endl;
		exit(1);
	}

	const double longAxisLength = m_radius[iSphere];
	const double shortAxisLength = longAxisLength * ( 1.0 - m_oblatenessHorizontal[iSphere] );
	const double verticalLength = longAxisLength * ( 1.0 - m_oblateness[iType][iSphere] );
	const double eps = 1.0e-9;

	double lengthH(-1.0);
	if( fabs(angleH - PI*0.5) < eps || fabs(angleH + PI*0.5) < eps ){
		lengthH = shortAxisLength;
	}
	else{
		const double constValH = longAxisLength * shortAxisLength / hypot( shortAxisLength, longAxisLength * tan(angleH) );
		lengthH = hypot( constValH, constValH * tan(angleH) );
	}
	
	if( fabs(angleV - PI*0.5) < eps || fabs(angleV + PI*0.5) < eps ){
		return verticalLength;
	}
	else{
		const double constValV = lengthH * verticalLength / hypot( verticalLength, lengthH * tan(angleV) );
		return hypot( constValV, constValV * tan(angleV) );
	}

}

double calcEdgeLengthNearSite( const CommonParameters::XYZ& coord ){

	double edgeLength(1.0e+20);

	for( int iObs = 0; iObs < m_numObsSite; ++iObs ){

		//const CommonParameters::XY coordXY = { coord.X, coord.Y };

		const double length = m_obsSite[iObs].calcMaximumLengthOfPoint( coord );

		if( length < edgeLength ){
			edgeLength = length;
		}

	}

	for( int iObs = 0; iObs < m_numObsLine; ++iObs ){
		const double length = m_obsLine[iObs].calcMaximumLengthOfPoint( coord );

		if( length < edgeLength ){
			edgeLength = length;
		}

	}

	return edgeLength;

}

// Debug write nodes data to VTK file
void debugWriteNodeToVTK( const std::string& fileName ){

	std::ofstream ofsVTK( fileName.c_str() );
	if( !ofsVTK ) {
		std::cout << " Error : Cannot open file " << fileName.c_str() << std::endl;
		exit(1);
	}

	ofsVTK << "# vtk DataFile Version 2.0" << std::endl;
	ofsVTK << "NodeList" << std::endl;
	ofsVTK << "ASCII" << std::endl;
	ofsVTK << "DATASET UNSTRUCTURED_GRID" << std::endl;
	
	ofsVTK.precision(9);
	ofsVTK << std::fixed;

	ofsVTK << "POINTS " << m_numNodes << " float" << std::endl;
	for( int iNode = 0; iNode < m_numNodes; ++iNode ){
		ofsVTK << m_nodeCoords[iNode].X << " " << m_nodeCoords[iNode].Y << " " << m_nodeCoords[iNode].Z << std::endl;
	}
	
	ofsVTK << "CELLS " << m_numNodes << " " << m_numNodes * 2 << std::endl;
	for( int i = 0; i < m_numNodes; ++i ){
		ofsVTK << 1 << " " << i << std::endl;
	}

	ofsVTK << "CELL_TYPES " << m_numNodes << std::endl;
	for( int i = 0; i < m_numNodes; ++i ){
		ofsVTK << "1" << std::endl;
	}

	ofsVTK << "CELL_DATA " << m_numNodes << std::endl;
	ofsVTK << "SCALARS Location int" << std::endl;
	ofsVTK << "LOOKUP_TABLE default" << std::endl;
	for( int iNode = 0; iNode < m_numNodes; ++iNode ){
		ofsVTK << ( locateOnEarthSurface(iNode) ? 1 : 0 ) << std::endl;
	}

	ofsVTK << "SCALARS NumLand int" << std::endl;
	ofsVTK << "LOOKUP_TABLE default" << std::endl;
	for( int iNode = 0; iNode < m_numNodes; ++iNode ){
		int icountLand(0);
		for( std::vector<int>::iterator itr = m_node2elem[iNode].begin(); itr != m_node2elem[iNode].end(); ++itr ){
			if( m_regionAttrAirOrSea.find( m_elements[ (*itr) - 1 ].attribute ) == m_regionAttrAirOrSea.end() ){
				++icountLand;
			}
		}
		ofsVTK << icountLand << std::endl;
	}

	ofsVTK << "SCALARS NumSeaOrAir int" << std::endl;
	ofsVTK << "LOOKUP_TABLE default" << std::endl;
	for( int iNode = 0; iNode < m_numNodes; ++iNode ){
		int icountSea(0);
		for( std::vector<int>::iterator itr = m_node2elem[iNode].begin(); itr != m_node2elem[iNode].end(); ++itr ){
			if( m_regionAttrAirOrSea.find( m_elements[ (*itr) - 1 ].attribute ) != m_regionAttrAirOrSea.end() ){
				++icountSea;
			}
		}
		ofsVTK << icountSea << std::endl;
	}

	ofsVTK.close();

}


