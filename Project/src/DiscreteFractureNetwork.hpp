#ifndef DISCRETE_FRACTURE_NETWORK_H
#define DISCRETE_FRACTURE_NETWORK_H

#include <iostream>
#include "Eigen/Eigen"
#include <vector>
#include <list>

using namespace std;
using namespace Eigen;

namespace DFNLibrary{

struct DFN
{
    // FRATTURE
    unsigned int NumberFractures = 0; // Numero fratture
    vector<unsigned int> IdFractures = {} ; // Identificatori fratture --> intero positivo (dimensione 1)
    vector<Matrix3Xd> VerticesFractures = {}; // Matrici con vertici fratture (3xN) con N=numero vertici frattura

    // TRACCE
    unsigned int NumberTraces = 0;
    vector<unsigned int> IdTraces = {} ; // Identificatori fratture --> intero positivo (dimensione 1)
    vector<Vector2i> FractureTraces = {}; // Fratture associate a traccia
    vector<Vector<bool,2>> TipsTraces = {} ; // Tips booleano false= passante, true= non passante
    vector<Matrix<double,3,2>> VerticesTraces = {}; // vettore con estremi traccia
    vector<double> LengthTraces = {}; // vettore con lunghezza tracce

    vector<list<unsigned int>> P_Traces = {}; // Lista di identificatori tracce passanti di frattura ordinati per lunghezza
    vector<list<unsigned int>> NP_Traces = {}; // Lista di identificatori tracce NON passanti di frattura ordinati per lunghezza

    double tolerance = 100*numeric_limits<double>::epsilon();
};

struct PolygonalMesh
{
    unsigned int NumberCell0D = 0; // Numero totale Cell0D
    vector<unsigned int> IdCell0D = {}; // identificatori Cell0D --> intero positivo (dimensione 1)
    vector<Vector3d> CoordinatesCell0D = {}; // coordinate Cell0D --> (x,y,z) doubles (dimensione 3)

    unsigned int NumberCell1D = 0; // Numero totale Cell1D
    vector<unsigned int> IdCell1D = {}; // identificatori Cell1D --> intero positivo (dimensione 1)
    vector<Vector2i> VerticesCell1D = {}; // vertici Cell1D --> (id origine,id fine) (dimensione 2)

    unsigned int NumberCell2D = 0; // Numero totale Cell2D
    vector<unsigned int> IdCell2D = {}; // identificatori Cell2D --> intero positivo (dimensione 1)
    vector<list<unsigned int>> VerticesCell2D = {}; // vertici Cell2D --> lista di vertici in senso antiorario
    vector<list<unsigned int>> EdgesCell2D = {}; // lati Cell2D --> lista di lati in senso antiorario

    double tolerance = 1000*numeric_limits<double>::epsilon();
};


}


#endif
