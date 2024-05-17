// DEFINISCO STRUTTURA DFN

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
    vector<unsigned int> IdTraces = {} ; // Identificatori fratture --> intero positivo (dimensione 1) //passaggio da lista??
    vector<Vector2i> FractureTraces = {}; // Fratture associate a traccia // passaggio da lista??
    vector<Vector<bool,2>> TipsTraces = {} ; // Tips booleano false= passante, true= non passante
    vector<Matrix<double,3,2>> VerticesTraces = {}; // vettore con estremi traccia
    vector<double> LengthTraces = {}; // vettore con lunghezza tracce

    vector<list<unsigned int>> P_Traces = {}; // Lista di identificatori tracce passanti di frattura ordinati per lunghezza
    vector<list<unsigned int>> NP_Traces = {}; // Lista di identificatori tracce NON passanti di frattura ordinati per lunghezza

    // per P_traces e NP_Traces vanno bene anche forwardlist??

    double tolerance = 100*numeric_limits<double>::epsilon(); // tolleranza che utente può inserire, aumento il default??
};


}


#endif
