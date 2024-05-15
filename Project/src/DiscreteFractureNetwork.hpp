// DEFINISCO STRUTTURA DFN

#ifndef DISCRETE_FRACTURE_NETWORK_H
#define DISCRETE_FRACTURE_NETWORK_H


namespace DFNLibrary{



/**struct PolygonalMesh
{
    unsigned int NumberCell0D = 0; // Numero totale Cell0D
    vector<unsigned int> IdCell0D = {}; // identificatori Cell0D --> intero positivo (dimensione 1)
    vector<Vector2d> CoordinatesCell0D = {}; // coordinate Cell0D --> (x,y) doubles (dimensione 2)
    map<unsigned int, list<unsigned int>> MarkersCell0D = {}; // markers Cell0D --> chiave = marker, valore = lista punti con quel marker

    unsigned int NumberCell1D = 0; // Numero totale Cell1D
    vector<unsigned int> IdCell1D = {}; // identificatori Cell1D --> intero positivo (dimensione 1)
    vector<Vector2i> VerticesCell1D = {}; // vertici Cell1D --> (id origine,id fine) (dimensione 2)
    map<unsigned int, list<unsigned int>> MarkersCell1D = {}; // markers Cell1D --> chiave = marker, valore = lista lati con quel marker

    unsigned int NumberCell2D = 0; // Numero totale Cell2D
    vector<unsigned int> IdCell2D = {}; // identificatori Cell2D --> intero positivo (dimensione 1)
    vector<VectorXi> VerticesCell2D = {}; // vertici Cell2D --> vettore di vertici in senso antiorario (dim variabile)
    vector<VectorXi> EdgesCell2D = {}; // lati Cell2D --> vettore di lati in senso antiorario (dim variabile)

    double tolerance= 10*numeric_limits<double>::epsilon(); // tolleranza per test su lunghezze e aree della mesh
};
**/




}

#endif
