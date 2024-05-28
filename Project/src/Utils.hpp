#ifndef DFN_UTILS_H
#define DFN_UTILS_H

#include "DiscreteFractureNetwork.hpp"

namespace DFNLibrary{

struct DFN_functions
{
/************************************************************* SECONDARY FUNCTIONS ***************************************************/
/** Restituisce vettore normale al piano passante per i 3 punti **/
    Vector3d NormalToPlane(Vector3d& p0,Vector3d& p1,Vector3d& p2);

/** Restituisce un vettore con terza componente= 1 se la frattura interseca la retta r: x=P0+st. Le prime due componenti sono le ascisse curvilinee delle intersezioni con r
 *  terza componente = 0 se frattura NON interseca r **/
    Vector3d IntersectionFractureWithLine(DFNLibrary::DFN& dfn, const unsigned int & idFrac, Vector3d& P0, Vector3d& t, Vector3d& n);

/** Inserisce l'id della traccia nella lista delle tracce passanti o non per ciascuna delle fratture coinvolte (già ordinate per lunghezza descrescente **/
    void InsertSortedTraces(DFNLibrary::DFN& dfn, const unsigned int & frac, const unsigned int & id_tr, const bool & Tips, const double & length);

/************************************************************* MAIN FUNCTIONS ***************************************************/

 /** Importa fratture dal filepath (path + file name) ( and test if it's correct__ AGGIUNGI)
 *  dfn: DFN struct
 *  restituisce risultato della lettura (true if is success, false otherwise) **/
    bool ImportFractures(const string& filepath, DFN& dfn);


 /** Calcola tutte le tracce presenti nel DFN e le salva nello struct dfn.
 *  Per ogni frattura salva separatamente l'insieme delle tracce passanti e di quelle non passanti, per ordine descrescente di lunghezza
 *  dfn: DFN struct **/
    void calculateTraces(DFN& dfn);


 /** Stampa su un file di output il numero totale delle tracce e per ciscuna l'id, le due fratture coinvolte e le coordinate 3D degli estremi.
 *  dfn: DFN struct
 *  outputFile: stringa con nome file di output **/
    void PrintTraces(const string& outputFile, DFN& dfn);


/** Stampa su un file di output per ogni frattura il numero totale delle tracce, seguita dalle tracce (prima le passanti e poi le non passanti,
 *  ordinate per lunghezza decrescente), corredate di Id, Tips (false= passante, true= non passante) e lunghezza.
 *  dfn: DFN struct
 *  outputFile: stringa con nome file di output **/
    void PrintSortedFractureTraces(const string& outputFile, DFN& dfn);

};

struct PolygonalMesh_functions
{
/************************************************************* SECONDARY FUNCTIONS ***************************************************/

/** Restituisce l'id della cell1D su cui giace il punto. Se non trova alcuna cell1D restitisce -1 **/
int edge_to_traceExtreme(Vector3d& ext_tr,DFNLibrary::PolygonalMesh& frac);

/** Restitusce lista ordinata per ascissa curvilinea crescente (punti da ext1_tr a ext2_tr) delle intersezioni della traccia con estremi (ext1_tr, ext2_tr) con i lati interni.
 *  Elemento i-esimo = (lato intersecato,ascissa intersezione)**/
list<Vector2d> IntersectTraceWithInternalEdges(Vector3d& ext1_tr,Vector3d& ext2_tr,DFNLibrary::PolygonalMesh& frac,list<unsigned int>& internal_edges);


/** Crea nuova cell0D nella PolygonalMesh con coordinate date
 *  frac: PolygonalMesh struct
 *  point: Vector3d con coordinate nuovo punto **/
unsigned int NewCell0D(DFNLibrary::PolygonalMesh& frac,Vector3d& point);

/** Crea nuova cell1D nella PolygonalMesh con estremi dati
 *  frac: PolygonalMesh struct
 *  ver1,ver2: unsigned int--> id vertici della cell1D **/
unsigned int NewCell1D(DFNLibrary::PolygonalMesh& frac, unsigned int&  ver1,unsigned int& ver2);

/** Inserisce id di nuova cell1D (generata a partire dal lato edge) nella lista dei lati interni o esterni
 *  id_NEW_E: unsigned int --> id di nuova cell1D
 *  edge: unsigned int --> id di cell1D originaria
 *  external_edges,internal_edges: liste di unsigned int contenenti gli id dei lati**/
void InternalExternalEdge(unsigned int& id_NEW_E,unsigned int& edge,list<unsigned int>& external_edges,list<unsigned int>& internal_edges);

/** Effettua sulla frattura frac il taglio lungo la traccia di estremi ext1_tr e ext2_tr, inserenso nuove cell0D, cell1D, cell2D in frac. Restituisce false se errore nel processo, true altrimenti.
 * frac: PolygonalMesh struct
 * intersezioni: lista di coppie (id lato intersecato, ascissa curvilinea sulla retta ext1_tr + (ext2_tr-ext1_tr)t)
 * ext1_tr, ext2_tr: Vector3d--> estremi traccia
 * external_edges, internal_edges: liste di unsigned int contenenti gli id dei lati esterni e interni**/
bool cut_divided_trace(DFNLibrary::PolygonalMesh& frac,list<Vector2d>& intersezioni, Vector3d& ext1_tr, Vector3d& ext2_tr,list<unsigned int>& external_edges,list<unsigned int>& internal_edges);

/************************************************************* MAIN FUNCTIONS ***************************************************/


/** Calcola per la frattura di vertici dati la Polygonal Mesh ottenuta compiendo i tagli lungo le sue tracce.
 *  In caso di errore restituisce la frattura parzialmente tagliata, come si trova al momento dell'errore.
 *  L'ordine di taglio seguito è: prima tracce passanti poi non passanti, entrambe gli insiemi ordinati per lugnhezza decrescente.
 *  frac_vertices: matrice con vertici frattura (in senso antiorario)
 *  p_traces: lista di identificativi di tracce passanti ordinate per lunghezza decrescente
 *  np_traces: lista di identificativi di tracce non passanti ordinate per lunghezza decrescente
 *  traces_extremes: matrice con estremi delle tracce
 *  tol: tolleranza da command line (se non è inserita è 0) **/
PolygonalMesh calculate_fracture_cuts(Matrix3Xd& frac_vertices, list<unsigned int>& p_traces, list<unsigned int>& np_traces,
                                      vector<Matrix<double,3,2>>& traces_extremes, double tol);
};

}

#endif
