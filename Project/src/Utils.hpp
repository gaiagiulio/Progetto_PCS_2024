#ifndef DFN_UTILS_H
#define DFN_UTILS_H

#include "DiscreteFractureNetwork.hpp"

namespace DFNLibrary{

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

/** Calcola per la frattura di vertici dati la Polygonal Mesh ottenuta compiendo i tagli lungo le sue tracce.
 *  L'ordine di taglio seguito è: prima tracce passanti poi non passanti, entrambe gli insiemi ordinati per lugnhezza decrescente.
 *  frac_vertices: matrice con vertici frattura (in senso antiorario)
 *  p_traces: lista di identificativi di tracce passanti ordinate per lunghezza decrescente
 *  np_traces: lista di identificativi di tracce non passanti ordinate per lunghezza decrescente
 *  traces_extremes: matrice con estremi delle tracce
 *  tol: tolleranza da command line (se non è inserita è 0) **/
PolygonalMesh calculate_fracture_cuts(Matrix3Xd& frac_vertices, list<unsigned int>& p_traces, list<unsigned int>& np_traces,
                                      vector<Matrix<double,3,2>>& traces_extremes, double tol);

}

#endif
