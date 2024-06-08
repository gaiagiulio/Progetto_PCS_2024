#ifndef DFN_UTILS_H
#define DFN_UTILS_H

#include "DiscreteFractureNetwork.hpp"

namespace DFNLibrary{

struct DFN_functions
{
/************************************************************* SECONDARY FUNCTIONS **************************************************************************/

/**+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ GENERAL INLINE FUNCTIONS ++++++++++++++++++++++++++++++++++++++++**/

/** Restituisce vettore normale al piano passante per i 3 punti **/
    Vector3d NormalToPlane(Vector3d& p0,Vector3d& p1,Vector3d& p2);

/** Restituisce ascissa curvinea su P0+st del punto V
 *  V_P0: V-P0, differenza dei due punti**/
    double ascissa_curvilinea(Vector3d& V_P0,Vector3d& t);

/** Calcola l'intersezione tra le due rette r1: p1+t1*s1 e r2: p2+t2*s2.
 *  Restituisce come terza componentente 0 se le due rette NON si intersecano (sennò dà 1).
 *  Restituisce come prime componenti rispettivamente s1 e s2, ascisse curvilinee dell'intersezione su r1 e r2 (se si intersecano)
 *  tol: è la tolleranza usata per valutare intersezione **/
    Vector3d IntersectionBetweenLines(Vector3d& t1,Vector3d& t2, Vector3d& p1, Vector3d& p2, double& tol);

/**+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ PART 1 FUNCTIONS +++++++++++++++++++++++++++++++++++++++++++++**/

/** Restituisce un vettore con terza componente= 1 se la frattura interseca la retta r: x=P0+st. Le prime due componenti sono le ascisse curvilinee delle intersezioni con r
 *  terza componente = 0 se frattura NON interseca r
 *  quarta componente = 0 se intersezione è un lato completamente sulla retta (sennò =1)  **/
    Vector4d IntersectionFractureWithLine(DFNLibrary::DFN& dfn, const unsigned int & idFrac, Vector3d& P0, Vector3d& t, Vector3d& n);

/** Inserisce l'id della traccia nella lista delle tracce passanti o non per ciascuna delle fratture coinvolte (già ordinate per lunghezza descrescente **/
    void InsertSortedTraces(DFNLibrary::DFN& dfn, const unsigned int & frac, const unsigned int & id_tr, const bool & Tips, const double & length);


/**+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ PART 2 FUNCTIONS +++++++++++++++++++++++++++++++++++++++++++++**/

/** Restituisce l'id della cell1D esterna (lato esterno) su cui giace il punto.
 *  Se non trova alcuna cell1D esterna, restituisce -1 come primo elemento.
 *  Se il punto coincide con una cell0D esistente, restituisce prima componente -2 e come seconda componente l'id della cell0D**/
    Vector2i edge_to_traceExtreme(Vector3d& ext_tr,list<unsigned int>& external_edges,DFNLibrary::PolygonalMesh& frac);

/** Crea nuova cell0D nella PolygonalMesh con coordinate date (e inserisce in posizione iesima di ver_to_cells una lista vuota)
 *  frac: PolygonalMesh struct
 *  point: Vector3d con coordinate nuovo punto
 *  ver_to_cells: vector<list<unsigned int>> elemento i-esimo è la lista degli id delle celle associate alla cell0D di id i**/
    unsigned int NewCell0D(DFNLibrary::PolygonalMesh& frac,Vector3d& point, vector<list<unsigned int>>& ver_to_cells);

/** Crea nuova cell1D nella PolygonalMesh con estremi dati
 *  frac: PolygonalMesh struct
 *  ver1,ver2: unsigned int--> id vertici della cell1D
 *  edge_to_cells: vector<Vector2i> elemento i-esimo è il vettore degli id delle celle associate alla cell1D di id i (se la cell1D è esterna il secondo elemento è -1)**/
    unsigned int NewCell1D(DFNLibrary::PolygonalMesh& frac, unsigned int&  ver1,unsigned int& ver2, vector<Vector2i>& edge_to_cells);

/** Inserisce id di nuova cell1D (generata a partire dal lato edge) nella lista dei lati interni o esterni
 *  id_NEW_E: unsigned int --> id di nuova cell1D
 *  edge: unsigned int --> id di cell1D originaria
 *  external_edges,internal_edges: liste di unsigned int contenenti gli id dei lati esterni e interni**/
    void InternalExternalEdge(unsigned int& id_NEW_E,unsigned int& edge,list<unsigned int>& external_edges,list<unsigned int>& internal_edges);

/** Cerca intersezioni della traccia con i lati interni della cella corrente.
 *  Restituisce in ordine: l2 (id lato intersecato: unsigned int),  s2 (ascissa punto intersezione: double), v2 (id vertice intersecato (se interseca in vertice):  unsigned int),
 *  edge_found2 (bool),  cell_found (bool).
 *  edge_found2= true (1) se intersezione 2 è su lato, false (0) se è in un vertice.
 *  cell_found= true (1) se ho trovato il lato/vertice cercato, false (0) altrimenti.
 *  Parametri richiesti
 *  frac: frattura
 *  internal_edges: lista lati inteni
 *  t_T: vettore associato alla retta rT su cui giace la traccia rT: ext1_tr+t_T*s
 *  ext1_tr: estremo iniziale traccia
 *  c2D= id cella corrente
 *  edge_found0= true se intersez prec su lato
 *  l10,l20= id dei due lati coinvolti da iteraz prec (se edge_found0= true)
 *  v0= id vertice intersez precedente (se edge_found0= false **/
    Vector<double,5> IntersectCellEdges(DFNLibrary::PolygonalMesh& frac,list<unsigned int>& internal_edges,Vector3d& t_T,Vector3d& ext1_tr, unsigned int c2D, bool edge_found0, unsigned int l10, unsigned int l20, unsigned int v0, double s0);

 /** Crea i 2 nuovi vertici e i 2 nuovi lati, li inserisce nella mesh e ne restituisce gli id in ordine di inserimento sulla cella c (ordine antiorario).
 *  Nella funzione si aggiornano anche le variabili associate ai lati e ai vertici usate da calculate_fracture_cuts e viene aggiornata (se esiste) la cella adiacente.
 *  frac: PolygonalMesh struct
 *  external_edges,internal_edges: liste di unsigned int contenenti gli id dei lati esterni e interni
 *  edge_to_cells: vector<Vector2i> elemento i-esimo è il vettore degli id delle celle associate alla cell1D di id i (se la cell1D è esterna il secondo elemento è -1)
 *  ver_to_cells: vector<list<unsigned int>> elemento i-esimo è la lista degli id delle celle associate alla cell0D di id i
 *  c: id cella corrente
 *  l: id lato corrente
 *  ext1_tr, ext2_tr: coordinate dei due estremi della traccia (in ordine di incontro ciclando in senso antiorario sui lati di c)
 *  it_l: puntatore al lato corrente in frac.EdgesCell2D[c]
 *  it_ver: puntatore al vertice corrente in frac.VerticesCell2D[c] **/
    Vector<unsigned int,4> BookSpecialCase_EE(DFNLibrary::PolygonalMesh& frac,list<unsigned int>& external_edges,list<unsigned int>& internal_edges, vector<Vector2i>& edge_to_cells, vector<list<unsigned int>>& ver_to_cells, unsigned int& c,
                                               unsigned int& l, Vector3d& ext1_tr, Vector3d& ext2_tr, list<unsigned int>::iterator& it_l,list<unsigned int>::iterator& it_ver);

 /** Crea il nuovo vertice e il nuovo lato, li inserisce nella mesh e ne restituisce gli id in ordine di inserimento sulla cella c (ordine antiorario).
 *  Nella funzione si aggiornano anche le variabili associate ai lati e ai vertici usate da calculate_fracture_cuts e viene aggiornata (se esiste) la cella adiacente..
 *  frac: PolygonalMesh struct
 *  external_edges,internal_edges: liste di unsigned int contenenti gli id dei lati esterni e interni
 *  edge_to_cells: vector<Vector2i> elemento i-esimo è il vettore degli id delle celle associate alla cell1D di id i (se la cell1D è esterna il secondo elemento è -1)
 *  ver_to_cells: vector<list<unsigned int>> elemento i-esimo è la lista degli id delle celle associate alla cell0D di id i
 *  c: id cella corrente
 *  l: id lato corrente
 *  v_in_extr1: booleano. True se il vertice è il primo del lato (ciclando in senso antiorario sui lati della cella), false se è il secondo.
 *  ext1_tr: coordinate del nuovo estremo della traccia (cell0D non ancora esistente)
 *  v: id vertice con estremo traccia
 *  it_l: puntatore al lato corrente in frac.EdgesCell2D[c]
 *  it_ver: puntatore al vertice corrente in frac.VerticesCell2D[c] **/
    Vector<unsigned int, 2> BookSpecialCase_VE(DFNLibrary::PolygonalMesh& frac,list<unsigned int>& external_edges,list<unsigned int>& internal_edges, vector<Vector2i>& edge_to_cells, vector<list<unsigned int>>& ver_to_cells, unsigned int& c,
                                               unsigned int& l, bool& v_in_extr1, Vector3d& ext1_tr, unsigned int& v_end,list<unsigned int>::iterator& it_l,list<unsigned int>::iterator& it_ver);

 /** Restituisce vettore di due booleani.
 *  Prima componente True se si è nel caso di traccia sovrapposta a lato esistente con un vertice esistente e uno nuovo nel lato (book case)
 *  Seconda componente True se si è nel book case e il taglio è avvenuto con successo
 *  frac: PolygonalMesh struct
 *  external_edges,internal_edges: liste di unsigned int contenenti gli id dei lati esterni e interni
 *  edge_to_cells: vector<Vector2i> elemento i-esimo è il vettore degli id delle celle associate alla cell1D di id i (se la cell1D è esterna il secondo elemento è -1)
 *  ver_to_cells: vector<list<unsigned int>> elemento i-esimo è la lista degli id delle celle associate alla cell0D di id i
 *  l10, l2: id lati (servirà solo quello per cui edge_found=true)
 *  v0, v2: id lati (servirà solo quello per cui edge_found=false)
 *  ext_tr: coordinate nuovo estremo su lato (cell=d non ancora esistente)
 *  edge_found0, edge_found2: bool True se l'estremo della traccia nella cella è sul lato, false se è in vertice **/
    Vector<bool,2> GeneralBookCase(DFNLibrary::PolygonalMesh& frac,list<unsigned int>& external_edges,list<unsigned int>& internal_edges, vector<Vector2i>& edge_to_cells, vector<list<unsigned int>>& ver_to_cells, unsigned int& l10, unsigned int& l2,
                                    Vector3d& ext_tr, unsigned int& v0, unsigned int& v2, bool& edge_found0, bool& edge_found2);
/************************************************************* PRIMARY FUNCTIONS ******************************************************************************/

/**+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ PART 1 FUNCTIONS +++++++++++++++++++++++++++++++++++++++++++++**/

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

/**+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ PART 2 FUNCTIONS +++++++++++++++++++++++++++++++++++++++++++++**/

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


/** Genera i file inp per le cell 0D, 1D e 2D della PolygonalMesh frac, che rappresenta la frattura di indice id_frac tagliata.
 *  Questi file possono essere usati per visualizzare la mesh su Paraview.
 *  I file sono nominati Geometry0Ds_frac(id_frac).inp, Geometry1Ds_frac(id_frac).inp, Geometry2Ds_frac(id_frac).inp. **/
void CreateMeshFiles(DFNLibrary::PolygonalMesh& frac, const unsigned int& id_frac);

}

#endif
