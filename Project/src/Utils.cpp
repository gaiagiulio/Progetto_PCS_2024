#include "Utils.hpp"
#include "DiscreteFractureNetwork.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

using namespace std;
using namespace Eigen;

/** Restituisce vettore normale al piano passante per i 3 punti **/
Vector3d NormalToPlane(Vector3d& p0,Vector3d& p1,Vector3d& p2)
{
    Vector3d v1 = p1-p0;
    Vector3d v2 = p2-p0;
    Vector3d n(v1[1]*v2[2]-v1[2]*v2[1],-v1[0]*v2[2]+v1[2]*v2[0], v1[0]*v2[1]-v1[1]*v2[0]);
    return n;
}

/** Restituisce un vettore con terza componente= 1 se la frattura interseca la retta r: x=P0+st. Le prime due componenti sono le ascisse curvilinee delle intersezioni con r
 *  terza componente = 0 se frattura NON interseca r **/
Vector3d IntersectionFractureWithLine(DFNLibrary::DFN& dfn, const unsigned int & idFrac, Vector3d& P0, Vector3d& t, Vector3d& n)
{
    Vector3d result(0,0,1);
    bool no_intersect = false;
    bool val_q1 = false;
    bool val_q2 = false;
    bool sign_zero_first = false;
    bool sign_zero = false;
    bool sign;
    bool sign_first;

    Matrix3Xd& ver = dfn.VerticesFractures[idFrac];
    unsigned int numVertices= ver.cols();
    // Valuto il primo vertice rispetto alla retta d'intersezione
    Vector3d V_P0 = ver(all,0)-P0;
    double prod = n[0]*(V_P0[1]*t[2]-V_P0[2]*t[1])+n[1]*(-V_P0[0]*t[2]+V_P0[2]*t[0])+n[2]*(V_P0[0]*t[1]-V_P0[1]*t[0]); //prodotto misto n, V_P0 (vettore che da P0 va nel vertice V) e t
    if (abs(prod) < dfn.tolerance) // primo vertice su r
    {
        bool sign_zero_first = true;
        bool sign_zero = true;
        val_q1 = true;
        result[0] = (V_P0[0]/t[0] + V_P0[1]/t[1] + V_P0[2]/t[2])/3 ; //calcola ascissa curvilinea di v su retta r:  st=  V_P0. NB: calcolo s in tutte e tre le componenti e faccio la media per avere un risultato robusto
    }
    else
    {
        sign_first = signbit(prod); // true se prod negativo
        sign = sign_first;
    }
    // Valuto gli altri vertici della frattura per vedere se è intersecata da r e trovare i punti d'intersezione
    for (unsigned int k=1; k<numVertices; k++)
    {
        V_P0 = ver(all,k)-P0;
        prod = n[0]*(V_P0[1]*t[2]-V_P0[2]*t[1])+n[1]*(-V_P0[0]*t[2]+V_P0[2]*t[0])+n[2]*(V_P0[0]*t[1]-V_P0[1]*t[0]);
        if (abs(prod) < dfn.tolerance) // vertice sulla retta r
        {
            if (sign_zero == true) // vertice precedente su retta
            {
                val_q2 = true;
                result[1] = (V_P0[0]/t[0] + V_P0[1]/t[1] + V_P0[2]/t[2])/3 ; // ascissa curvilinea di v su r
                break; // ho trovato le due intersezioni: esco dal ciclo sui vertici
            }
            else
            {
                if (val_q1 == false)
                {
                    val_q1 = true;
                    result[0] = (V_P0[0]/t[0] + V_P0[1]/t[1] + V_P0[2]/t[2])/3 ;
                    sign_zero = true;
                }
                else
                {
                    val_q2 = true;
                    result[1]= (V_P0[0]/t[0] + V_P0[1]/t[1] + V_P0[2]/t[2])/3 ; // ascissa curvilinea di v su r
                    break; // ho trovato le due intersezioni: esco dal ciclo sui vertici
                }
            }
        }
        else // vertice NON è sulla retta r
        {
            bool s = signbit(prod);
            if (sign == s) // vertice sullo STESSO lato di r rispetto al precedente
            {
                if (sign_zero == true) // vertice = unica intersezione con r
                {
                    no_intersect=true;
                    break;
                }
                continue;
            }
            else // vertice sul lato OPPOSTO di r rispetto al precedente
            {
                if (sign_zero == true) // vertice precedente su retta
                    continue;
                else // vertice precedente non su retta
                {
                    // calcolo intersezione tra r e retta r2 passante per il vertice e il precedente con Mx=b con M=(t,Pk-P(k-1)) e b= P(k-1)-P0
                    Vector3d t2=ver(all,k)-ver(all,k-1);
                    Matrix<double,3,2> M{{t[0],t2[0]},{t[1],t2[1]},{t[2],t2[2]}};
                    Vector3d b2= ver(all,k-1)-P0;
                    Vector2d x = M.fullPivLu().solve(b2); // vettore con prima componente ascissa curvilinea su r, seconda componente ascissa curvilinea su r2
                    if (val_q1 == false)
                    {
                        val_q1 = true;
                        result[0] = x[0];
                    }
                    else
                    {
                        val_q2 = true;
                        result[1] = x[0] ; // ascissa curvilinea di v su r
                        break; // ho trovato le due intersezioni: esco dal ciclo sui vertici
                    }
                }
            }
            sign_zero= false;
            sign = s;
        }
    }

    // gestione dei casi dove unico punto di intersezione con r è il primo vertice (A) e dove la seconda intersezione è tra il primo e l'ultimo vertice (B)
    if ((val_q1==true) && (no_intersect==false))
        if (val_q2==false)
        {
            if (sign_zero_first==true) // caso A
                no_intersect=true;
            else // caso B
            {
                // calcolo intersezione tra r e retta r2 passante per l'ultimo vertice e il primo // con Mx=b con M=(t,Pn-P1) e b= P1-P0
                Vector3d t2=ver(all,numVertices-1)-ver(all,0);
                Matrix<double,3,2> M{{t[0],t2[0]},{t[1],t2[1]},{t[2],t2[2]}};
                Vector3d b2= ver(all,0)-P0;
                Vector2d x = M.fullPivLu().solve(b2);
                val_q2 = true;
                result[1]= x[0];
            }
        }
    if (val_q2==false)
        no_intersect = true;
    if (no_intersect == true)
        result[2]=0;
    return result;
}

/** Inserisce l'id della traccia nella lista delle tracce passanti o non per ciascuna delle fratture coinvolte (già ordinate per lunghezza descrescente **/
void InsertSortedTraces(DFNLibrary::DFN& dfn, const unsigned int & frac, const unsigned int & id_tr, const bool & Tips, const double & length)
{
    bool inserted = false;
    if (Tips) // traccia NON passante su frattura
    {
        //inserisco in NP_Traces[frac] in base a lunghezza
        for (auto it = dfn.NP_Traces[frac].begin(); it != dfn.NP_Traces[frac].end();it++)
        {
            if ((length > dfn.LengthTraces[*(it)]) && (not inserted))
            {
                dfn.NP_Traces[frac].insert(it,id_tr);
                inserted = true;
            }
        }
        if (not inserted)
            dfn.NP_Traces[frac].push_back(id_tr);
    }
    else // traccia passante su frattura
    {
        //inserisco in P_Traces[frac] in base a lunghezza
        for (auto it = dfn.P_Traces[frac].begin(); it != dfn.P_Traces[frac].end();it++)
        {
            if ((length > dfn.LengthTraces[*(it)]) && (not inserted))
            {
                dfn.P_Traces[frac].insert(it,id_tr);
                inserted = true;
            }
        }
        if (not inserted)
            dfn.P_Traces[frac].push_back(id_tr);
    }
}

/** Restituisce l'id della cell1D su cui giace il punto. Se non trova alcuna cell1D restitisce -1 **/
unsigned int edge_to_traceExtreme(Vector3d& ext_tr,DFNLibrary::PolygonalMesh& frac)
{
    bool l_found= false;
    unsigned int l;
    for (unsigned int i=0; i<frac.NumberCell1D; i++)
    {
        unsigned int v0 = frac.VerticesCell1D[i][0];  //id vertice 1
        unsigned int v1 = frac.VerticesCell1D[i][1];  //id vertice 2

        Vector3d v1_v0 = frac.CoordinatesCell0D[v1]-frac.CoordinatesCell0D[v0]; // vettore da v0 a v1
        Vector3d ext_v0 = ext_tr-frac.CoordinatesCell0D[v0]; // vettore da v0 a estremo
        Vector3d n_ext_edge(v1_v0[1]*ext_v0[2]-v1_v0[2]*ext_v0[1],-v1_v0[0]*ext_v0[2]+v1_v0[2]*ext_v0[0], v1_v0[0]*ext_v0[1]-v1_v0[1]*ext_v0[0]); // prodotto vettoriale (cioè normale al piano passante per v0,v1 e estremo)
        if (n_ext_edge[0]*n_ext_edge[0]+n_ext_edge[1]*n_ext_edge[1]+n_ext_edge[2]*n_ext_edge[2] <= frac.tolerance) // estremo sulla retta contenente il lato (perchè n_ext_edge norma nulla)
        {
            double s_ext = (ext_v0[0]/v1_v0[0] + ext_v0[1]/v1_v0[1] + ext_v0[2]/v1_v0[2])/3; // ascissa curvilinea dell'estremo sulla retta R: V0+(V1-V0)s (retta su cui poggia il lato)
            if ((s_ext>(-frac.tolerance)) && (s_ext<(1+frac.tolerance))) // s in [0,1], cioè estremo nel lato
            {
                l=i;
                l_found= true;
            }
        }

        if (l_found) // esco dal ciclo quando trovo il lato
            break;
    }
    if (not l_found)
    {
        cerr << "Considered point in NOT on an edge (cell1D)" << endl;
        return -1;
    }
    else
        return l;
}

/** Restitusce lista ordinata per ascissa curvilinea crescente (punti da ext1_tr a ext2_tr) delle intersezioni della traccia con estremi (ext1_tr, ext2_tr) con i lati interni.
 *  Elemento i-esimo = (lato intersecato,ascissa intersezione)**/
list<Vector2d> IntersectTraceWithInternalEdges(Vector3d& ext1_tr,Vector3d& ext2_tr,DFNLibrary::PolygonalMesh& frac,list<unsigned int>& internal_edges)
{
    list<Vector2d> intersezioni={}; //lista ordinata di intersezioni traccia (allontanandosi da ext1_tr verso ext2_tr). Ogni vettore primo elemento: lato intersecante, secondo elemento: ascissa su traccia
    Vector3d t_T = ext2_tr - ext1_tr; // vettore tangente a rT: ext1_tr + t_T*s (retta su cui giace la traccia)
    for (auto it_edge = internal_edges.begin(); it_edge != internal_edges.end();it_edge++)
    {
        unsigned int id_edge = *(it_edge); // identificatore lato corrente

        // Calcolo intersezione retta rL associata al lato e retta rT associata alla traccia
        Vector3d t_L = frac.CoordinatesCell0D[frac.VerticesCell1D[id_edge][1]] - frac.CoordinatesCell0D[frac.VerticesCell1D[id_edge][0]] ; // vettore tangente a rL
        Vector3d pv_TL(t_T[1]*t_L[2]-t_T[2]*t_L[1],t_T[0]*t_L[2]+t_T[2]*t_L[0], t_T[0]*t_L[1]-t_T[1]*t_L[0]); // prodotto vettoriale tra t_T e t_L
        if ((pv_TL[0]*pv_TL[0]+pv_TL[1]*pv_TL[1]+pv_TL[2]*pv_TL[2])>frac.tolerance) // se t_T e t_L NON sono paralleli (norma pv_TL NON nulla) calcolo possibile intersezione
        {
            Matrix<double,3,2> M{{t_T[0],t_L[0]},{t_T[1],t_L[1]},{t_T[2],t_L[2]}};
            Vector3d b = frac.CoordinatesCell0D[frac.VerticesCell1D[id_edge][0]] - ext1_tr ;
            Vector2d sol_intersez = M.fullPivLu().solve(b);
            double s = sol_intersez[0]; // ascissa intersezione su rT: ext1_tr + t_T*s
            if ((s>frac.tolerance) && (s<(1-frac.tolerance))) //intersezione nel segmento della traccia
            {
                bool inserted = false;
                for (auto it_intsz = intersezioni.begin(); it_intsz != intersezioni.end();it_intsz++) // ciclo per inserirlo in ordine
                {
                    Vector2d elemento = *(it_intsz);
                    if ((s < elemento[1]) && (not inserted))
                    {
                        intersezioni.insert(it_intsz,{id_edge,s}); // inserisci per ascissa crescente
                        inserted = true;
                    }
                }
                if (not inserted)
                    intersezioni.push_back({id_edge,s}); // aggiungi alla fine della lista
            }
        }
    }
    return intersezioni;
}

/** Crea nuova cell0D nella PolygonalMesh con coordinate date
 *  frac: PolygonalMesh struct
 *  point: Vector3d con coordinate nuovo punto **/
unsigned int NewCell0D(DFNLibrary::PolygonalMesh& frac,Vector3d& point)
{
    unsigned int id_NEW_V = frac.NumberCell0D;
    frac.NumberCell0D += 1; // aumento numero cell0D
    frac.IdCell0D.push_back(id_NEW_V); // inserisco il nuovo id nella mesh
    frac.CoordinatesCell0D.push_back(point);
    return id_NEW_V;
}

/** Crea nuova cell1D nella PolygonalMesh con estremi dati
 *  frac: PolygonalMesh struct
 *  ver1,ver2: unsigned int--> id vertici della cell1D **/
unsigned int NewCell1D(DFNLibrary::PolygonalMesh& frac, unsigned int&  ver1,unsigned int& ver2)
{
    unsigned int id_NEW_E = frac.NumberCell1D;
    frac.NumberCell1D += 1; // aumento numero cell1D
    frac.IdCell1D.push_back(id_NEW_E); // inserisco il nuovo id nella mesh
    frac.VerticesCell1D.push_back({ver1,ver2});
    return id_NEW_E;
}

namespace DFNLibrary{

bool ImportFractures(const string& filepath, DFN& dfn)
{
    ifstream file;
    file.open(filepath);

    if(file.fail())
    {
        cerr << "Error while opening file Fractures" << endl;
        return false;
    }

    string header;
    getline(file,header);

    string line;
    getline(file,line);
    istringstream convert(line);
    convert >> dfn.NumberFractures; // leggo numero totale fratture

    if (dfn.NumberFractures==0)
    {
        cerr << "There is no fracture" << endl;
        return false;
    }

    // Imposto lunghezza dei vettori associati alle fratture
    dfn.IdFractures.reserve(dfn.NumberFractures);
    dfn.VerticesFractures.reserve(dfn.NumberFractures);

    char sep;
    unsigned int id;

    for (unsigned int i=0; i < dfn.NumberFractures; i++)
    {
        getline(file,header);

        unsigned int numVertices;
        getline(file,line);
        istringstream convert(line);
        convert >> id >> sep >> numVertices;

        getline(file,header);

        MatrixXd vertices(3,numVertices);
        for (unsigned int j=0; j <3; j++)
        {
            getline(file,line);
            istringstream convert(line);
            convert >> vertices(j,0);
            for (unsigned int v=1 ; v < numVertices; v++)
            {
                convert >> sep >> vertices(j,v);
            }
        }

        dfn.IdFractures.push_back(id); // inserisco id frattura
        dfn.VerticesFractures.push_back(vertices); // // inserisco vertici frattura
    };

    file.close();

    // Controllo banale da togliere!! o modificare
    for (unsigned int i=0; i < dfn.NumberFractures; i++)
    {
        cout << "id frattura: " << dfn.IdFractures[i] << endl;
    }

    return true;
}


void calculateTraces(DFN& dfn)
{
    // CALCOLARE TRACCE incrociando tutte le fratture a due a due
    unsigned int num_fractures = dfn.NumberFractures;
    // Reserve dei vettori associati alle tracce: Scelgo numero grande perchè non so quante saranno. (a seguire metodo banale per regolare N grande in base al numero di fratture nel DFN)
    unsigned int N=0;
    if (num_fractures <= 20)
        N = ceil(num_fractures*num_fractures/2) ; //numero di fratture se ciascuna si intersecasse con tutte le altre
    else if (num_fractures <= 50)
        N = ceil(num_fractures*num_fractures/4); //numero di fratture se ciascuna si intersecasse con metà delle altre
    else if (num_fractures <= 200)
        N = ceil(num_fractures*num_fractures/10); //numero di fratture se ciascuna si intersecasse con 1/5 delle altre
    else
        N = 5000;

    dfn.IdTraces.reserve(N);
    dfn.FractureTraces.reserve(N);
    dfn.TipsTraces.reserve(N);
    dfn.VerticesTraces.reserve(N);
    dfn.LengthTraces.reserve(N);

    dfn.P_Traces.resize(num_fractures); // NB: check se serve inizializzare con le liste vuote
    dfn.NP_Traces.resize(num_fractures);

    // Variabile interna alla funzione, di supporto, con dati da NON ricalcolare per ogni frattura (calcolati solo per i=0)
    // normale al piano n (Vector3d), d (double), centroide (Vector3d), raggio max (double) --> Vector8d
    vector<Vector<double,8>> FractureUsefulData;
    FractureUsefulData.resize(num_fractures);

    for (unsigned int i=0; i <num_fractures; i++)
    {
        unsigned int frac1 = dfn.IdFractures[i] ;
        Matrix3Xd& ver1 = dfn.VerticesFractures[frac1]; // matrice dei vertici della frattura 1
        //Cerco il piano della frattura 1 : calcolo prendendo i primi tre vertici (per la prima iterazione, sennò lo leggo da FractureUsefulData)
        Vector3d n1;
        double d1;
        if (i==0)
        {
            Vector3d p10 = ver1(all,0);
            Vector3d p11 = ver1(all,1);
            Vector3d p12 = ver1(all,2);

            n1 = NormalToPlane(p10,p11,p12); // vettore normale a piano 1
            d1 = n1.dot(p10); // piano 1 x: dot(n1,x)=d1
        }
        else
        {
            n1 << FractureUsefulData[i][0],FractureUsefulData[i][1],FractureUsefulData[i][2];
            d1 =  FractureUsefulData[i][3];
        }

        // Calcolo centroide e raggio per frattura 1 (se i=0 lo calcolo, sennò lo leggo da FractureUsefulData[i])
        Vector3d c1; // centroide frattura 1 (baricentro geometrico)
        double radius1=0; // raggio massimo di sfera contenente la frattura 1, di centro c1

        if (i==0)
        {
            unsigned int numVertices= ver1.cols();
            c1 << 0,0,0;
            for (unsigned int vert=0; vert < numVertices; vert++)
            {
                for (unsigned int coor=0; coor<3; coor++)
                    c1[coor] += ver1(coor,vert);
            }
            c1 = c1/numVertices;

            for (unsigned int vert=0; vert < numVertices; vert++)
            {
                double r_vert=(ver1(0,vert)-c1[0])*(ver1(0,vert)-c1[0])+(ver1(1,vert)-c1[1])*(ver1(1,vert)-c1[1])+(ver1(2,vert)-c1[2])*(ver1(2,vert)-c1[2]); // quadrato della distanza del vertice da c1
                if (r_vert > radius1)
                    radius1 = r_vert;
            }
            radius1 = sqrt(radius1); // raggio massimo
        }
        else
        {
            c1 << FractureUsefulData[i][4],FractureUsefulData[i][5],FractureUsefulData[i][6];
            radius1 = FractureUsefulData[i][7];
        }

        for (unsigned int j=i+1; j <num_fractures; j++)
        {
            unsigned int frac2 = dfn.IdFractures[j];
            Matrix3Xd& ver2 = dfn.VerticesFractures[frac2]; // matrice dei vertici della frattura 2

            // CONTROLLO PER ESCLUSIONE CASI (lo attivo solo per i>0, perchè devo calcolare tutti i piani)

            // Calcolo centroide e raggio per frattura 2 (se i=0 lo calcolo, sennò lo leggo da FractureUsefulData[j])
            Vector3d c2; // centroide frattura 2
            double radius2=0; // raggio massimo di sfera contenente la frattura 2, di centro c2

            if (i==0)
            {
                unsigned int numVertices= ver2.cols();
                c2 << 0,0,0;
                for (unsigned int vert=0; vert < numVertices; vert++)
                {
                    for (unsigned int coor=0; coor<3; coor++)
                        c2[coor] += ver2(coor,vert);
                }
                c2 = c2/numVertices;

                for (unsigned int vert=0; vert < numVertices; vert++)
                {
                    double r_vert=(ver2(0,vert)-c2[0])*(ver2(0,vert)-c2[0])+(ver2(1,vert)-c2[1])*(ver2(1,vert)-c2[1])+(ver2(2,vert)-c2[2])*(ver2(2,vert)-c2[2]); // quadrato della distanza del vertice da c1
                    if (r_vert > radius2)
                        radius2 = r_vert;
                }
                radius2 = sqrt(radius2); // raggio massimo
                // Inserisco in FractureUsefulData[j] (per iterazioni successive)
                for (unsigned int index=4; index<7; index++)
                    FractureUsefulData[j][index] = c2[index-4];
                FractureUsefulData[j][7] = radius2;

            }
            else
            {
                c2 << FractureUsefulData[j][4],FractureUsefulData[j][5],FractureUsefulData[j][6];
                radius2 = FractureUsefulData[j][7];
            }

            // Meccanismo di esclusione: passo a iterazione successiva se le due sfere non si intersecano
            // NB: lo attivo solo per i>0, perchè devo calcolare tutti i piani. Per 0 lo attivo in un secondo momento.
            if (i>0)
            {
                double dist_centroids = (c1[0]-c2[0])*(c1[0]-c2[0])+(c1[1]-c2[1])*(c1[1]-c2[1])+(c1[2]-c2[2])*(c1[2]-c2[2]);
                if (dist_centroids > radius1*radius1 + radius2*radius2 + 2*radius1*radius2) // d>r1+r2
                    continue;
            }

            //Cerco il piano della frattura 2 : calcolo prendendo i primi tre vertici (per i=0, sennò lo leggo da FractureUsefulData)

            Vector3d n2; // vettore normale a piano 2
            double d2; // piano 2 x: dot(n2,x)=d2
            if (i==0)
            {
                Vector3d p20 = ver2(all,0);
                Vector3d p21 = ver2(all,1);
                Vector3d p22 = ver2(all,2);

                n2 = NormalToPlane(p20,p21,p22);
                d2 = n2.dot(p20);
                // Inserisco in FractureUsefulData[j] (per iterazioni successive)
                for (unsigned int index=0; index<3; index++)
                    FractureUsefulData[j][index] = n2[index];
                FractureUsefulData[j][3] = d2;
            }
            else
            {
                n2 << FractureUsefulData[j][0],FractureUsefulData[j][1],FractureUsefulData[j][2];
                d2 =  FractureUsefulData[j][3];
            }

            if (i==0) // controllo per esclusione alla prima iterazione
            {
                double dist_centroids = (c1[0]-c2[0])*(c1[0]-c2[0])+(c1[1]-c2[1])*(c1[1]-c2[1])+(c1[2]-c2[2])*(c1[2]-c2[2]);
                if (dist_centroids > radius1*radius1 + radius2*radius2 + 2*radius1*radius2) // d>r1+r2
                    continue;
            }


            // Calcola retta d'intersezione tra piani r: x = P0 +st

            Vector3d t(n1[1]*n2[2]-n1[2]*n2[1],-n1[0]*n2[2]+n1[2]*n2[0], n1[0]*n2[1]-n1[1]*n2[0]); // t: vettore tangente a retta

            // Esclusione di piani paralleli (t ha norma nulla)
            double t_norm = t[0]*t[0] + t[1]*t[1] + t[2]*t[2];
            if (t_norm < dfn.tolerance)
                continue;

            Matrix3d A{{n1[0],n1[1],n1[2]},{n2[0],n2[1],n2[2]},{t[0],t[1],t[2]}};
            Vector3d b(d1,d2,0);
            Vector3d P0 = A.fullPivLu().solve(b); // P0

            // Calcolo se le fratture sono intersecate da r: x = P0 +st
            // Frattura 1
            Vector3d int1= IntersectionFractureWithLine(dfn, frac1, P0, t, n1); // vettore con (q1,q2,flag) con q1 e q2 ascisse di intersezioni se esistono
            if (int1[2]<numeric_limits<double>::epsilon()) // uso tol con epsiolon di macchina perchè è in vettore di double (anche se non dovrei avere problemi)
                continue; // passo a considerare altra coppia di fratture se una non interseca la retta
            //Frattura 2
            Vector3d int2= IntersectionFractureWithLine(dfn, frac2, P0, t, n2); // vettore con (q1,q2,flag) con q1 e q2 ascisse di intersezioni se esistono
            if (int2[2]<numeric_limits<double>::epsilon()) // uso tol con epsiolon di macchina perchè è in vettore di double (anche se non dovrei avere problemi)
                continue;

            // CALCOLO TRACCE (calcolo estremi tracce e per ciascuna frattura se siano passanti o meno) --> uso le ascisse curvilinee

            // I1=[q1,q2] intervallo di intersezione della retta con la frattura 1
            double q1 = min(int1[0],int1[1]);
            double q2 = max(int1[0],int1[1]);
            // I2=[q3,q4] intervallo di intersezione della retta con la frattura 2
            double q3 = min(int2[0],int2[1]);
            double q4 = max(int2[0],int2[1]);

            bool Tips1 = false; // false= traccia passante per la frattura 1, false= traccia NON passante per la frattura 1
            bool Tips2 = false; // false= traccia passante per la frattura 2, false= traccia NON passante per la frattura 2

            // Calcolo (se esiste) traccia di estremi P1 e P2 (ascisse curvilinee su retta r)
            double P1;
            double P2;

            if (abs(min(q1,q3)-q1)/max(max(abs(min(q1,q3)),abs(q1)),1.) < dfn.tolerance) // min(q1,q3)=q1
            {
                if (abs(min(q2,q3)-q2)/max(max(abs(min(q2,q3)),abs(q2)),1.) < dfn.tolerance) // min(q2,q3)=q2
                    continue;
                else // min(q2,q3)=q3
                {
                    P1 = q3;
                    if (abs(q1-q3)/max(max(abs(q1),abs(q3)),1.)<dfn.tolerance) // q1 = q3
                    {
                        bool test3 = (abs(q2-q4)/max(max(abs(q2),abs(q4)),1.)<dfn.tolerance); // q2 = q4
                        if (abs(min(q2,q4)-q2)/max(max(abs(min(q2,q4)),abs(q2)),1.) < dfn.tolerance) // min(q2,q4)=q2
                        {
                            P2 = q2;
                            if (not test3) // q2 != q4
                                Tips2 = true;
                        }
                        else // min(q2,q4)=q4
                        {
                            P2 = q4;
                            if (not test3) // q2 != q4
                                Tips1 = true;
                        }
                    }
                    else // q1 != q3
                    {
                        Tips1 = true;
                        if (abs(min(q2,q4)-q2)/max(max(abs(min(q2,q4)),abs(q2)),1.) < dfn.tolerance) // min(q2,q4)=q2
                        {
                            if (abs(q2-q3)/max(max(abs(q2),abs(q3)),1.)<dfn.tolerance) // q2 = q3
                                continue;
                            else // q2 != q3
                            {
                                P2 = q2;
                                if (abs(q2-q4)/max(max(abs(q2),abs(q4)),1.)>dfn.tolerance) // q2 != q4
                                    Tips2 = true;
                            }
                        }
                        else // min(q2,q4)=q4
                            P2 = q4;
                    }
                }
            }
            else // min(q1,q3)=q3
            {
                if (abs(min(q1,q4)-q4)/max(max(abs(min(q1,q4)),abs(q4)),1.) < dfn.tolerance) // min(q1,q4)=q4
                    continue;
                else // min(q1,q4)=q1
                {
                    P1 = q1;
                    if (abs(q1-q3)/max(max(abs(q1),abs(q3)),1.)<dfn.tolerance) // q1 = q3
                    {
                        bool test3 = (abs(q2-q4)/max(max(abs(q2),abs(q4)),1.)<dfn.tolerance); // q2 = q4
                        if (abs(min(q2,q4)-q4)/max(max(abs(min(q2,q4)),abs(q4)),1.) < dfn.tolerance) // min(q2,q4)=q4
                        {
                            P2 = q4;
                            if (not test3) // q2 != q4
                                Tips1 = true;
                        }
                        else // min(q2,q4)=q2
                        {
                            P2 = q2;
                            if (not test3) // q2 != q4
                                Tips2 = true;
                        }
                    }
                    else // q1 != q3
                    {
                        Tips2 = true;
                        if (abs(min(q2,q4)-q4)/max(max(abs(min(q2,q4)),abs(q4)),1.) < dfn.tolerance) // min(q2,q4)=q4
                        {
                            if (abs(q2-q3)/max(max(abs(q2),abs(q3)),1.)<dfn.tolerance) // q2 = q3
                                continue;
                            else // q2 != q3
                            {
                                P2 = q4;
                                if (abs(q2-q4)/max(max(abs(q2),abs(q4)),1.)>dfn.tolerance) // q2 != q4
                                    Tips1 = true;
                            }
                        }
                        else // min(q2,q4)=q2
                            P2 = q2;
                    }
                }
            }

            // Inserisci traccia trovata (se non è stata trovata si è incappati in continue)
            dfn.NumberTraces += 1;
            unsigned int id_trac = dfn.NumberTraces-1;
            dfn.IdTraces.push_back(id_trac); // id di traccia i-esima è i-1 (indice di traccia nelle strutture associate)
            dfn.FractureTraces.push_back({frac1,frac2});
            dfn.TipsTraces.push_back({Tips1,Tips2});
            Vector3d T1 = P0 +P1*t;
            Vector3d T2 = P0 +P2*t;
            Matrix<double,3,2> estremi_traccia{{T1[0],T2[0]},{T1[1],T2[1]},{T1[2],T2[2]}};
            dfn.VerticesTraces.push_back(estremi_traccia);
            double length = (T1[0]-T2[0])*(T1[0]-T2[0])+(T1[1]-T2[1])*(T1[1]-T2[1])+(T1[2]-T2[2])*(T1[2]-T2[2]);
            dfn.LengthTraces.push_back(length);

            InsertSortedTraces(dfn, frac1, id_trac, Tips1, length);
            InsertSortedTraces(dfn, frac2, id_trac, Tips2, length);

            // NB: controlla se conviene con funzione
        }



    }

}


void PrintTraces(const string& outputFile, DFN& dfn)
{
    ofstream output(outputFile);

    output << "# Number of Traces" << "\n" << dfn.NumberTraces << endl;
    output << "# TraceId; FractureId1; FractureId2; X1; Y1; Z1; X2; Y2; Z2 " << endl;

    string sep = "; " ;
    for (unsigned int i=0; i <dfn.NumberTraces; i++)
        output << dfn.IdTraces[i] << sep << (dfn.FractureTraces[i])[0] << sep << (dfn.FractureTraces[i])[1] << sep <<
            (dfn.VerticesTraces[i])(0,0) << sep << (dfn.VerticesTraces[i])(1,0) << sep << (dfn.VerticesTraces[i])(2,0) << sep <<
            (dfn.VerticesTraces[i])(0,1) << sep << (dfn.VerticesTraces[i])(1,1) << sep << (dfn.VerticesTraces[i])(2,1) << endl ;

    output.close();
}


void PrintSortedFractureTraces(const string& outputFile, DFN& dfn)
{
    ofstream output(outputFile);

    string sep = "; " ;
    for (unsigned int i=0; i <dfn.NumberFractures; i++)
    {
        unsigned int num_tot = dfn.P_Traces[i].size() + dfn.NP_Traces[i].size(); // #tracce= #passanti + #NONpassanti
        unsigned int id_frac = dfn.IdFractures[i];
        output << "# FractureId; NumTraces \n" << id_frac << sep << num_tot << endl;
        output << "# TraceId; Tips; Length" << endl;

        // Stampa fratture passanti
        for (auto it = dfn.P_Traces[id_frac].begin(); it != dfn.P_Traces[id_frac].end();it++)
            output << *(it) << sep << false << sep << dfn.LengthTraces[*(it)] << endl ;

        // Stampa fratture non passanti
        for (auto it = dfn.NP_Traces[id_frac].begin(); it != dfn.NP_Traces[id_frac].end();it++)
            output << *(it) << sep << true << sep << dfn.LengthTraces[*(it)] << endl ;
    }

    output.close();
}


PolygonalMesh calculate_fracture_cuts(Matrix3Xd& frac_vertices, list<unsigned int>& p_traces, list<unsigned int>& np_traces,vector<Matrix<double,3,2>>& traces_extremes, double tol)
{
    PolygonalMesh frac;
    frac.tolerance = max(frac.tolerance, tol);

    frac.NumberCell0D = frac_vertices.cols(); // numero vertici frattura
    frac.NumberCell1D = frac.NumberCell0D; // numero lati frattura
    frac.NumberCell2D = 1; // inizialmente frattura unica cell2D

    unsigned int num_tot_traces = p_traces.size() + np_traces.size();
    // Faccio reserve iniziale vettori cell0D
    frac.IdCell0D.reserve(4*num_tot_traces+frac.NumberCell0D);
    frac.CoordinatesCell0D.reserve(4*num_tot_traces+frac.NumberCell0D);
    // Faccio reserve iniziale vettori cell1D
    frac.IdCell1D.reserve(4*num_tot_traces+frac.NumberCell1D);
    frac.VerticesCell1D.reserve(4*num_tot_traces+frac.NumberCell1D);
    // Faccio reserve iniziale vettori cell2D
    frac.IdCell2D.reserve(4*num_tot_traces);
    frac.VerticesCell2D.reserve(4*num_tot_traces);
    frac.VerticesCell2D.reserve(4*num_tot_traces);

    list<unsigned int> external_edges={}; // lista contenente id lati esterni (da escludere per intersezione tracce con lati interni)
    list<unsigned int> internal_edges={}; // lista contenente id lati interni

    // Inserisco frattura come prima cell2D
    frac.IdCell2D.push_back(0);
    frac.VerticesCell2D.push_back({});
    frac.EdgesCell2D.push_back({});

    // Inserisco i vertici frattura --> cell0D iniziali
    for (unsigned int i=0; i<frac.NumberCell0D; i++)
    {
        frac.IdCell0D.push_back(i);
        frac.CoordinatesCell0D.push_back(frac_vertices(all,i));
        frac.VerticesCell2D[0].push_back(i); // verso antiorario
    }
    // Inserisco i lati frattura --> cell1D iniziali
    for (unsigned int i=0; i<frac.NumberCell1D; i++)
    {
        frac.IdCell1D.push_back(i);
        if (i==frac.NumberCell1D-1)
            frac.VerticesCell1D.push_back({i,0}); //ultimo lato unisce primo e ultimo vertice
        else
            frac.VerticesCell1D.push_back({i,i+1});
        frac.EdgesCell2D[0].push_back(i);
        external_edges.push_back(i); // tutti i lati della frattura sono esterni
    }

    // CICLO SU TRACCE PASSANTI
    for (auto it_tr = p_traces.begin(); it_tr != p_traces.end();it_tr++)
    {
        unsigned int id_tr = *(it_tr); // identificatore traccia
        Vector3d ext1_tr = traces_extremes[id_tr](all,0); // coordinate estremo 1 traccia
        Vector3d ext2_tr = traces_extremes[id_tr](all,1); // coordinate estremo 2 traccia
        // Associo ai due estremi i due lati su cui giacciono
        unsigned int l1 = edge_to_traceExtreme(ext1_tr,frac);
        unsigned int l2 = edge_to_traceExtreme(ext2_tr,frac);
        // CASO LIBRO l1=l2

        // Calcolo intersezioni traccia corrente con lati interni
        list<Vector2d> intersezioni= IntersectTraceWithInternalEdges(ext1_tr,ext2_tr,frac,internal_edges);

        // generalizza se passante aggiungi estremi traccia, se non lo è funzione che li calcola
        intersezioni.push_front({l1,0});
        intersezioni.push_back({l2,1});

        // Tagli lungo traccia
        unsigned int counter_intsz=1; //utile per alcuni controlli
        unsigned int sz_int= intersezioni.size(); // per evitare di ciclare quando ho solo più un intersezione nella lista

        unsigned int id_intersez_prec; // per evitare di salvare nuove cell0D di intersezioni interne due volte
        unsigned int id_NEW_V;
        unsigned int id_OLD_V;
        unsigned int id_NEW_E;
        unsigned int id_NEW_E_T;
        unsigned int id_NEW_C;

        for (auto it_intsz = intersezioni.begin(); it_intsz != intersezioni.end();it_intsz++)
        {
            if (counter_intsz<sz_int) // quando mi trovo su ultimo elemento ho finito
            {
                // TAGLIO DA INTERSEZIONE CORRENTE A SUCCESSIVA
                Vector2d intsz1=*(it_intsz);
                Vector2d intsz2=*(next(it_intsz)); // considera elemento successivo in lista intersezioni

                // Cerco cella2D che presenta entrambe i lati coinvolti (cell_2D)
                unsigned int edge1= static_cast<unsigned int>(std::round(intsz1[0])) ; // riporto ad unsigned int per non avere problemi con double
                unsigned int edge2= static_cast<unsigned int>(std::round(intsz2[0])) ; // riporto ad unsigned int per non avere problemi con double
                bool found_cell2D = false;
                unsigned int cell_2D; // fare controllo se non la trova
                for (unsigned int c2D=0; c2D < frac.NumberCell2D; c2D++)
                {
                    if (not found_cell2D)
                    {
                        bool found_l1 = false;
                        for (auto it_L_2D = frac.EdgesCell2D[c2D].begin(); it_L_2D != frac.EdgesCell2D[c2D].end();it_L_2D++)
                        {
                            unsigned int current_edge=*(it_L_2D); // lo salvo in double per compararlo senza problemi con lato salvato come double nel vettore dentro intersezioni
                            if (current_edge== edge1)
                                found_l1 = true; // nella cella corrente è presente edge1
                            if (found_l1)
                            {
                                for (auto it_L2_2D = frac.EdgesCell2D[c2D].begin(); it_L2_2D != frac.EdgesCell2D[c2D].end();it_L2_2D++)
                                {
                                    unsigned int current_edge2=*(it_L2_2D);
                                    if (current_edge2== edge2)
                                    {
                                        found_cell2D = true; // nella cella corrente è presente anche edge2 (è cella cercata)
                                        cell_2D = c2D;
                                        break; // smetto di ciclare sui lati di questa cella per trovare edge2
                                    }
                                }
                                break; // smetto di ciclare sui lati di questa cella per trovare edge1
                            }
                        }
                    }
                }
                if (not found_cell2D)
                    cerr << "ERROR while looking for cell2d to cut"<< endl;
                else
                {
                    // EFFETTUO TAGLIO su cella2D trovata (cell_2D) congiungente intersezione corrente e successiva
                    Vector3d point1 = ext1_tr +(ext2_tr-ext1_tr)*intsz1[1]; // coordinate punto di intersezione corrente (uso ascissa su retta associata alla traccia)
                    Vector3d point2 = ext1_tr +(ext2_tr-ext1_tr)*intsz2[1]; // coordinate punto di intersezione successiva (uso ascissa su retta associata alla traccia)

                    bool seconda_cella = false; // true se sto lavorando su nuova cella
                    bool l1_found = false; // true se ho già trovato l1
                    bool l2_found = false; // true se ho già trovato l2
                    list<unsigned int> lati_iterazione = frac.EdgesCell2D[cell_2D]; //copia lati cella su cui iterare
                    list<unsigned int> vertici_iterazione = frac.VerticesCell2D[cell_2D]; //copia vertici cella su cui iterare

                    auto iter_ver = vertici_iterazione.begin();
                    unsigned int iter_old_cell=1;
                    //NB: STO ESCLUDENDO IL CASO DI DUE VERTICI SU STESSO LATO--- DA FARE
                    // NB. VALUTARE COSA SUCCEDE CON LIBRO
                    for (auto iter_edge = lati_iterazione.begin(); iter_edge != lati_iterazione.end(); iter_edge++)
                    {
                        unsigned int edge=*(iter_edge); // lato corrente
                        unsigned int vert=*(iter_ver); // vertice corrente
                        if (not seconda_cella) // sto lavorando su cella vecchia
                        {
                            if ((edge==edge1) || (edge==edge2) ) // lato corrente è uno dei dei lati con estremo traccia
                            {
                                seconda_cella = true;
                                frac.EdgesCell2D[cell_2D].resize(iter_old_cell); // Cancello dalla lista dei lati della vecchia cella tutti quelli successivi a lato corrente
                                frac.VerticesCell2D[cell_2D].resize(iter_old_cell); // Cancello dalla lista dei vertici della vecchia cella tutti quelli successivi a vertice corrente

                                if (edge==edge1) //lato corrente=lato1
                                {
                                    l1_found=true;
                                    if (counter_intsz==1) // primo taglio (anche il primo estremo non è ancora inserito tra le cell0D)
                                        id_NEW_V = NewCell0D(frac,point1); // creo nuova cell0D per interesezione1
                                    else
                                        id_NEW_V = id_intersez_prec; // già salvata in taglio precedente
                                }
                                else //lato corrente=lato2
                                {
                                    l2_found=true;
                                    id_NEW_V = NewCell0D(frac,point2); // creo nuova cell0D per interesezione2
                                }

                                frac.VerticesCell2D[cell_2D].push_back(id_NEW_V); // aggiungo nuovo vertice a lista di cella originaria

                                if (vert == frac.VerticesCell1D[edge][0]) // vertice corrente è il primo del lato corrente
                                {
                                    id_OLD_V = frac.VerticesCell1D[edge][1];
                                    frac.VerticesCell1D[edge][1] = id_NEW_V;
                                }
                                else // vert == frac.VerticesCell1D[edge][1] // vertice corrente è il secondo del lato corrente
                                {
                                    id_OLD_V = frac.VerticesCell1D[edge][0];
                                    frac.VerticesCell1D[edge][0] = id_NEW_V;
                                }
                                id_NEW_E = NewCell1D(frac, id_NEW_V,id_OLD_V); // creo nuovo lato (sottolato di corrente) da estremo traccia a secondo vertice di lato originario
                                //Inserire in extrnals o internals FAREEEEEEEEEEEE

                                // Creo nuova cell2D
                                id_NEW_C = frac.NumberCell2D;
                                frac.NumberCell2D +=1; // aumento numero cell2D
                                frac.IdCell2D.push_back(id_NEW_C); // inserisco il nuovo id nella mesh
                                frac.VerticesCell2D.push_back({id_NEW_V}); // inserisco primo vertice (estremo traccia su lato corrente)
                                frac.EdgesCell2D.push_back({id_NEW_E}); // inserisco primo lato (ottenuto da lato corrente, tratto da estremo traccia su lato a estremo originario lato)

                            }
                            else // lato corrente NON è coinvolto in traccia
                            {
                                if ((l1_found) && (l2_found)) //lavoro su nuova traccia terminato (si deve concludere la vecchia)
                                {
                                    frac.EdgesCell2D[cell_2D].push_back(edge);
                                    frac.VerticesCell2D[cell_2D].push_back(vert);
                                }
                                else // ancora non creata cella nuova
                                {
                                    iter_old_cell +=1;
                                }
                            }
                        }
                        else // sto lavorando su cella nuova
                        {
                            if ((edge==edge1) || (edge==edge2) ) // lato corrente è uno dei dei lati con estremo traccia
                            {
                                seconda_cella = false;
                                if (l2_found) // ho trovato edge1 (perchè edge2 già trovato)
                                {
                                    l1_found = true;
                                    if (counter_intsz==1) // primo taglio (anche il primo estremo non è ancora inserito tra le cell0D)
                                        id_NEW_V = NewCell0D(frac,point1); // creo nuova cell0D per interesezione1
                                    else
                                        id_NEW_V = id_intersez_prec; // già salvata in taglio precedente
                                }
                                else // ho trovato edge2 (perchè edge1 già trovato)
                                {
                                    l2_found = true;
                                    id_NEW_V = NewCell0D(frac,point2); // creo nuova cell0D per interesezione2
                                }

                                frac.VerticesCell2D[id_NEW_C].push_back(vert); // inserisco vertice in lista di vertici della nuova cella
                                frac.EdgesCell2D[id_NEW_C].push_back(edge); // inserisco lato in lista di lati della nuova cella
                                frac.VerticesCell2D[id_NEW_C].push_back(id_NEW_V); // inserisco nuovo vertice (secondo estremo traccia) in lista di vertici della nuova cella

                                if (vert == frac.VerticesCell1D[edge][0]) // vertice corrente è il primo del lato corrente
                                {
                                    id_OLD_V = frac.VerticesCell1D[edge][1];
                                    frac.VerticesCell1D[edge][1] = id_NEW_V;
                                }
                                else // vert == frac.VerticesCell1D[edge][1] // vertice corrente è il secondo del lato corrente
                                {
                                    id_OLD_V = frac.VerticesCell1D[edge][0];
                                    frac.VerticesCell1D[edge][0] = id_NEW_V;
                                }
                                id_NEW_E = NewCell1D(frac, id_NEW_V,id_OLD_V);// creo nuovo lato (sottolato di corrente) da estremo traccia a secondo vertice di lato originario
                                //Inserire in extrnals o internals FAREEEEEEEEEEEE

                                // Creo lato traccia
                                unsigned int id_other_extr_tr= id_NEW_V-1; //id dell'altra cell0D estremo della traccia (quella già visitata)
                                id_NEW_E_T = NewCell1D(frac, id_NEW_V,id_other_extr_tr);
                                //Inserire in extrnals o internals FAREEEEEEEEEEEE

                                frac.EdgesCell2D[id_NEW_C].push_back(id_NEW_E_T); // inserisco lato traccia nella lista dei lati della cella nuova (e così la concludo)
                                frac.VerticesCell2D[cell_2D].push_back(id_NEW_V); // inserisco secondo estremo traccia in lista vertici cella originaria
                                frac.EdgesCell2D[cell_2D].push_back(id_NEW_E_T); // inserisco lato traccia nella lista dei lati della cella originaria
                                frac.EdgesCell2D[cell_2D].push_back(id_NEW_E); // inserisco nuovo lato (ottenuto a partire da lato corrente) nella lista dei lati della cella originaria

                            }
                            else // lato corrente NON è coinvolto in traccia
                            {
                                frac.VerticesCell2D[id_NEW_C].push_back(vert); // inserisco vertice in lista di vertici della nuova cella
                                frac.EdgesCell2D[id_NEW_C].push_back(edge); // inserisco lato in lista di lati della nuova cella
                            }
                        }
                        iter_ver++;  // incremento iteratore vertici
                    }
                }
                id_intersez_prec = max(frac.VerticesCell1D[edge2][0],frac.VerticesCell1D[edge2][1]); //id di seconda intersezione come cella0D (è l'id maggiore tra i due vertici del lato2
            } // fine sottotaglio

            counter_intsz +=1;
        }

    }


}


}



