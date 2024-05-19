#include "Utils.hpp"
#include "DiscreteFractureNetwork.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

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
            if (length > dfn.LengthTraces[*(it)])
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
            if (length > dfn.LengthTraces[*(it)])
            {
                dfn.P_Traces[frac].insert(it,id_tr);
                inserted = true;
            }
        }
        if (not inserted)
            dfn.P_Traces[frac].push_back(id_tr);
    }
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

    // Reserve dei vettori associati alle tracce: Scelgo numero grande perchè non so quante saranno.
    const unsigned int N = 5000; // crea funzione sensata che stimi una costante da cui partire in base al numero di fratture

    dfn.IdTraces.reserve(N);
    dfn.FractureTraces.reserve(N);
    dfn.TipsTraces.reserve(N);
    dfn.VerticesTraces.reserve(N);
    dfn.LengthTraces.reserve(N);

    dfn.P_Traces.resize(dfn.NumberFractures); // NB: check se serve inizializzare con le liste vuote
    dfn.NP_Traces.resize(dfn.NumberFractures);

    for (unsigned int i=0; i < dfn.NumberFractures; i++)
    {
        unsigned int frac1 = dfn.IdFractures[i] ;
        Matrix3Xd& ver1 = dfn.VerticesFractures[frac1]; // matrice dei vertici della frattura 1
        //Cerco il piano della frattura 1 : calcolo prendendo i primi tre vertici
        Vector3d p10 = ver1(all,0);
        Vector3d p11 = ver1(all,1);
        Vector3d p12 = ver1(all,2);

        Vector3d n1 = NormalToPlane(p10,p11,p12); // vettore normale a piano 1
        double d1=n1.dot(p10); // piano 1 x: dot(n1,x)=d1

        for (unsigned int j=0; (j < dfn.NumberFractures) && (j>i); j++)
        {
            unsigned int frac2 = dfn.IdFractures[j];
            Matrix3Xd& ver2 = dfn.VerticesFractures[frac2]; // matrice dei vertici della frattura 2

            // CONTROLLO PER ESCLUSIONE CASI

            //Cerco il piano della frattura 2 : calcolo prendendo i primi tre vertici
            Vector3d p20 = ver2(all,0);
            Vector3d p21 = ver2(all,1);
            Vector3d p22 = ver2(all,2);

            Vector3d n2 = NormalToPlane(p20,p21,p22); // vettore normale a piano 2
            double d2=n2.dot(p20); // piano 2 x: dot(n2,x)=d2

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
                break; // passo a considerare altra coppia di fratture se una non interseca la retta
            //Frattura 2
            Vector3d int2= IntersectionFractureWithLine(dfn, frac2, P0, t, n2); // vettore con (q1,q2,flag) con q1 e q2 ascisse di intersezioni se esistono
            if (int2[2]<numeric_limits<double>::epsilon()) // uso tol con epsiolon di macchina perchè è in vettore di double (anche se non dovrei avere problemi)
                break;

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
                    break;
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
                                break;
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
                    break;
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
                                break;
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

            // Inserisci traccia trovata (se non è stata trovata si è incappati in break)
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
            output << *(it) << sep << false << sep << dfn.LengthTraces[*(it)] ;

        // Stampa fratture non passanti
        for (auto it = dfn.NP_Traces[id_frac].begin(); it != dfn.NP_Traces[id_frac].end();it++)
            output << *(it) << sep << true << sep << dfn.LengthTraces[*(it)] ;
    }

    output.close();
}



}



