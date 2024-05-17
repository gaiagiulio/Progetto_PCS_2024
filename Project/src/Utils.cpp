#include "Utils.hpp"
#include "DiscreteFractureNetwork.hpp"

#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;
using namespace Eigen;

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

    // gestione dei casi dove uiìnico punto di intersezione con r è il primo vertice (A) e dove la seconda intersezione è tra il primo e l'ultimo vertice (B)
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

            // Calcolo tracce



            // se entrambe ok calcola traccia usando ascissa curvilinea
            // classifica
        }
    }



}

}



