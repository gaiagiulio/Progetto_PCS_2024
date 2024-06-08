#include "DiscreteFractureNetwork.hpp"
#include "Utils.hpp"


using namespace DFNLibrary;

int main(int argc, char ** argv)
{

    DFN dfn;
    DFN_functions fun_dfn;

    // Lettura di tolleranza da command line
    double tol = 0.0;
    if(argc == 2)
    {
        istringstream str(argv[1]);
        str >> tol;
        cout << "Tolleranza da command line: "<< tol << endl;

        dfn.tolerance = max(dfn.tolerance, tol);
    }

    //Importazione fratture
    bool operazione = fun_dfn.ImportFractures("C:/Users/utente/PROSEGUIMENTO_PROGETTOPCS/Progetto_PCS_2024/Project/DFN/FR10_data.txt", dfn);
    if (not operazione)
    {
        cerr << "Error while importing fractures" << endl;
        return 1;
    }
    std::cout << "N fratture: " << dfn.NumberFractures << std::endl;

    fun_dfn.calculateTraces(dfn); // Calcolo tracce



    // Stampa output tracce su file
    fun_dfn.PrintTraces("PROVA_output_1.txt", dfn);
    fun_dfn.PrintSortedFractureTraces("PROVA_output_2.txt", dfn);

    vector<PolygonalMesh> cutted_fractures(dfn.NumberFractures); //--> elemento i-esimo= PolygonalMesh di frattura con id=i tagliata

    // Taglio fratture con le tracce
    for (unsigned int i=0; i<dfn.NumberFractures; i++)
    {
        unsigned int id_frac = dfn.IdFractures[i];
        // nb: controlla che tol sia ancora tolleranza da command line (e se vuoi cambiala)
        double tol=1.e-10;
        cout << "Lavoro su frattura : "<< id_frac << endl;
        cutted_fractures[i] = fun_dfn.calculate_fracture_cuts(dfn.VerticesFractures[id_frac], dfn.P_Traces[id_frac], dfn.NP_Traces[id_frac], dfn.VerticesTraces, tol);
        cout << "Frattura " << id_frac<< " restituita \n" << endl;
    }

    // Crea file da esportare in Paraview per visualizzare la PolygonalMesh generata con il taglio della frattura
    // Togliere il commentato se si vogliono creare i file
    // Lettura da terminale di id delle fratture per cui creare i file DA FAREEEEEEEEEEEEEEEEEEEEE
    list<unsigned int> fracture_to_print = {0,1,2};
    for (auto it_frac = fracture_to_print.begin(); it_frac != fracture_to_print.end(); it_frac++ )
    {
        unsigned int id_frac= *(it_frac); // id frattura per cui creare i file
        CreateMeshFiles(cutted_fractures[id_frac],id_frac); // Crea file per la frattura id_frac
    }


    return 0;
}
/**list<int> l={1,2,3,4,5};
    auto it= l.begin();
    l.insert(++it,8);
    l.insert(it,9);
    for (auto itl= l.begin(); itl!= l.end(); itl++)
    {
        cout << *(itl) << endl;
    }**/

/**
    std::list<int> l={1,2,3,4,5,6};
    auto it = l.begin();
    while (it != l.end())
    {
        std::cout << *(it) << std::endl;
        it++;
    }
 **/

/**std::list<int> l={1,2,3,4,5,6};
    for (auto it = l.begin(); it != l.end();it++)
    {
        std::cout << *(it) << std::endl;
    }**/
/**std::list<int> l={1,2,3,4,5,6};
    unsigned int j=1;
    unsigned int sz_l= l.size();
    for (auto it = l.begin(); it != l.end();it++)
    {
        if (j<sz_l)
            std::cout << *(it) << " e "<< *(next(it))<<std::endl;
        j +=1;
    }**/
/**for (int n : l)
    {
        std::cout << n << std::endl;
    }; **/
