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
    bool operazione = fun_dfn.ImportFractures("./DFN/FR10_data.txt", dfn);

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

    // Lettura da terminale di id delle fratture per cui creare i file (imposta "Run in terminal" in Projects per poter scrivere in input da terminale)
    /**cout << "Inserisci gli id delle fratture per cui creare i file da esportare su Paraview (andando a capo dopo ogni id e termmina con F): " << endl;
    unsigned int ID;
    list<unsigned int> fracture_to_print={};
    while (cin >> ID )
        fracture_to_print.push_back(ID);

    for (auto it_frac = fracture_to_print.begin(); it_frac != fracture_to_print.end(); it_frac++ )
    {
        unsigned int id_frac= *(it_frac); // id frattura per cui creare i file
        cout << "Creo file per frattura " << id_frac <<endl;
        CreateMeshFiles(cutted_fractures[id_frac],id_frac); // Crea file per la frattura id_frac
    }**/

    return 0;
}
