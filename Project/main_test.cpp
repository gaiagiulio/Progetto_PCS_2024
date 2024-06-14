#include "TestDFN.hpp"
#include "UCD_test.hpp"
#include "Utils.hpp"
using namespace DFNLibrary;
using namespace std;

int main(int argc, char **argv)

{
/*    DFN dfn;
    DFN_functions fun_dfn;

    // Carica i dati delle fratture
    bool operazione1 = fun_dfn.ImportFractures("C:/Users/Giulia/Desktop/progetto pcs4/Progetto_PCS_2024/Debug/DFN/testDFN.txt", dfn);
    //Calcolo tracce
    fun_dfn.calculateTraces(dfn);
    // Verifica correttezza PrintTraces
    string outputFile1 = "outputTest.txt";
    // Chiama la funzione PrintTraces per stampare i dati nel file
    fun_dfn.PrintTraces(outputFile1, dfn);

    //Verifica correttezza PrintSortedFractureTraces
    string outputFile2 = "outputTest2.txt";
    // Chiama la funzione PrintSortedFractureTraces per stampare i dati nel file
    fun_dfn.PrintSortedFractureTraces(outputFile2, dfn);

    //PolygonalMesh frac;

    // //Verifica correttezza CreateMeshFiles
    // list<unsigned int> external_edges({});
    // vector<Vector2i> edge_to_cells({});
    // vector<list<unsigned int>> ver_to_cells({});
    // Matrix3Xd frac_vertices(3,3);
    // frac_vertices << 0, 7, 0,
    //                 0, 0, 5,
    //                 0, 0, 0;
    // list<unsigned int> p_traces({});
    // list<unsigned int> np_traces({});
    // fun_dfn.InitializeMesh(frac, external_edges, edge_to_cells, ver_to_cells, frac_vertices, p_traces, np_traces);
    // CreateMeshFiles(frac, 3);    //fun_dfn.CreateMeshFiles(frac, id_frac);
*/
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();

    return 0;

}
