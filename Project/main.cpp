#include "DiscreteFractureNetwork.hpp"
#include "Utils.hpp"


using namespace DFNLibrary;

int main()
{
    DFN dfn;
    // FAI INSERIRE TOLLERANZA
    // fai inserire a utente quale file leggere

    bool operazione = ImportFractures("D:/Gaia/Politecnico/Programmazione e calcolo scientifico/Progetto_PCS_2024/Project/DFN/FR10_data.txt", dfn);
    std::cout << "N fratture: " << dfn.NumberFractures << std::endl;

    calculateTraces(dfn);

    PrintTraces("PROVA_output_1.txt", dfn);
    PrintSortedFractureTraces("PROVA_output_2.txt", dfn);


    /**std::list<int> l={1,2,3,4,5,6};
    for (auto it = l.begin(); it != l.end();it++)
    {
        std::cout << *(it) << std::endl;
    }**/
    /**for (int n : l)
    {
        std::cout << n << std::endl;
    }; **/
    return 0;
}
