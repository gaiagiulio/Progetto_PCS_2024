#include "DiscreteFractureNetwork.hpp"
#include "Utils.hpp"


using namespace DFNLibrary;

int main()
{
    DFN dfn;
    // FAI INSERIRE TOLLERANZA

    bool operazione = ImportFractures("D:/Gaia/Politecnico/Programmazione e calcolo scientifico/Progetto_PCS_2024/Project/DFN/FR3_data.txt", dfn);
    std::cout << operazione << std::endl;
    return 0;
}
