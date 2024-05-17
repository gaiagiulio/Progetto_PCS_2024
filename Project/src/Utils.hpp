#ifndef DFN_UTILS_H
#define DFN_UTILS_H

#include "DiscreteFractureNetwork.hpp"

namespace DFNLibrary{

/** Import fractures from filepath (path + file name) ( and test if it's correct__ AGGIUNGI)
 *  dfn: a DFN struct
 *  return the result of the reading (true if is success, false otherwise) **/
bool ImportFractures(const string& filepath, DFN& dfn);


/** AGGIUNGI DESCRIZIONE PER CALCOLO TRACCE E COSA FA **/
void calculateTraces(DFN& dfn);

}

#endif
