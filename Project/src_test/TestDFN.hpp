#ifndef __unit_test_hpp
#define __unit_test_hpp
#include "Utils.hpp"
#include "DiscreteFractureNetwork.hpp"
#include <gtest/gtest.h>
using namespace DFNLibrary;
DFN_functions fun_dfn;
DFN dfn;
<<<<<<< Updated upstream

TEST(ImportFractures_Test,ImportFractures)  //verifico che l'apertura del file di input per importare le fratture vada a buon fine
{
    bool expected=true;
    bool result=fun_dfn.ImportFractures("./test.txt",dfn);
    ASSERT_TRUE(expected=result);
}

TEST(NormalToPlane_Test,Plane) //testo la funzione che calcola la normale al piano
=======
/**TEST(NormalToPlane_Test,Plane)
>>>>>>> Stashed changes
{
    Vector3d p0(6.7949650570084286e-01,5.1566886122886846e-01,1.9054542365205804e-01);
    Vector3d p1(2.0959413133064569e-01,9.9389350435296486e-01,1.9054542365205804e-01);
    Vector3d p2(7.7027229455623514e-02,8.6363358811981283e-01,5.9704177318758545e-01);
    Vector3d expected(1.943965716878755e-01,1.910135998195619e-01,1.246062032624707e-01);
    Vector3d result=fun_dfn.NormalToPlane(p0,p1,p2);
    ASSERT_TRUE(result.isApprox(expected,1e-10));
}

TEST(Acissa_Test,Ascissa)  //testo la funzione che calcola l'ascissa curvilinea
{
    Vector3d V_P0(2,2,4);
    Vector3d Ve_P0(0,0,0);
    Vector3d t1(1,1,2);
    Vector3d t2(-1,-1,-2);
    double expected1(2);
    double expected2(-2);
    double expected3(0);
    bool operazione = fun_dfn.ImportFractures("./test.txt", dfn);
    double s1=fun_dfn.ascissa_curvilinea(V_P0,t1);  //ascissa curvilinea positiva
    double s2=fun_dfn.ascissa_curvilinea(V_P0,t2);  //ascissa curvilinea negativa
    double s3=fun_dfn.ascissa_curvilinea(Ve_P0,t1);  //ascissa curvilinea nulla
    ASSERT_EQ(expected1,s1);
    ASSERT_EQ(expected2,s2);
    ASSERT_EQ(expected3,s3);
}

TEST(Intersection_Test, NoIntersection) //verifico il caso in cui la retta non interseca la frattura
{
    Vector3d P0(0,2,0);
    Vector3d t(1,2,0);
    Vector3d n(0,0,1);
    Vector4d expected(0,0,0,1);
    bool operazione = fun_dfn.ImportFractures("./test.txt", dfn);
    Vector4d result=fun_dfn.IntersectionFractureWithLine(dfn,0,P0,t,n);
    ASSERT_TRUE(result.isApprox(expected,1e-10));
}

TEST(Intersection_Test, IntersectionInVertice)  //verifico il caso in cui la retta interseca la frattura in un suo vertice
{
    Vector3d P0(1/2,1/4,0);
    Vector3d t(1/2,1/4,0);
    Vector3d n(0,0,1);
    Vector4d expected(0,0,1,1);
    bool operazione = fun_dfn.ImportFractures("./test.txt", dfn);
    Vector4d result=fun_dfn.IntersectionFractureWithLine(dfn,1,P0,t,n);
    ASSERT_EQ(expected[2],result[2]);
}

TEST(Intersection_Test, IntersectionTwoEdges)  //verifico il caso in cui la retta interseca la frattura su due lati
{
    Vector3d P0(0,0,0);
    Vector3d t(1,1,0);
    Vector3d n(0,0,1);
    Vector4d expected(-1,2,1,1);
    bool operazione = fun_dfn.ImportFractures("./test.txt", dfn);
    Vector4d result=fun_dfn.IntersectionFractureWithLine(dfn,2,P0,t,n);
    ASSERT_TRUE(result.isApprox(expected,1e-10));
}**/

// TEST(Intersection_Test, BookCase)  //verifico il "caso libro", cioè quello in cui la retta passa per un lato della frattura
// {
//     Vector3d P0(0,0,0);
//     Vector3d t(1,1,0);
//     Vector3d n(0,0,1);
//     Vector4d expected(0,2,1,0);
//     bool operazione = fun_dfn.ImportFractures("./test.txt", dfn);
//     Vector4d result=fun_dfn.IntersectionFractureWithLine(dfn,3,P0,t,n);
//     ASSERT_TRUE(result.isApprox(expected,1e-10));
// }

<<<<<<< Updated upstream
TEST(Intersection_Test, IntersectionVerticeEdge)  //verifico il caso in cui la retta interseca la frattura su un vertice e un lato
=======
/**TEST(Intersection_Test, IntersectionVerticeEdge)
>>>>>>> Stashed changes
{
    Vector3d P0(3,1,0);
    Vector3d t(0,1,0);
    Vector3d n(0,0,1);
    Vector4d expected(-1,1,1,1);
    bool operazione = fun_dfn.ImportFractures("./test.txt", dfn);
    Vector4d result=fun_dfn.IntersectionFractureWithLine(dfn,4,P0,t,n);
    ASSERT_TRUE(result.isApprox(expected,1e-10));
}**/



#endif

//.inp

