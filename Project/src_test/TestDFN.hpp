#ifndef __unit_test_hpp
#define __unit_test_hpp
#include "Utils.hpp"
#include "DiscreteFractureNetwork.hpp"
#include <gtest/gtest.h>
using namespace DFNLibrary;
DFN_functions fun_dfn;
DFN dfn;
TEST(NormalToPlane_Test,Plane)
{
    Vector3d p0(6.7949650570084286e-01,5.1566886122886846e-01,1.9054542365205804e-01);
    Vector3d p1(2.0959413133064569e-01,9.9389350435296486e-01,1.9054542365205804e-01);
    Vector3d p2(7.7027229455623514e-02,8.6363358811981283e-01,5.9704177318758545e-01);
    Vector3d expected(1.943965716878755e-01,1.910135998195619e-01,1.246062032624707e-01);
    Vector3d result=fun_dfn.NormalToPlane(p0,p1,p2);
    ASSERT_TRUE(result.isApprox(expected,1e-10));
}

TEST(Intersection_Test, NoIntersection)
{
    Vector3d P0(0,2,0);
    Vector3d t(1,2,0);
    Vector3d n(0,0,1);
    Vector3d expected(0,0,0);
    bool operazione = fun_dfn.ImportFractures("./test.txt", dfn);
    Vector3d result=fun_dfn.IntersectionFractureWithLine(dfn,0,P0,t,n);
    ASSERT_TRUE(result.isApprox(expected,1e-10));
}

TEST(Intersection_Test, IntersectionInPoint)
{
    Vector3d P0(1/2,1/4,0);
    Vector3d t(1/2,1/4,0);
    Vector3d n(0,0,1);
    Vector3d expected(2,0,1);
    bool operazione = fun_dfn.ImportFractures("./test.txt", dfn);
    Vector3d result=fun_dfn.IntersectionFractureWithLine(dfn,1,P0,t,n);
    ASSERT_TRUE(result.isApprox(expected,1e-10));
}

// TEST(Intersection_Test, Intersection)
// {
//     Vector3d P0(-1,-1.33,0);
//     Vector3d t(1,2.67,0);
//     Vector3d n(0,0,1);
//     Vector3d expected(-1,2,1);
//     bool operazione = fun_dfn.ImportFractures("./test.txt", dfn);
//     Vector3d result=fun_dfn.IntersectionFractureWithLine(dfn,2,P0,t,n);
//     ASSERT_TRUE(result.isApprox(expected,1e-10));
// }
#endif

//.inp

