#ifndef unit_test_hpp
#define unit_test_hpp

#include <gtest/gtest.h>
#include "Utils.hpp"
#include "DiscreteFractureNetwork.hpp"
#include "Eigen/Eigen"
#include <iostream>
using namespace std;
using namespace Eigen;

TEST(NormalToPlane_Test,Plane)
{
    Vector3d p0(6.7949650570084286e-01,5.1566886122886846e-01,1.9054542365205804e-01);
    Vector3d p1(2.0959413133064569e-01,9.9389350435296486e-01,1.9054542365205804e-01);
    Vector3d p2(7.7027229455623514e-02,8.6363358811981283e-01,5.9704177318758545e-01);
    Vector3d expected(1.943965716878755e-01,1.910135998195619e-01,1.246062032624707e-01);
    DFNLibrary::DFN_functions dfn_fun;
    Vector3d result=dfn_fun.NormalToPlane(p0,p1,p2);
    ASSERT_TRUE(result.isApprox(expected, 1e-10));
}

#endif
