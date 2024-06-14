#ifndef __unit_test_hpp
#define __unit_test_hpp
#include "Utils.hpp"
#include "DiscreteFractureNetwork.hpp"
#include <gtest/gtest.h>
using namespace DFNLibrary;
DFN_functions fun_dfn;
DFN dfn;
PolygonalMesh frac;

TEST(ImportFractures_Test,ImportFractures)  //verifico che l'apertura del file di input per importare le fratture vada a buon fine
{
    bool expected=true;
    bool result=fun_dfn.ImportFractures("C:/Users/utente/PROJECT DEFINITIVO/Progetto_PCS_2024/Debug/test.txt",dfn);
    ASSERT_TRUE(expected=result);
}

TEST(NormalToPlane_Test,Plane) //testo la funzione che calcola la normale al piano
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
}

TEST(Intersection_Test, BookCase_firstVectice)  //verifico il "caso libro", cioè quello in cui la retta passa per un lato della frattura e il primo vertice si trova sulla retta
{
    Vector3d P0(0,0,0);
    Vector3d t(1,1,0);
    Vector3d n(0,0,1);
    Vector4d expected(0,2,1,0);
    bool operazione = fun_dfn.ImportFractures("./test.txt", dfn);
    Vector4d result=fun_dfn.IntersectionFractureWithLine(dfn,3,P0,t,n);
    ASSERT_TRUE(result.isApprox(expected,1e-10));
}

TEST(Intersection_Test, BookCase)  //verifico il "caso libro", cioè quello in cui la retta passa per un lato della frattura ma il primo vertice non si trova sulla retta
{
    Vector3d P0(0,-1,0);
    Vector3d t(1,1,0);
    Vector3d n(0,0,1);
    Vector4d expected(0,1,1,0);
    bool operazione = fun_dfn.ImportFractures("./test.txt", dfn);
    Vector4d result=fun_dfn.IntersectionFractureWithLine(dfn,4,P0,t,n);
    ASSERT_TRUE(result.isApprox(expected,1e-10));
}


TEST(Intersection_Test, IntersectionVerticeEdge)  //verifico il caso in cui la retta interseca la frattura in un suo lato e un vertice
{
    Vector3d P0(3,1,0);
    Vector3d t(0,1,0);
    Vector3d n(0,0,1);
    Vector4d expected(-1,1,1,1);
    bool operazione = fun_dfn.ImportFractures("./test.txt", dfn);
    Vector4d result=fun_dfn.IntersectionFractureWithLine(dfn,5,P0,t,n);
    ASSERT_TRUE(result.isApprox(expected,1e-10));
}

TEST(IntersectionLines_Test,IntersectionLines)
{
    Vector3d t1(1,1,0);
    Vector3d t2(-1,1,0);
    Vector3d p1(0,0,0);
    Vector3d p2(2,0,0);
    double tol=1e-10;
    Vector3d expected(1,1,1);
    Vector3d result=fun_dfn.IntersectionBetweenLines(t1,t2,p1,p2,tol);
}

TEST(CalculateTraces_Test, CalculateTraces)
{
    DFN dfn;
    bool operazione=fun_dfn.ImportFractures("C:/Users/utente/PROJECT DEFINITIVO/Progetto_PCS_2024/Debug/test2.txt",dfn);
    fun_dfn.calculateTraces(dfn);
    vector<unsigned int> IdTraces={{0},{1}}; // Identificatori tracce --> intero positivo (dimensione 1)
    vector<Vector2i> FractureTraces = {{0,1},{0,2}}; // Fratture associate a traccia
    vector<Vector<bool,2>> TipsTraces = {{false,true},{true,true}}; // Tips booleano false= passante, true= non passante
    vector<Matrix<double,3,2>> VerticesTraces; // vettore con estremi traccia
    Matrix<double,3,2> mat1;
    mat1 << 4, 4,
            5, 0,
            0, 0;
    Matrix<double,3,2> mat2;
    mat2 << 1, 1,
            2, 0,
            0, 0;
    VerticesTraces.push_back(mat1);
    VerticesTraces.push_back(mat2);
    vector<double> LengthTraces = {25,4}; // vettore con lunghezza tracce
    vector<list<unsigned int>> P_Traces = {{0}}; // Lista di identificatori tracce passanti di frattura ordinati per lunghezza
    vector<list<unsigned int>> NP_Traces = {{1}}; // Lista di identificatori tracce NON passanti di frattura ordinati per lunghezza
    ASSERT_EQ(IdTraces,dfn.IdTraces);
    ASSERT_EQ(FractureTraces,dfn.FractureTraces);
    ASSERT_EQ(TipsTraces,dfn.TipsTraces);
    ASSERT_EQ(VerticesTraces,dfn.VerticesTraces);
    ASSERT_EQ(LengthTraces,dfn.LengthTraces);
}

TEST(InsertSortedTraces_test,InsertSortedTraces)
{
    DFN dfn;
    bool result=fun_dfn.ImportFractures("C:/Users/utente/PROJECT DEFINITIVO/Progetto_PCS_2024/Debug/test3.txt",dfn);
    unsigned int frac=0;
    unsigned int id_tr=0;
    bool Tips=false;
    double length=3.51;
    fun_dfn.InsertSortedTraces(dfn,frac,id_tr,Tips,length);
    vector<list<unsigned int>> expected={{0}};
    ASSERT_EQ(expected,dfn.P_Traces);
}


//TEST PART 2
TEST(NewCell0D_test,NewCell0D) //verifico la creazione di una nuova cella 0D
{
    list<unsigned int> external_edges({});
    vector<Vector2i> edge_to_cells({});
    vector<list<unsigned int>> ver_to_cells({});
    Matrix3Xd frac_vertices(3,3);
    frac_vertices << 0, 7, 0,
        0, 0, 5,
        0, 0, 0;
    list<unsigned int> p_traces({});
    list<unsigned int> np_traces({});
    fun_dfn.InitializeMesh(frac, external_edges, edge_to_cells, ver_to_cells, frac_vertices, p_traces, np_traces);
    Vector3d point(7,5,0);
    unsigned int expected=3;
    unsigned int result=fun_dfn.NewCell0D(frac, point, ver_to_cells);
    ASSERT_EQ(expected,result);
}

TEST(NewCell1D_test,NewCell1D) //verifico la creazione di una nuova cella 1D
{
    list<unsigned int> external_edges({});
    vector<Vector2i> edge_to_cells({});
    vector<list<unsigned int>> ver_to_cells({});
    Matrix3Xd frac_vertices(3,3);
    frac_vertices << 0, 7, 0,
        0, 0, 5,
        0, 0, 0;
    list<unsigned int> p_traces({});
    list<unsigned int> np_traces({});
    fun_dfn.InitializeMesh(frac, external_edges, edge_to_cells, ver_to_cells, frac_vertices, p_traces, np_traces);
    unsigned int ver1=1;
    unsigned int ver2=2;
    unsigned int expected=3;
    unsigned int result=fun_dfn.NewCell1D(frac, ver1, ver2, edge_to_cells);
    ASSERT_EQ(expected,result);
}

TEST(InitializeMesh_test,InitializeMesh)
{
    PolygonalMesh frac;
    list<unsigned int> external_edges={0,1,2,3,4,5};
    vector<Vector2i> edge_to_cells={{0,-1},{1,-1},{1,-1},{1,-1},{0,-1},{0,-1},{0,1}};
    vector<list<unsigned int>> ver_to_cells={{0},{0,1},{1},{1},{0,1},{0}};
    Matrix3Xd frac_vertices(3,4);
    frac_vertices << 1, 3, 1, -2,
                    0, 0, 2, 0,
                    0, 2, 1, 1;
    list<unsigned int> p_traces={0};
    list<unsigned int> np_traces={};
    fun_dfn.InitializeMesh(frac,external_edges,edge_to_cells,ver_to_cells,frac_vertices,p_traces,np_traces);
    unsigned int NumberCell0D=4;
    vector<unsigned int> IdCell0D = {0,1,2,3};
    vector<Vector3d> CoordinatesCell0D = {{1,0,0},{3,0,2},{1,2,1},{-2,0,1}};

    unsigned int NumberCell1D = 4;
    vector<unsigned int> IdCell1D = {0,1,2,3};
    vector<Vector2i> VerticesCell1D = {{0,1},{1,2},{2,3},{3,0}};

    unsigned int NumberCell2D = 1;
    vector<unsigned int> IdCell2D = {0};
    vector<list<unsigned int>> VerticesCell2D = {{0,1,2,3}};
    vector<list<unsigned int>> EdgesCell2D = {{0,1,2,3}};

    double tolerance = 10000*numeric_limits<double>::epsilon();
    ASSERT_EQ(NumberCell0D,frac.NumberCell0D);
    ASSERT_EQ(IdCell0D,frac.IdCell0D);
    ASSERT_EQ(CoordinatesCell0D,frac.CoordinatesCell0D);
    ASSERT_EQ(NumberCell1D,frac.NumberCell1D);
    ASSERT_EQ(IdCell1D,frac.IdCell1D);
    ASSERT_EQ(VerticesCell1D,frac.VerticesCell1D);
    ASSERT_EQ(NumberCell2D,frac.NumberCell2D);
    ASSERT_EQ(IdCell2D,frac.IdCell2D);
    ASSERT_EQ(VerticesCell2D,frac.VerticesCell2D);
    ASSERT_EQ(EdgesCell2D,frac.EdgesCell2D);
}


TEST(Edge_test,Edge_to_traceExtreme) //verifico che l'estremo della traccia(in questo caso passante) sia associato ad un lato esterno della frattura
{
    list<unsigned int> external_edges({});
    vector<Vector2i> edge_to_cells({});
    vector<list<unsigned int>> ver_to_cells({});
    Matrix3Xd frac_vertices(3,3);
    frac_vertices << 0, 7, 0,
        0, 0, 5,
        0, 0, 0;
    list<unsigned int> p_traces({});
    list<unsigned int> np_traces({});
    fun_dfn.InitializeMesh(frac, external_edges, edge_to_cells, ver_to_cells, frac_vertices, p_traces, np_traces);
    Vector3d ext_tr(4,0,0);
    Vector2i expected(0,-1);
    Vector2i result(fun_dfn.edge_to_traceExtreme(ext_tr,external_edges,frac));
    ASSERT_TRUE(result.isApprox(expected,1e-10));
}

TEST(Edge_test,Edge_to_traceExtreme_notFound) //verifico che l'estremo della traccia non trovi nessuna cella1D esterna
{
    list<unsigned int> external_edges({});
    vector<Vector2i> edge_to_cells({});
    vector<list<unsigned int>> ver_to_cells({});
    Matrix3Xd frac_vertices(3,3);
    frac_vertices << 0, 7, 0,
        0, 0, 5,
        0, 0, 0;
    list<unsigned int> p_traces({});
    list<unsigned int> np_traces({});
    fun_dfn.InitializeMesh(frac, external_edges, edge_to_cells, ver_to_cells, frac_vertices, p_traces, np_traces);
    Vector3d ext_tr(4,2,0);
    Vector2i expected(-1,-1);
    Vector2i result(fun_dfn.edge_to_traceExtreme(ext_tr,external_edges,frac));
    ASSERT_TRUE(result.isApprox(expected,1e-10));
}

TEST(Edge_test,Edge_to_traceExtreme_cell0D) //verifico che l'estremo della traccia coincide con una cell0D esistente
{
    list<unsigned int> external_edges({});
    vector<Vector2i> edge_to_cells({});
    vector<list<unsigned int>> ver_to_cells({});
    Matrix3Xd frac_vertices(3,3);
    frac_vertices << 0, 7, 0,
        0, 0, 5,
        0, 0, 0;
    list<unsigned int> p_traces({});
    list<unsigned int> np_traces({});
    fun_dfn.InitializeMesh(frac, external_edges, edge_to_cells, ver_to_cells, frac_vertices, p_traces, np_traces);
    Vector3d ext_tr(0,0,0);
    Vector2i expected(-2,0);
    Vector2i result(fun_dfn.edge_to_traceExtreme(ext_tr,external_edges,frac));
    ASSERT_TRUE(result.isApprox(expected,1e-10));
}

TEST(IntersectCellsEdges_test,IntersectCellsEdges)
{
    PolygonalMesh frac;
    frac.NumberCell0D = 8; // Numero totale Cell0D
    for (unsigned int i=0; i<8; i++)
        frac.IdCell0D.push_back(i);
    frac.CoordinatesCell0D.push_back({0,0,0});
    frac.CoordinatesCell0D.push_back({2,0,0});
    frac.CoordinatesCell0D.push_back({4,0,0});
    frac.CoordinatesCell0D.push_back({7,0,0});
    frac.CoordinatesCell0D.push_back({7,5,0});
    frac.CoordinatesCell0D.push_back({4, 5 ,0});
    frac.CoordinatesCell0D.push_back({2,5,0});
    frac.CoordinatesCell0D.push_back({0,5,0});
    frac.NumberCell1D = 10; // Numero totale Cell1D
    for (unsigned int i=0; i<10; i++)
        frac.IdCell1D.push_back((i));
    frac.VerticesCell1D.push_back({0,1});
    frac.VerticesCell1D.push_back({1,2});
    frac.VerticesCell1D.push_back({2,3});
    frac.VerticesCell1D.push_back({3,4});
    frac.VerticesCell1D.push_back({4,5});
    frac.VerticesCell1D.push_back({5,6});
    frac.VerticesCell1D.push_back({6,7});
    frac.VerticesCell1D.push_back({7,0});
    frac.VerticesCell1D.push_back({6,1});
    frac.VerticesCell1D.push_back({5,2});
    frac.NumberCell2D = 3; // Numero totale Cell2D
    for (unsigned int i=0; i<3; i++)
        frac.IdCell2D.push_back(i);
    frac.VerticesCell2D.push_back({0,1,6,7});
    frac.VerticesCell2D.push_back({1,2,5,6});
    frac.VerticesCell2D.push_back({2,3,4,5});
    frac.EdgesCell2D.push_back({0,8,6,7});
    frac.EdgesCell2D.push_back({1,9,5,8});
    frac.EdgesCell2D.push_back({2,3,4,9});
    // lati Cell2D --> lista di lati in senso antiorario
    double tolerance = 10000*numeric_limits<double>::epsilon();
    list<unsigned int> internal_edges({8,9});
    Vector3d t_T(4,0,0);
    Vector3d ext1_tr(0,2,0);
    unsigned int c2D=0;
    bool edge_found0=true;
    unsigned int l10=7;
    unsigned int l20=7;
    unsigned int v0=0;
    double s0=0.0;
    Vector<double,5> expected(8,0.5,0,true,true);
    Vector<double,5> result(fun_dfn.IntersectCellEdges(frac, internal_edges, t_T, ext1_tr, c2D, edge_found0, l10, l20, v0, s0));
    ASSERT_TRUE(result.isApprox(expected,1e-10));
}

TEST(BookCase_test,BookCase_EE)
{
    PolygonalMesh frac;
    frac.NumberCell0D = 4; // Numero totale Cell0D
    for (unsigned int i=0; i<4; i++)
        frac.IdCell0D.push_back(i); // identificatori Cell0D --> intero positivo (dimensione 1)
    frac.CoordinatesCell0D.push_back({4,0,0});
    frac.CoordinatesCell0D.push_back({4,3,0});
    frac.CoordinatesCell0D.push_back({0,3,0});
    frac.CoordinatesCell0D.push_back({0,0,0});


    frac.NumberCell1D = 4; // Numero totale Cell1D
    for (unsigned int i=0; i<4; i++)
        frac.IdCell1D.push_back(i); // identificatori Cell1D --> intero positivo (dimensione 1)
    frac.VerticesCell1D.push_back({0,1});
    frac.VerticesCell1D.push_back({1,2});
    frac.VerticesCell1D.push_back({2,3});
    frac.VerticesCell1D.push_back({3,0});

    frac.NumberCell2D = 1; // Numero totale Cell2D
    frac.IdCell2D = {0}; // identificatori Cell2D --> intero positivo (dimensione 1)
    frac.VerticesCell2D.push_back({0,1,2,3});
    frac.EdgesCell2D.push_back({0,1,2,3});

    double tolerance = 10000*numeric_limits<double>::epsilon();

    list<unsigned int> internal_edges={};
    list<unsigned int> external_edges={0,1,2,3};
    unsigned int c=0;
    unsigned int l=3;
    Vector3d ext1_tr(1,0,0);
    Vector3d ext2_tr(3,0,0);
    vector<Vector2i> edge_to_cells={{0,-1}, {0,-1}, {0,-1}, {0,-1}};
    vector<list<unsigned int>> ver_to_cells={{0},{0},{0},{0}};
    list<unsigned int>::iterator it_l={frac.EdgesCell2D[c].begin()};
    list<unsigned int>::iterator it_ver={frac.VerticesCell2D[c].begin()};
    Vector<unsigned int,4> expected(4,5,4,5);
    Vector<unsigned int,4> result(fun_dfn.BookSpecialCase_EE(frac, external_edges, internal_edges, edge_to_cells, ver_to_cells, c,
                                                              l, ext1_tr, ext2_tr, it_l, it_ver));
    ASSERT_TRUE(result.isApprox(expected,1e-10));
}

TEST(BookCase_test,BookCase_VE)
{
    PolygonalMesh frac;
    frac.NumberCell0D = 4; // Numero totale Cell0D
    for (unsigned int i=0; i<4; i++)
        frac.IdCell0D.push_back(i); // identificatori Cell0D --> intero positivo (dimensione 1)
    frac.CoordinatesCell0D.push_back({4,0,0});
    frac.CoordinatesCell0D.push_back({4,3,0});
    frac.CoordinatesCell0D.push_back({0,3,0});
    frac.CoordinatesCell0D.push_back({0,0,0});


    frac.NumberCell1D = 4; // Numero totale Cell1D
    for (unsigned int i=0; i<4; i++)
        frac.IdCell1D.push_back(i); // identificatori Cell1D --> intero positivo (dimensione 1)
    frac.VerticesCell1D.push_back({0,1});
    frac.VerticesCell1D.push_back({1,2});
    frac.VerticesCell1D.push_back({2,3});
    frac.VerticesCell1D.push_back({3,0});

    frac.NumberCell2D = 1; // Numero totale Cell2D
    frac.IdCell2D = {0}; // identificatori Cell2D --> intero positivo (dimensione 1)
    frac.VerticesCell2D.push_back({0,1,2,3});
    frac.EdgesCell2D.push_back({0,1,2,3});

    double tolerance = 10000*numeric_limits<double>::epsilon();

    list<unsigned int> internal_edges={};
    list<unsigned int> external_edges={0,1,2,3};
    unsigned int c=0;
    unsigned int l=3;
    bool v_in_extr1=true;
    Vector3d ext1_tr(2,0,0);
    unsigned int v=3;
    vector<Vector2i> edge_to_cells={{0,-1}, {0,-1}, {0,-1}, {0,-1}};
    vector<list<unsigned int>> ver_to_cells={{0},{0},{0},{0}};
    list<unsigned int>::iterator it_l={frac.EdgesCell2D[c].begin()};
    list<unsigned int>::iterator it_ver={frac.VerticesCell2D[c].begin()};
    Vector<unsigned int, 2> expected(4,4);
    Vector<unsigned int, 2> result(fun_dfn.BookSpecialCase_VE(frac, external_edges, internal_edges, edge_to_cells, ver_to_cells, c,
                                                              l, v_in_extr1, ext1_tr, v, it_l, it_ver));
    ASSERT_TRUE(result.isApprox(expected,1e-10));
}

TEST(BookCase_test,GeneralBookCase)
{
    PolygonalMesh frac;
    frac.NumberCell0D = 4; // Numero totale Cell0D
    for (unsigned int i=0; i<4; i++)
        frac.IdCell0D.push_back(i); // identificatori Cell0D --> intero positivo (dimensione 1)
    frac.CoordinatesCell0D.push_back({4,0,0});
    frac.CoordinatesCell0D.push_back({4,3,0});
    frac.CoordinatesCell0D.push_back({0,3,0});
    frac.CoordinatesCell0D.push_back({0,0,0});


    frac.NumberCell1D = 4; // Numero totale Cell1D
    for (unsigned int i=0; i<4; i++)
        frac.IdCell1D.push_back(i); // identificatori Cell1D --> intero positivo (dimensione 1)
    frac.VerticesCell1D.push_back({0,1});
    frac.VerticesCell1D.push_back({1,2});
    frac.VerticesCell1D.push_back({2,3});
    frac.VerticesCell1D.push_back({3,0});

    frac.NumberCell2D = 1; // Numero totale Cell2D
    frac.IdCell2D = {0}; // identificatori Cell2D --> intero positivo (dimensione 1)
    frac.VerticesCell2D.push_back({0,1,2,3});
    frac.EdgesCell2D.push_back({0,1,2,3});

    double tolerance = 10000*numeric_limits<double>::epsilon();
    list<unsigned int> internal_edges={};
    list<unsigned int> external_edges={0,1,2,3};
    vector<Vector2i> edge_to_cells={{0,-1}, {0,-1}, {0,-1}, {0,-1}};
    vector<list<unsigned int>> ver_to_cells={{0},{0},{0},{0}};
    unsigned int l10=0;
    unsigned int l2=3;
    Vector3d ext_tr(3,0,0);
    unsigned int v0=3;
    unsigned int v2=0;
    bool edge_found0=false;
    bool edge_found2=true;
    Vector<bool,2> expected(true,true);
    Vector<bool,2> result(fun_dfn.GeneralBookCase(frac, external_edges, internal_edges, edge_to_cells, ver_to_cells, l10, l2,
                                                   ext_tr, v0, v2, edge_found0, edge_found2));
    ASSERT_TRUE(result.isApprox(expected,1e-10));
}

TEST(BookCase_test,GeneralBookCaseVV)
{
    PolygonalMesh frac;
    frac.NumberCell0D = 4; // Numero totale Cell0D
    for (unsigned int i=0; i<4; i++)
        frac.IdCell0D.push_back(i); // identificatori Cell0D --> intero positivo (dimensione 1)
    frac.CoordinatesCell0D.push_back({4,0,0});
    frac.CoordinatesCell0D.push_back({4,3,0});
    frac.CoordinatesCell0D.push_back({0,3,0});
    frac.CoordinatesCell0D.push_back({0,0,0});


    frac.NumberCell1D = 4; // Numero totale Cell1D
    for (unsigned int i=0; i<4; i++)
        frac.IdCell1D.push_back(i); // identificatori Cell1D --> intero positivo (dimensione 1)
    frac.VerticesCell1D.push_back({0,1});
    frac.VerticesCell1D.push_back({1,2});
    frac.VerticesCell1D.push_back({2,3});
    frac.VerticesCell1D.push_back({3,0});

    frac.NumberCell2D = 1; // Numero totale Cell2D
    frac.IdCell2D = {0}; // identificatori Cell2D --> intero positivo (dimensione 1)
    frac.VerticesCell2D.push_back({0,1,2,3});
    frac.EdgesCell2D.push_back({0,1,2,3});

    double tolerance = 10000*numeric_limits<double>::epsilon();
    list<unsigned int> internal_edges={};
    list<unsigned int> external_edges={0,1,2,3};
    vector<Vector2i> edge_to_cells={{0,-1}, {0,-1}, {0,-1}, {0,-1}};
    vector<list<unsigned int>> ver_to_cells={{0},{0},{0},{0}};
    unsigned int l10=0;
    unsigned int l2=0;
    Vector3d ext_tr(4,0,0);
    unsigned int v0=3;
    unsigned int v2=0;
    bool edge_found0=false;
    bool edge_found2=false;
    Vector<bool,2> expected(false,false);
    Vector<bool,2> result(fun_dfn.GeneralBookCase(frac, external_edges, internal_edges, edge_to_cells, ver_to_cells, l10, l2,
                                                   ext_tr, v0, v2, edge_found0, edge_found2));
    ASSERT_TRUE(result.isApprox(expected,1e-10));
}

TEST(InternalExternal_test,ExternalEdge)
{
    PolygonalMesh frac;
    frac.NumberCell0D = 4; // Numero totale Cell0D
    for (unsigned int i=0; i<4; i++)
        frac.IdCell0D.push_back(i); // identificatori Cell0D --> intero positivo (dimensione 1)
    frac.CoordinatesCell0D.push_back({5,0,0});
    frac.CoordinatesCell0D.push_back({6,3,0});
    frac.CoordinatesCell0D.push_back({3,3,0});
    frac.CoordinatesCell0D.push_back({2,0,0});


    frac.NumberCell1D = 4; // Numero totale Cell1D
    for (unsigned int i=0; i<4; i++)
        frac.IdCell1D.push_back(i); // identificatori Cell1D --> intero positivo (dimensione 1)
    frac.VerticesCell1D.push_back({0,1});
    frac.VerticesCell1D.push_back({1,2});
    frac.VerticesCell1D.push_back({2,3});
    frac.VerticesCell1D.push_back({3,0});

    frac.NumberCell2D = 1; // Numero totale Cell2D
    frac.IdCell2D = {0}; // identificatori Cell2D --> intero positivo (dimensione 1)
    frac.VerticesCell2D.push_back({0,1,2,3});
    frac.EdgesCell2D.push_back({0,1,2,3});

    double tolerance = 10000*numeric_limits<double>::epsilon();
    list<unsigned int> internal_edges={};
    list<unsigned int> external_edges={0,1,2,3};
    unsigned int id_NEW_E=4;
    unsigned int edge=2;
    fun_dfn.InternalExternalEdge(id_NEW_E,edge,external_edges,internal_edges);
    list<unsigned int> new_edge=external_edges;
    list<unsigned int> expected={0,1,2,3,4};
    ASSERT_EQ(new_edge,expected);
}

TEST(InternalExternal_test,InternalEdge)
{
    PolygonalMesh frac;
    frac.NumberCell0D = 4; // Numero totale Cell0D
    for (unsigned int i=0; i<4; i++)
        frac.IdCell0D.push_back(i); // identificatori Cell0D --> intero positivo (dimensione 1)
    frac.CoordinatesCell0D.push_back({1,-1,0});
    frac.CoordinatesCell0D.push_back({1,1,0});
    frac.CoordinatesCell0D.push_back({-1,1,0});
    frac.CoordinatesCell0D.push_back({-1,-1,0});

    frac.NumberCell1D = 5; // Numero totale Cell1D
    for (unsigned int i=0; i<5; i++)
        frac.IdCell1D.push_back(i); // identificatori Cell1D --> intero positivo (dimensione 1)
    frac.VerticesCell1D.push_back({0,1});
    frac.VerticesCell1D.push_back({1,2});
    frac.VerticesCell1D.push_back({2,3});
    frac.VerticesCell1D.push_back({3,0});
    frac.VerticesCell1D.push_back({3,1});

    frac.NumberCell2D = 2; // Numero totale Cell2D
    frac.IdCell2D.push_back(0);
    frac.IdCell2D.push_back(1);    // identificatori Cell2D --> intero positivo (dimensione 1)
    frac.VerticesCell2D.push_back({0,1,3});
    frac.VerticesCell2D.push_back({1,2,3});
    frac.EdgesCell2D.push_back({0,4,3});
    frac.EdgesCell2D.push_back({1,2,4});

    double tolerance = 10000*numeric_limits<double>::epsilon();
    list<unsigned int> internal_edges={4};
    list<unsigned int> external_edges={0,1,2,3};
    unsigned int id_NEW_E=5;
    unsigned int edge=4;
    fun_dfn.InternalExternalEdge(id_NEW_E,edge,external_edges,internal_edges);
    list<unsigned int> new_edge=internal_edges;
    list<unsigned int> expected={4,5};
    ASSERT_EQ(new_edge,expected);
}

TEST(calculate_fracture_cuts_test,calculate_fracture_cuts_casobase)
{
    Matrix3Xd frac_vertices(3,4);
    frac_vertices << 4, 4, 0, 0,
                    0, 5, 5, 0,
                    0, 0, 0, 0;
    list<unsigned int> p_traces={0};
    list<unsigned int> np_traces={};
    vector<Matrix<double,3,2>> traces_extremes;
    Matrix<double, 3, 2> mat1;
    mat1 << 4, 0,
            3, 3,
            0, 0;
    traces_extremes.push_back(mat1);
    double tol=10000*numeric_limits<double>::epsilon();
    PolygonalMesh res=fun_dfn.calculate_fracture_cuts(frac_vertices,p_traces,np_traces,traces_extremes, tol);
    unsigned int NumberCell0D=6;
    vector<unsigned int> IdCell0D = {0,1,2,3,4,5};
    vector<Vector3d> CoordinatesCell0D = {{4,0,0},{4,5,0},{0,5,0},{0,0,0},{4,3,0},{0,3,0}};

    unsigned int NumberCell1D = 7;
    vector<unsigned int> IdCell1D = {0,1,2,3,4,5,6};
    vector<Vector2i> VerticesCell1D = {{0,4},{1,2},{2,5},{3,0},{4,1},{5,4},{5,3}};

    unsigned int NumberCell2D = 2;
    vector<unsigned int> IdCell2D = {0,1};
    vector<list<unsigned int>> VerticesCell2D = {{0,4,5,3}, {4,1,2,5}};
    vector<list<unsigned int>> EdgesCell2D = {{0,5,6,3}, {4,1,2,5}};
    ASSERT_EQ(NumberCell0D,res.NumberCell0D);
    ASSERT_EQ(IdCell0D,res.IdCell0D);
    ASSERT_EQ(CoordinatesCell0D,res.CoordinatesCell0D);
    ASSERT_EQ(NumberCell1D,res.NumberCell1D);
    ASSERT_EQ(IdCell1D,res.IdCell1D);
    ASSERT_EQ(VerticesCell1D,res.VerticesCell1D);
    ASSERT_EQ(NumberCell2D,res.NumberCell2D);
    ASSERT_EQ(IdCell2D,res.IdCell2D);
    ASSERT_EQ(VerticesCell2D,res.VerticesCell2D);
    ASSERT_EQ(EdgesCell2D,res.EdgesCell2D);
}


#endif

//.inp

