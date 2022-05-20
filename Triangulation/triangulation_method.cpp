/**
 * Copyright (C) 2015 by Liangliang Nan (liangliang.nan@gmail.com)
 * https://3d.bk.tudelft.nl/liangliang/
 *
 * This file is part of Easy3D. If it is useful in your research/work,
 * I would be grateful if you show your appreciation by citing it:
 * ------------------------------------------------------------------
 *      Liangliang Nan.
 *      Easy3D: a lightweight, easy-to-use, and efficient C++
 *      library for processing and rendering 3D data. 2018.
 * ------------------------------------------------------------------
 * Easy3D is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License Version 3
 * as published by the Free Software Foundation.
 *
 * Easy3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "triangulation.h"
#include "matrix_algo.h"
#include <easy3d/optimizer/optimizer_lm.h>


using namespace easy3d;


/**
 * TODO: Finish this function for reconstructing 3D geometry from corresponding image points.
 * @return True on success, otherwise false. On success, the reconstructed 3D points must be written to 'points_3d'
 *      and the recovered relative pose must be written to R and t.
 */
bool Triangulation::triangulation(
        double fx, double fy,     /// input: the focal lengths (same for both cameras)
        double cx, double cy,     /// input: the principal point (same for both cameras)
        const std::vector<Vector2D> &points_0,  /// input: 2D image points in the 1st image.
        const std::vector<Vector2D> &points_1,  /// input: 2D image points in the 2nd image.
        std::vector<Vector3D> &points_3d,       /// output: reconstructed 3D points
        Matrix33 &R,   /// output: 3 by 3 matrix, which is the recovered rotation of the 2nd camera
        Vector3D &t    /// output: 3D vector, which is the recovered translation of the 2nd camera
) const
{
    /// NOTE: there might be multiple workflows for reconstructing 3D geometry from corresponding image points.
    ///       This assignment uses the commonly used one explained in our lecture.
    ///       It is advised to define a function for the sub-tasks. This way you have a clean and well-structured
    ///       implementation, which also makes testing and debugging easier. You can put your other functions above
    ///       triangulation(), or put them in one or multiple separate files.

    std::cout << "\nTODO: I am going to implement the triangulation() function in the following file:" << std::endl
              << "\t    - triangulation_method.cpp\n\n";

    std::cout << "[Liangliang]:\n"
                 "\tFeel free to use any provided data structures and functions. For your convenience, the\n"
                 "\tfollowing three files implement basic linear algebra data structures and operations:\n"
                 "\t    - Triangulation/matrix.h  Matrices of arbitrary dimensions and related functions.\n"
                 "\t    - Triangulation/vector.h  Vectors of arbitrary dimensions and related functions.\n"
                 "\t    - Triangulation/matrix_algo.h  Determinant, inverse, SVD, linear least-squares...\n"
                 "\tPlease refer to the above files for a complete list of useful functions and their usage.\n\n"
                 "\tIf you choose to implement the non-linear method for triangulation (optional task). Please\n"
                 "\trefer to 'Tutorial_NonlinearLeastSquares/main.cpp' for an example and some explanations.\n\n"
                 "\tIn your final submission, please\n"
                 "\t    - delete ALL unrelated test or debug code and avoid unnecessary output.\n"
                 "\t    - include all the source code (and please do NOT modify the structure of the directories).\n"
                 "\t    - do NOT include the 'build' directory (which contains the intermediate files in a build step).\n"
                 "\t    - make sure your code compiles and can reproduce your results without ANY modification.\n\n" << std::flush;

    /// Below are a few examples showing some useful data structures and APIs.

    /// define a 2D vector/point
    Vector2D b(1.1, 2.2);

    /// define a 3D vector/point
    Vector3D a(1.1, 2.2, 3.3);

    /// get the Cartesian coordinates of a (a is treated as Homogeneous coordinates)
    Vector2D p = a.cartesian();

    /// get the Homogeneous coordinates of p
    Vector3D q = p.homogeneous();

    /// define a 3 by 3 matrix (and all elements initialized to 0.0)
    Matrix33 A;

    /// define and initialize a 3 by 3 matrix
    Matrix33 T(1.1, 2.2, 3.3,
               0, 2.2, 3.3,
               0, 0, 1);

    /// define and initialize a 3 by 4 matrix
    Matrix34 M(1.1, 2.2, 3.3, 0,
               0, 2.2, 3.3, 1,
               0, 0, 1, 1);

    /// set first row by a vector
    M.set_row(0, Vector4D(1.1, 2.2, 3.3, 4.4));

    /// set second column by a vector
    M.set_column(1, Vector3D(5.5, 5.5, 5.5));

    /// define a 15 by 9 matrix (and all elements initialized to 0.0)
    Matrix W(15, 9, 0.0);
    /// set the first row by a 9-dimensional vector
    W.set_row(0, {0, 1, 2, 3, 4, 5, 6, 7, 8}); // {....} is equivalent to a std::vector<double>

    /// get the number of rows.
    int num_rows = W.rows();

    /// get the number of columns.
    int num_cols = W.cols();

    /// get the the element at row 1 and column 2
    double value = W(1, 2);

    /// get the last column of a matrix
    Vector last_column = W.get_column(W.cols() - 1);

    /// define a 3 by 3 identity matrix
    Matrix33 I = Matrix::identity(3, 3, 1.0);

    /// matrix-vector product
    Vector3D v = M * Vector4D(1, 2, 3, 4); // M is 3 by 4

    ///For more functions of Matrix and Vector, please refer to 'matrix.h' and 'vector.h'

    // TODO: delete all above example code in your final submission

    //--------------------------------------------------------------------------------------------------------------
    // implementation starts ...

    // TODO: check if the input is valid (always good because you never known how others will call your function).
    if (size(points_0)==size(points_1)){
        std::cout << "Input is valid, triangulation will be applied" << std::endl;
    } else {
        std::cout << "Input is invalid" << std::endl;
    }

    // TODO: Estimate relative pose of two views. This can be subdivided into
    //      - estimate the fundamental matrix F;
    //      - compute the essential matrix E;
    //      - recover rotation R and t.

    //------------estimate the fundamental matrix F ---------------
    //---normalization---
    //calculation of the new centroids
    float centr_x0=0;
    float centr_y0=0;
    float centr_x1=0;
    float centr_y1=0;

    for (int i=0; i<size(points_0); i++){
        centr_x0 = centr_x0 + points_0[i][0];
        centr_y0 = centr_y0 + points_0[i][1];
        centr_x1 = centr_x1 + points_1[i][0];
        centr_y1 = centr_y1 + points_1[i][1];
    }

    centr_x0 = centr_x0/size(points_0);
    centr_y0 = centr_y0/size(points_0);
    centr_x1 = centr_x1/size(points_1);
    centr_y1 = centr_y1/size(points_1);


    //translation and scaling
    double mean_dist_x0 = 0;
    double mean_dist_y0 = 0;
    double mean_dist_x1 = 0;
    double mean_dist_y1 = 0;

    Matrix33 M_trans_0; //initializing translation matrix for image 0
    M_trans_0.set_column(0, {1, 0, 0});
    M_trans_0.set_column(1, {0, 1, 0});
    M_trans_0.set_column(2, {(-centr_x0), (-centr_y0), 1});

    Matrix33 M_trans_1; //initializing translation matrix for image 1
    M_trans_1.set_column(0, {1, 0, 0});
    M_trans_1.set_column(1, {0, 1, 0});
    M_trans_1.set_column(2, {(-centr_x1), (-centr_y1), 1});

    //initializing vectors to store the points after translation
    std::vector<Vector2D> vect_trans_pt_0(points_0);
    std::vector<Vector2D> vect_trans_pt_1(points_1);

    for (int i=0; i<size(points_0); i++){
        //translation for image 0
        Matrix pt_0(3,1,0.0);
        pt_0.set_column(0,{points_0[i][0],points_0[i][1],1});
        Matrix pt_0_trans = M_trans_0 * pt_0;

        mean_dist_x0 = mean_dist_x0 + sqrt (pow(pt_0_trans(0,0),2)); //calculating mean distance
        mean_dist_y0 = mean_dist_y0 + sqrt (pow(pt_0_trans(1,0),2));

        vect_trans_pt_0[i][0] = (pt_0_trans(0,0)); //filing the translated points vector
        vect_trans_pt_0[i][1] = (pt_0_trans(1,0));


        //translation for image 0
        Matrix pt_1(3,1,0.0);
        pt_1.set_column(0,{points_1[i][0],points_0[i][1],1});
        Matrix pt_1_trans = M_trans_0 * pt_1;

        mean_dist_x1 = mean_dist_x1 + sqrt (pow(pt_1_trans(0,0),2));
        mean_dist_y1 = mean_dist_y1 + sqrt (pow(pt_1_trans(1,0),2));

        vect_trans_pt_1[i][0] = (pt_1_trans(0,0));
        vect_trans_pt_1[i][1] = (pt_1_trans(1,0));

    }

    //scaling
    //calculating scaling factors
    mean_dist_x0 = mean_dist_x0 / size(points_0);
    mean_dist_y0 = mean_dist_y0 / size(points_0);
    mean_dist_x1 = mean_dist_x1 / size(points_1);
    mean_dist_y1 = mean_dist_y1 / size(points_1);

    double sc_fc_x0 = sqrt(2) / mean_dist_x0;
    double sc_fc_y0 = sqrt(2) / mean_dist_y0;
    double sc_fc_x1 = sqrt(2) / mean_dist_x1;
    double sc_fc_y1 = sqrt(2) / mean_dist_y1;

    Matrix m_sf_0(2,2,0.0); //scaling matrix for image 0
    //m_sf_0.set_column(0, {sc_fc_x0, 0});
    //m_sf_0.set_column(1, {0, sc_fc_y0});
    m_sf_0(0,0) = sc_fc_x0;
    m_sf_0(0,1) = 0;
    m_sf_0(1,0) = 0;
    m_sf_0(1,1) = sc_fc_y0;

    Matrix m_sf_1(2,2,0.0); //scaling matrix for image 1
    //m_sf_1.set_column(0, {sc_fc_x1, 0});
    //m_sf_1.set_column(1, {0, sc_fc_y1});
    m_sf_1(0,0) = sc_fc_x1;
    m_sf_1(0,1) = 0;
    m_sf_1(1,0) = 0;
    m_sf_1(1,1) = sc_fc_y1;

    std::vector<Vector2D> vect_scal_pt_0(points_0); //initializing vector to store the scaled points for image 0
    std::vector<Vector2D> vect_scal_pt_1(points_1);

    for (int i=0; i<size(points_0); i++){
        //applying the scaling
        Matrix sc_fc_pt0(2,1,0.0); //initializing a matrix with the translated points for image 0
        //sc_fc_pt0.set_column(0,{vect_trans_pt_0[i][0],vect_trans_pt_0[i][1]});
        sc_fc_pt0(0,0) = vect_trans_pt_0[i][0];
        sc_fc_pt0(1,0) = vect_trans_pt_0[i][1];

        Matrix pt_0_scaled =  m_sf_0 * sc_fc_pt0;
        vect_scal_pt_0[i][0] = (pt_0_scaled(0,0));
        vect_scal_pt_0[i][1] = (pt_0_scaled(1,0));

        Matrix sc_fc_pt1(2,1,0.0);
        //sc_fc_pt1.set_column(0,{vect_trans_pt_1[i][0],vect_trans_pt_1[i][1]});
        sc_fc_pt1(0,0) = vect_trans_pt_1[i][0];
        sc_fc_pt1(1,0) = vect_trans_pt_1[i][1];
        Matrix pt_1_scaled =  m_sf_1 * sc_fc_pt1;
        vect_scal_pt_1[i][0] = (pt_1_scaled(0,0));
        vect_scal_pt_1[i][1] = (pt_1_scaled(1,0));
        }

    //crating matrix W
    Matrix m_W(size(points_0), 8);
    for (int i=0; i<size(points_0); i++){
        m_W.set_row(i,{vect_scal_pt_0[i][0]*vect_scal_pt_1[i][0], vect_scal_pt_0[i][1]* vect_scal_pt_1[i][0], vect_scal_pt_1[i][0], vect_scal_pt_0[i][0]*vect_scal_pt_1[i][1], vect_scal_pt_0[i][1]*vect_scal_pt_1[i][1], vect_scal_pt_1[i][1], vect_scal_pt_0[i][0], 1});
    }

    std::cout << m_W << "m_W" << std::endl;



















    // TODO: Reconstruct 3D points. The main task is
    //      - triangulate a pair of image points (i.e., compute the 3D coordinates for each corresponding point pair)

    // TODO: Don't forget to
    //          - write your recovered 3D points into 'points_3d' (so the viewer can visualize the 3D points for you);
    //          - write the recovered relative pose into R and t (the view will be updated as seen from the 2nd camera,
    //            which can help you check if R and t are correct).
    //       You must return either 'true' or 'false' to indicate whether the triangulation was successful (so the
    //       viewer will be notified to visualize the 3D points and update the view).
    //       There are a few cases you should return 'false' instead, for example:
    //          - function not implemented yet;
    //          - input not valid (e.g., not enough points, point numbers don't match);
    //          - encountered failure in any step.
    return points_3d.size() > 0;
}