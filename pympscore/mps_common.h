#pragma once

// namespace
#define MPS_NAMESPACE_BEGIN namespace mps {
#define MPS_NAMESPACE_END }

//
#include <string>


// linear algebra
#include <Eigen/Dense>

////////////////////////////////////////////////////////////////////////////////
MPS_NAMESPACE_BEGIN;

// vectors and matrices
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::VectorXi;
using Eigen::MatrixXi;

using Eigen::Vector3d;
using Eigen::Matrix3d;
using Eigen::Vector3i;
using Eigen::Matrix3i;

// use Eigen::Ref to pass vectors between Cpp and Python
template<typename EigenType>
using EigenRef = Eigen::Ref<EigenType>;

template<typename EigenType>
using EigenConstRef = Eigen::Ref<const EigenType>;

//----------------------------------------

//extern int SpaceDim;

int GetSpaceDim();
void SetSpaceDim(int ndim);

VectorXd MakeSpaceDimVec();
MatrixXd MakeSpaceDimMat();


MPS_NAMESPACE_END;
////////////////////////////////////////////////////////////////////////////////

