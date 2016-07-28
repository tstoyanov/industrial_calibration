#include <stdio.h>
#include "ceres/ceres.h"
#include "ceres/types.h"
#include "ceres/rotation.h"
#include "intrinsic_cal/curve_fit.hpp"

using std::string;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;

main()
{
  int num_pts =100;
  Problem problem; 
  int num_pnts = 10;
  int num_data_sets = 1;
  double P1[2]; // zeta and omega_n
  double P2[num_data_sets];// K, one for each data set
  double t[num_pnts];
  double v[num_pnts];
  double dt = .001;
  double wn = 10.1;
  double zeta = 0.707;
  double K = 1.3;
  double phi = atan(zeta);
  for(int i=0;i<num_pts;i++){
    t[i] = i*dt;
    v[i] = K*(1-exp(-zeta*wn*t[i])*sin(wn*t[i] + phi)/sin(phi));
  }
  double *params;
  CostFunction* cost_function[num_pnts]; // not sure I need a new cost function each time, but this just uses memory
  for(int i=0;i<num_pts;i++){
    cost_function[i] = SecondOrderStep::Create(t[i],v[i]);
    problem.AddResidualBlock(cost_function[i], NULL, params);
  }

  Solver::Options options;
  Solver::Summary summary;
  options.linear_solver_type = ceres::DENSE_SCHUR;
  options.minimizer_progress_to_stdout = true;
  options.max_num_iterations = 2000;
  ceres::Solve(options, &problem, &summary);
  if(summary.termination_type != ceres::NO_CONVERGENCE){
    printf("No convergence\n");
  }
  else{
    printf("convergence\n");
  }

}
