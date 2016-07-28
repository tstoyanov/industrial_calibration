#ifndef CURVE_FIT_HPP
#define CURVE_FIT_HPP

#include "ceres/ceres.h"

  class  SecondOrderStep
  {
  public:
    SecondOrderStep(double time, double value) :
      time_(time), value_(value)
    {
    }
    template<typename T>
    bool operator()(	    const T* const k,  /** gain k*/
		    const T* const terms,  /** zeta, omega_n*/
		    T* residual) const
    {
      T K      = k[0];
      T zeta = terms[0];
      T wn   = terms[1];
      T phi = atan(zeta);
      T sterm = sqrt(T(1)-zeta*zeta);
      T t(time_);
      T v(value_);
      residual[0] =  v -  K*(T(1)-exp(-zeta*wn*t)*sin(wn*sterm*t + phi)/sin(phi));
      return true;
    } /** end of operator() */

    /** Factory to hide the construction of the CostFunction object from */
    /** the client code. */
    static ceres::CostFunction* Create(const double time, const double value)
    {
      return (new ceres::AutoDiffCostFunction<SecondOrderStep, 1, 1, 2>(new SecondOrderStep(time, value)));
    }
    double value_; /** value observed */
    double time_; /** time it was observed */
  }; 

#endif
