// What peira's quasi-experiments share with tyche: least squares and their robust covariances (tyche::regression, its sandwich) and
// panels with absorbed effects (tyche::panel), named as the estimators here use them.

#pragma once

#include <tyche/panel/panel.hpp>
#include <tyche/regression/linear.hpp>
#include <tyche/regression/sandwich.hpp>

#include <calx/linalg/rectangular.hpp>

#include <stdexcept>
#include <string>

namespace peira::quasi {

namespace tr = tyche::regression;
using calx::linalg::Matrix;
using calx::linalg::Rectangular;
using tyche::panel::PanelCovariance;
using tyche::panel::PanelEffects;
using tyche::panel::PanelOptions;
using tyche::panel::within;
using tyche::panel::WithinFit;

namespace detail {
inline void require(bool ok, const char* what) {
  if (!ok) throw std::invalid_argument(std::string("peira::quasi: ") + what);
}
using tyche::regression::join_columns;
using tyche::regression::sandwich;
using tyche::regression::sandwich_meat;
using tyche::regression::solve_spd;
}  // namespace detail

}  // namespace peira::quasi
