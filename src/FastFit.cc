/**
 * Thomas Keck 2016
 */

#include <FastFit.h>


Helix::Helix(double magnetic_field, const Eigen::Matrix<double, 3, 1> &x, const Eigen::Matrix<double, 3, 1> &p) {
  
 /**
  * Kappa as defined in Data Analysis Techniques for High Energy Physics p. 258,
  * but divided by 100 to convert from m to cm.
  *
  * It is the same number that KFit uses as well (see KFitConst.h)
  *
  * It depends on the units we are using:
  * Magnetic field in T
  * Distances in cm
  * Momenta in GeV/c
  * Charges in C/e
  */
  const double kappa = -0.00299792458; 
  const double alpha = magnetic_field * kappa * charge;

  // Working with the inverse transverse momentum should be safe,
  // because tracks with a transverse momentum of 0 shouldn't be in
  // the detector acceptance!
  const double ptinv = 1.0 / std::sqrt(p(X)*p(X) + p(Y)*p(Y));
  if (not std::isfinite(ptinv)) {
    std::cerr << "Invalid inverse transverse momentum. This should happen." << std::endl;
    return false;
  }

  m_dP = (-x(X) * p(X) - x(Y) * p(Y)) * ptinv;
  m_dO = (-x(Y) * p(X) + x(X) * p(Y)) * ptinv;
  m_dR = std::sqrt(x(X)*x(X) + x(Y)*x(Y));
  m_area = 2 * m_dO + alpha  * m_ptinv * m_dR * m_dR;
  m_denom 1 + std::sqrt(1 + alpha * ptinv * m_area);

  m_chi = 0;
  m_arcLength = m_dP;
  
  if (alpha != 0) {
      m_chi = - std::atan2(alpha * ptinv * m_dP, 1 + alpha * ptinv * m_dO);
      m_arcLength = - m_chi / (alpha * ptinv);
  } 

  m_parametrisation(PtInverse) = ptinv;
  m_parametrisation(TanLambda) = p(Z) * ptinv;
  m_parametrisation(PhiAngle0) = std::atan2(p(Y), p(X)) + m_chi;
  m_parametrisation(Vertical0) = x(Z) + m_parametrisation(TanLambda) * m_arcLength;
  m_parametrisation(Distance0) = m_area / m_denom;

}


FastFit::FastFit(unsigned int numberOfDaughters) {

  m_numberOfDaughters = numberOfDaughters;

  m_charges.resize(m_numberOfDaughters);
  m_momenta.resize(m_numberOfDaughters);
  m_smoothed_momenta.resize(m_numberOfDaughters);
  m_positions.resize(m_numberOfDaughters);
  m_variances.resize(m_numberOfDaughters);

  m_A.resize(m_numberOfDaughters);
  m_B.resize(m_numberOfDaughters);
  m_c.resize(m_numberOfDaughters);
  m_G.resize(m_numberOfDaughters);
  m_GB.resize(m_numberOfDaughters);
  m_S.resize(m_numberOfDaughters);
  m_measurement.resize(m_numberOfDaughters);
  
  m_vertex = Eigen::Matrix<double, 3, 1>::Zero();
}

bool FastFit::fit(unsigned int maximumNumberOfFitIterations, double magnetic_field)
{

  Eigen::Matrix<double, 5, 6> J;
  Eigen::Matrix<double, 5, 5> G_inv;
  Eigen::Matrix<double, 3, 3> S_inv;

  unsigned int iteration;
  for (iteration = 0; iteration < maximumNumberOfFitIterations; ++iteration) {

    Eigen::Matrix<double, 3, 1> current_vertex = m_vertex;

    // Initial guess for vertex is just the IP
    // with a huge uncertainty of 1 m
    m_vertex = Eigen::Matrix<double, 3, 1>::Zero();
    m_C << 100.0,   0.0,   0.0,
        0.0, 100.0,   0.0,
        0.0,   0.0, 100.0;
    m_C_inv << 0.01, 0.0, 0.0,
            0.0, 0.01, 0.0,
            0.0, 0.0, 0.01;

    for (unsigned int i = 0; i < m_numberOfDaughters; ++i) {

      auto charge = m_charges[i];
      auto& p = m_momenta[i];
      auto& p_s = m_smoothed_momenta[i];
      auto& x = m_positions[i];
      auto& V = m_variances[i];

      auto& A = m_A[i];
      auto& B = m_B[i];
      auto& c = m_c[i];

      auto& G = m_G[i];
      auto& GB = m_GB[i];

      auto& S = m_S[i];
      auto& measurement = m_measurement[i];

      Helix helix(magnetic_field, x, p);

      measurement << helix.get_parametrisation();

      A << 1, 0, - p_s(X) / p_s(Z),
      0, 1, - p_s(Y) / p_s(Z),
      0, 0, alpha* p_s(Y) / p_s(Z),
      0, 0, - alpha* p_s(X) / p_s(Z),
      0, 0, 0;

      B << 0, 0, 0,
      0, 0, 0,
      1, 0, 0,
      0, 1, 0,
      0, 0, 1;

      c << reference_z* p_s(X) / p_s(Z),
      reference_z* p_s(Y) / p_s(Z),
      - alpha* reference_z* p_s(Y) / p_s(Z),
      alpha* reference_z* p_s(X) / p_s(Z),
      0;

      J.block<5, 3>(0, 0) = A;
      J.block<5, 3>(0, 3) = B;
      G_inv = J * V * J.transpose();

      if (invertAndCheck(G_inv, G) == false) {
        std::cout << "Affected Matrix was G" << std::endl;
        std::cout << G_inv << std::endl << std::endl;
        std::cout << J << std::endl << std::endl;
        std::cout << V << std::endl << std::endl;
        return false;
      }

      S_inv = B.transpose() * G * B;

      if (invertAndCheck(S_inv, S) == false) {
        std::cout << "Affected Matrix was S" << std::endl;
        return false;
      }

      GB = G - G * B * S * B.transpose() * G; // G is symmetric, doesn't need to be transposed

      m_vertex += A.transpose() * GB * (measurement - c);
      m_C_inv += A.transpose() * GB * A;
    }

    if (invertAndCheck(m_C_inv, m_C) == false) {
      std::cout << "Affected Matrix was m_C" << std::endl;
      return false;
    }

    m_vertex = m_C * m_vertex;

    for (unsigned int i = 0; i < m_numberOfDaughters; ++i) {

      auto& p_s = m_smoothed_momenta[i];

      auto& A = m_A[i];
      auto& B = m_B[i];
      auto& c = m_c[i];

      auto& G = m_G[i];
      auto& S = m_S[i];

      auto& measurement = m_measurement[i];

      // New smoothed momentum at the final vertex position
      p_s = S * B.transpose() * G * (measurement - c - A * m_vertex);

    }
  }

  m_chi2 = 0.0;
  m_ndf = 0;

  Eigen::Matrix<double, 5, 1> r;
  Eigen::Matrix<double, 3, 3> E;

  for (unsigned int i = 0; i < m_numberOfDaughters; ++i) {

    auto& p_s = m_smoothed_momenta[i];
    auto& V = m_variances[i];

    auto& A = m_A[i];
    auto& B = m_B[i];
    auto& c = m_c[i];

    auto& G = m_G[i];
    auto& S = m_S[i];

    auto& measurement = m_measurement[i];

    E = - m_C * A.transpose() * G * B * S;
    // TODO m_C has to be retransformed to the original X variance position!
    V.block<3, 3>(0, 0) = m_C;
    V.block<3, 3>(3, 3) = S + E.transpose() * m_C_inv * E;
    V.block<3, 3>(3, 0) = E;
    V.block<3, 3>(0, 3) = E;

    r = measurement - c - A * m_vertex - B * p_s;
    m_chi2 += (r.transpose() * G * r)(0, 0);

  }

  m_ndf = 2 * m_numberOfDaughters - 3;

  return true;
}


double FastFit::GetDaughterVariance(unsigned int i, unsigned int component_i, unsigned int component_j) const {
    if(i >= m_numberOfDaughters) {
      throw std::runtime_error("Invalid index passed to GetDaughterMomentum");
    }
    if(component_i >= 6) {
      throw std::runtime_error("Invalid index passed to GetDaughterMomentum");
    }
    if(component_j >= 6) {
      throw std::runtime_error("Invalid index passed to GetDaughterMomentum");
    }
    const auto &V = m_variances[i];
    return V(component_i, component_j);
}
    

double FastFit::GetDaughterMomentum(unsigned int i, unsigned int component) const {
    if(i >= m_numberOfDaughters) {
      throw std::runtime_error("Invalid index passed to GetDaughterMomentum");
    }
    if(component >= 3) {
      throw std::runtime_error("Invalid index passed to GetDaughterMomentum");
    }
    const auto &p = m_smoothed_momenta[i];
    return p(component);
}


double FastFit::GetVertex(unsigned int component) const {
    if(component >= 3) {
      throw std::runtime_error("Invalid index passed to GetVertex");
    }
    return m_vertex(component);
}
