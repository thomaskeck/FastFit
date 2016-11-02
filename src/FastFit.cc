/**
 * Thomas Keck 2016
 */

#include <FastFit.h>
#include <cmath>

Helix::Helix(const double alpha, const Eigen::Matrix<double, 3, 1> &x, const Eigen::Matrix<double, 3, 1> &p) : m_alpha(alpha) {
  
  // Transverse momentum squared
  const double T = p(X)*p(X) + p(Y)*p(Y);
  
  // Working with the inverse transverse momentum should be safe,
  // because tracks with a transverse momentum of 0 shouldn't be in
  // the detector acceptance!
  if (T <= 0) {
    std::cerr << "Invalid inverse transverse momentum. This should happen." << std::endl;
    return;
  }

  // Transverse distance to origin squared
  const double V = x(X)*x(X) + x(Y)*x(Y);
  // Orthogonal projection of x and p
  const double P = p(X) * x(Y) - p(Y) * x(X);
  // Discriminant
  const double D = std::sqrt(T + alpha * (2 * P + alpha * V));

  // 0. Inverse transverse momentum
  m_parametrisation(0) = 1.0 / std::sqrt(T);

  // 1. Slope angle between p_z and p_t
  m_parametrisation(1) = p(Z) / std::sqrt(T);

  // 2. Direction angle between p_t and x-axis
  m_parametrisation(2) = std::atan2(p(Y) - x(X) * alpha, p(X) + x(Y) * alpha);

  // 3. 2d Distance of the perigee to the origin in the x-y plane
  m_parametrisation(3) = (2 * P + alpha * V) / (T + D);

  // 4. vertical distance of perigee to the x-y plane
  double chi = 0;
  if (alpha != 0) {
      chi = std::atan2(alpha *(p(X) * x(X) + p(Y) * x(Y)), T + alpha * P) / alpha;
  } else {
      chi = (p(X) * x(X) + p(Y) * x(Y)) / T;
  }
  m_parametrisation(4) = x(Z) - p(Z) * chi;

  // Some temporary constants to keep the complicated Jacobian clean
  const double a1 = (2 * D * (D + T) - alpha * (2 * P + alpha * V)) / (D * (D + T) * (D + T));
  const double a2 = alpha / ((alpha * x(X) - p(Y)) * (alpha * x(X) - p(Y)) + (alpha * x(Y) + p(X)) * (alpha * x(Y) + p(X)));
  const double a3 = p(Z) / (2 * P * alpha + T + V * alpha * alpha);
  m_A <<
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
       -(alpha * x(Y) + p(X)) * a2, (alpha * x(X) - p(Y)) * a2, 0.0,
        (alpha * x(X) - p(Y)) * a1, (alpha * x(Y) + p(X)) * a1, 0.0,
       -(alpha * x(Y) + p(X)) * a3, (alpha * x(X) - p(Y)) * a3, 1.0;

  const double b1 = 1 / (T * std::sqrt(T));
  const double b2 = 1.0 / ((alpha * x(X) - p(Y)) * (alpha * x(X) - p(Y)) + (alpha * x(Y) + p(X)) * (alpha * x(Y) + p(X)));
  const double b3 = p(Z) / (alpha * alpha * ( p(X) * x(X) + p(Y) * x(Y) ) * ( p(X) * x(X) + p(Y) * x(Y) ) + (P * alpha + T) * (P * alpha + T));
  const double b4 = 1.0 / (D * (D + T) * (D + T));
  m_B <<
       -p(X) * b1, -p(Y) * b1, 0.0,
       -p(X) * p(Z) * b1, -p(Y) * p(Z) * b1, 1.0 / std::sqrt(T),
        (alpha * x(X) - p(Y)) * b2, (alpha * x(Y) + p(X)) * b2, 0.0,
        (2 * D * x(Y) * (D + T) - (2 * P + V * alpha) * (2 * D * p(X) + alpha * x(Y) + p(X))) * b4,
       -(2 * D * x(X) * (D + T) + (2 * P + V * alpha) * (2 * D * p(Y) - alpha * x(X) + p(Y))) * b4,
       0.0,
       -(x(X)*(P * alpha + T) - (alpha * x(Y) + 2 * p(X)) * (p(X) * x(X) + p(Y) * x(Y))) * b3,
       -(x(Y)*(P * alpha + T) + (alpha * x(X) - 2 * p(Y)) * (p(X) * x(X) + p(Y) * x(Y))) * b3,
       -chi;

  m_c0 = m_parametrisation - m_A * x - m_B * p; 

}

/**
 * Centers the angle of the 5d - parametrisation around the given center.
 * This means we use the 2*pi periodicy of p(2) to ensure that max(|p(2) - center|) <= pi
 */
Eigen::Matrix<double, 5, 1> centerAngle(Eigen::Matrix<double, 5, 1> p, double center) {

    p(2) = std::fmod(p(2) - center + M_PI, (2 * M_PI)) + center - M_PI;
    return p;

}

FastFit::FastFit(unsigned int numberOfDaughters, double magnetic_field) {

  m_numberOfDaughters = numberOfDaughters;
  m_magnetic_field = magnetic_field;

  m_smoothed_momenta.resize(m_numberOfDaughters);
  m_original_helix.resize(m_numberOfDaughters);
  m_current_helix.resize(m_numberOfDaughters);
  m_variances.resize(m_numberOfDaughters);

  m_G.resize(m_numberOfDaughters);
  m_GB.resize(m_numberOfDaughters);
  m_S.resize(m_numberOfDaughters);

  // Initial guess for vertex is just the IP
  // with a very large uncertainty (this makes the fit more stable)
  m_use_ip_constraint = false;
  m_ip_vertex = Eigen::Matrix<double, 3, 1>::Zero();
  m_ip_variance << 1000.0,   0.0,   0.0,
                       0.0, 1000.0,   0.0,
                       0.0,   0.0, 1000.0;
  m_ip_variance_inv << 0.001,   0.0,   0.0,
                       0.0, 0.001,   0.0,
                       0.0,   0.0, 0.001;
  
  m_vertex = Eigen::Matrix<double, 3, 1>::Zero();
}

bool FastFit::fit(unsigned int maximumNumberOfFitIterations)
{

  Eigen::Matrix<double, 5, 6> J;
  Eigen::Matrix<double, 5, 5> G_inv;
  Eigen::Matrix<double, 3, 3> S_inv;
 
  unsigned int iteration;
  for (iteration = 0; iteration < maximumNumberOfFitIterations; ++iteration) {

    Eigen::Matrix<double, 3, 1> current_vertex = m_vertex;

    m_vertex = m_ip_variance_inv * m_ip_vertex;
    m_C_inv = m_ip_variance_inv;

    for (unsigned int i = 0; i < m_numberOfDaughters; ++i) {

      const auto& p_s = m_smoothed_momenta[i];
      const auto& V = m_variances[i];
      auto& G = m_G[i];
      auto& GB = m_GB[i];
      auto& S = m_S[i];
     
      if(iteration == 0) {
        m_current_helix[i] = m_original_helix[i];
      } else {
        m_current_helix[i] = Helix(m_original_helix[i].GetAlpha(), current_vertex, p_s);
      }
      const auto &measurement = m_original_helix[i].GetParametrisation();
      const auto &A = m_current_helix[i].GetVertexJacobian();
      const auto &B = m_current_helix[i].GetMomentumJacobian();
      const auto &c = m_current_helix[i].GetConstantOffset();

      J.block<5, 3>(0, 0) = A;
      J.block<5, 3>(0, 3) = B;
      G_inv = J * V * J.transpose();

      if (invertAndCheck(G_inv, G) == false) {
        std::cout << "Affected Matrix was G" << std::endl;
        return false;
      }

      S_inv = B.transpose() * G * B;

      if (invertAndCheck(S_inv, S) == false) {
        std::cout << "Affected Matrix was S" << std::endl;
        return false;
      }

      GB = G - G * B * S * B.transpose() * G; // G is symmetric, doesn't need to be transposed

      m_vertex += A.transpose() * GB * centerAngle(measurement - c, measurement(2));
      m_C_inv += A.transpose() * GB * A;
    }

    if (invertAndCheck(m_C_inv, m_C) == false) {
      std::cout << "Affected Matrix was m_C" << std::endl;
      return false;
    }

    m_vertex = m_C * m_vertex;

    for (unsigned int i = 0; i < m_numberOfDaughters; ++i) {

      auto& p_s = m_smoothed_momenta[i];

      const auto &measurement = m_original_helix[i].GetParametrisation();
      const auto &A = m_current_helix[i].GetVertexJacobian();
      const auto &B = m_current_helix[i].GetMomentumJacobian();
      const auto &c = m_current_helix[i].GetConstantOffset();

      const auto& G = m_G[i];
      const auto& S = m_S[i];

      // New smoothed momentum at the final vertex position
      p_s = S * B.transpose() * G * centerAngle(measurement - c - A * m_vertex, measurement(2));
    }
  }

  m_chi2 = 0.0;
  m_ndf = 0;

  Eigen::Matrix<double, 5, 1> r;
  Eigen::Matrix<double, 3, 3> E;

  for (unsigned int i = 0; i < m_numberOfDaughters; ++i) {

    const auto& p_s = m_smoothed_momenta[i];
    auto& V = m_variances[i];

    const auto& G = m_G[i];
    const auto& S = m_S[i];

    const auto &measurement = m_original_helix[i].GetParametrisation();
    const auto &A = m_current_helix[i].GetVertexJacobian();
    const auto &B = m_current_helix[i].GetMomentumJacobian();
    const auto &c = m_current_helix[i].GetConstantOffset();

    E = - m_C * A.transpose() * G * B * S;
    V.block<3, 3>(0, 0) = m_C;
    V.block<3, 3>(3, 3) = S + E.transpose() * m_C_inv * E;
    V.block<3, 3>(3, 0) = E;
    V.block<3, 3>(0, 3) = E;

    r = centerAngle(measurement - c - A * m_vertex - B * p_s, measurement(2));
    m_chi2 += (r.transpose() * G * r)(0, 0);

  }

  m_ndf = 2 * m_numberOfDaughters - 3;

  if (m_use_ip_constraint) {
    m_ndf += 3;
  }

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

double FastFit::GetVariance(unsigned int component_i, unsigned int component_j) const {
    if(component_i >= 6) {
      throw std::runtime_error("Invalid index passed to GetDaughterMomentum");
    }
    if(component_j >= 6) {
      throw std::runtime_error("Invalid index passed to GetDaughterMomentum");
    }

    double variance = 0;
    if(component_i >= 3 && component_j >= 3) {
      variance = m_C(component_i, component_j);
    } else {
        for (unsigned int i = 0; i < m_numberOfDaughters; ++i) {
          variance += m_variances[i](component_i, component_j);
        }
    }
    return variance;
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
