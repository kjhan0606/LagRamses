#define main weak_temperature_reference_cli_main
#include "weak_temperature_dsmc_reference.cpp"
#undef main

#include <fstream>
#include <iomanip>
#include <iostream>

int main(int argc, char** argv) {
  if (argc != 5) {
    std::cerr << "usage: test_common_ic_initializer ic_part N seed profile\n";
    return 2;
  }
  const std::filesystem::path path = argv[1];
  Params p;
  p.particles = static_cast<std::size_t>(std::stoull(argv[2]));
  p.seed = static_cast<std::uint64_t>(std::stoull(argv[3]));
  p.rho = 1.0;
  p.sigma_over_m = 1.0;
  p.temperature = 1.0;
  p.epsilon = 0.02;
  p.k_lambda = 0.25;
  p.energy_matched_uniform_control =
      std::string(argv[4]) == "energy_matched_uniform_control";
  if (!p.energy_matched_uniform_control && std::string(argv[4]) != "perturbed") {
    std::cerr << "invalid profile name\n";
    return 2;
  }

  const double lambda = 1.0 / (p.rho * p.sigma_over_m);
  const double wave_number = p.k_lambda / lambda;
  const double length = 2.0 * kPi / wave_number;
  Rng rng(p.seed);
  const auto expected = initialize(p, rng, length, wave_number);
  const double expected_mass = p.rho * length * length * length /
                               static_cast<double>(p.particles);
  const double x_limit = 8.0 * std::numeric_limits<double>::epsilon() * length;

  std::ifstream input(path);
  if (!input) {
    std::cerr << "cannot open common IC: " << path << '\n';
    return 2;
  }
  std::string line;
  std::size_t row_count = 0;
  double max_x_error = 0.0;
  double max_velocity_error = 0.0;
  while (std::getline(input, line)) {
    if (row_count >= expected.size()) {
      std::cerr << "common IC has too many rows\n";
      return 1;
    }
    std::istringstream row_stream(line);
    double row[7]{};
    if (!(row_stream >> row[0] >> row[1] >> row[2] >> row[3] >>
          row[4] >> row[5] >> row[6])) {
      std::cerr << "malformed row " << row_count + 1 << '\n';
      return 1;
    }
    std::string extra;
    if (row_stream >> extra) {
      std::cerr << "unexpected extra field on row " << row_count + 1 << '\n';
      return 1;
    }
    const double x_after_reader = row[0] + length / 2.0;
    const double x_error = std::abs(x_after_reader - expected[row_count].x);
    max_x_error = std::max(max_x_error, x_error);
    if (x_error > x_limit) {
      std::cerr << std::setprecision(17) << "x mismatch at row " << row_count + 1
                << ": error=" << x_error << " limit=" << x_limit << '\n';
      return 1;
    }
    const double expected_velocity[3] = {expected[row_count].v.x,
                                         expected[row_count].v.y,
                                         expected[row_count].v.z};
    for (int component = 0; component < 3; ++component) {
      const double error = std::abs(row[3 + component] - expected_velocity[component]);
      const double limit = 8.0 * std::numeric_limits<double>::epsilon() *
                           std::max(1.0, std::abs(expected_velocity[component]));
      max_velocity_error = std::max(max_velocity_error, error);
      if (error > limit) {
        std::cerr << std::setprecision(17) << "velocity draw mismatch at row "
                  << row_count + 1 << " component=" << component
                  << " error=" << error << " limit=" << limit << '\n';
        return 1;
      }
    }
    if (row[6] != expected_mass) {
      std::cerr << "mass mismatch at row " << row_count + 1 << '\n';
      return 1;
    }
    ++row_count;
  }
  if (row_count != expected.size()) {
    std::cerr << "common IC row count mismatch: got " << row_count
              << ", expected " << expected.size() << '\n';
    return 1;
  }
  std::cout << "COMMON_IC_WITHIN_8_EPS_OF_REFERENCE_INITIALIZER profile=" << argv[4]
            << " particles=" << row_count << " max_x_roundtrip_error="
            << std::setprecision(17) << max_x_error << " max_velocity_error="
            << max_velocity_error << '\n';
  return 0;
}
