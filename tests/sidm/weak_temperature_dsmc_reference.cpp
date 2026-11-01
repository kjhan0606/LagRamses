// Independent 1D-in-space / 3D-in-velocity DSMC reference for a weak
// periodic temperature mode.  This file deliberately does not import or
// reproduce RAMSES pairing, probability, angular-scattering, or RNG code.
//
// Model: equal represented masses, constant sigma/m, isotropic elastic
// scattering, periodic x, and A=L^2 transverse area.  Within each frozen
// collision bin, a continuous-time null-collision process samples unordered
// pairs uniformly and accepts with g/g_max.  Streaming uses symmetric
// half-step / collision / half-step splitting.

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(__linux__)
#include <sys/resource.h>
#endif

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
using Complex = std::complex<double>;

struct Vec3 {
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
};

Vec3 operator+(Vec3 a, Vec3 b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
Vec3 operator-(Vec3 a, Vec3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
Vec3 operator*(double a, Vec3 b) { return {a * b.x, a * b.y, a * b.z}; }
Vec3 operator/(Vec3 b, double a) { return {b.x / a, b.y / a, b.z / a}; }
double dot(Vec3 a, Vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
double norm(Vec3 a) { return std::sqrt(dot(a, a)); }

struct Particle {
  double x = 0.0;
  Vec3 v{};
  double y = 0.0;
  double z = 0.0;
  std::uint64_t id = 0;
};

// SplitMix64 plus explicit open-interval uniforms and Box-Muller normals.
// This avoids implementation-dependent std::distribution algorithms.
class Rng {
 public:
  explicit Rng(std::uint64_t seed) : state_(seed) {}

  std::uint64_t integer() {
    std::uint64_t z = (state_ += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
  }

  double uniform() {
    const std::uint64_t q = integer() >> 13;  // 51 random bits
    return (static_cast<double>(q) + 0.5) * 0x1.0p-51;
  }

  std::size_t index(std::size_t n) {
    if (n == 0) throw std::runtime_error("cannot draw from an empty range");
    return static_cast<std::size_t>(uniform() * static_cast<double>(n));
  }

  double normal() {
    if (has_spare_) {
      has_spare_ = false;
      return spare_;
    }
    const double r = std::sqrt(-2.0 * std::log(uniform()));
    const double theta = 2.0 * kPi * uniform();
    spare_ = r * std::sin(theta);
    has_spare_ = true;
    return r * std::cos(theta);
  }

 private:
  std::uint64_t state_;
  bool has_spare_ = false;
  double spare_ = 0.0;
};

std::array<std::size_t, 2> choose_unordered_pair(std::size_t n, Rng& rng) {
  if (n < 2) throw std::runtime_error("pair selection needs at least two particles");
  const std::size_t first = rng.index(n);
  std::size_t second = rng.index(n - 1);
  if (second >= first) ++second;
  return {first, second};
}

struct Params {
  std::size_t particles = 4096;
  std::size_t cells = 64;
  std::size_t steps = 100;
  std::size_t sample_every = 10;
  double rho = 1.0;
  double sigma_over_m = 1.0;
  double temperature = 1.0;
  double epsilon = 0.01;
  double k_lambda = 0.25;
  double dt = std::numeric_limits<double>::quiet_NaN();
  double bulk_boost_x = 0.0;
  std::uint64_t seed = 20261027;
  bool collisions = true;
  bool energy_matched_uniform_control = false;
  bool export_endpoint_states = false;
  std::filesystem::path input_ic;
  std::filesystem::path output_dir;
};

void usage(std::ostream& out) {
  out << "Usage: weak_temperature_dsmc --output-dir PATH [options]\n"
      << "  --particles N       represented particles (default 4096)\n"
      << "  --cells NCELL       uniform x bins (default 64)\n"
      << "  --steps N           time steps (default 100)\n"
      << "  --sample-every N    output cadence in steps (default 10)\n"
      << "  --rho R             mean mass density (default 1)\n"
      << "  --sigma-over-m S    constant cross section per mass (default 1)\n"
      << "  --temperature T     T0 with m=kB=1 (default 1)\n"
      << "  --epsilon E         initial temperature-mode amplitude (default .01)\n"
      << "  --energy-matched-control  use uniform rho and T0*sqrt(1-epsilon^2)\n"
      << "  --k-lambda K        k/(rho*sigma/m) (default .25)\n"
      << "  --bulk-boost U      add uniform x velocity U (default 0)\n"
      << "  --dt DT             code-time step; default .005 collision times\n"
      << "  --seed N            SplitMix64 collision seed (default 20261027)\n"
      << "  --input-ic PATH     read the common seven-column RAMSES ASCII IC\n"
      << "  --export-endpoint-states  write ID-keyed t=0 and final states\n"
      << "  --collisionless     disable collisions for the exact free-stream control\n";
}

std::size_t parse_size(const std::string& text, const char* name) {
  if (text.empty() || text[0] == '-')
    throw std::runtime_error(std::string("invalid positive integer for ") + name);
  std::size_t used = 0;
  const unsigned long long value = std::stoull(text, &used);
  if (used != text.size() || value == 0ULL)
    throw std::runtime_error(std::string("invalid positive integer for ") + name);
  return static_cast<std::size_t>(value);
}

std::uint64_t parse_seed(const std::string& text) {
  if (text.empty() || text[0] == '-') throw std::runtime_error("invalid --seed value");
  std::size_t used = 0;
  const unsigned long long value = std::stoull(text, &used);
  if (used != text.size()) throw std::runtime_error("invalid --seed value");
  return static_cast<std::uint64_t>(value);
}

double parse_double(const std::string& text, const char* name) {
  std::size_t used = 0;
  const double value = std::stod(text, &used);
  if (used != text.size() || !std::isfinite(value))
    throw std::runtime_error(std::string("invalid finite value for ") + name);
  return value;
}

Params parse_args(int argc, char** argv) {
  Params p;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--help" || arg == "-h") {
      usage(std::cout);
      std::exit(0);
    }
    if (arg == "--collisionless") {
      p.collisions = false;
      continue;
    }
    if (arg == "--energy-matched-control") {
      p.energy_matched_uniform_control = true;
      continue;
    }
    if (arg == "--export-endpoint-states") {
      p.export_endpoint_states = true;
      continue;
    }
    if (i + 1 >= argc) throw std::runtime_error("missing value after " + arg);
    const std::string value = argv[++i];
    if (arg == "--output-dir") p.output_dir = value;
    else if (arg == "--input-ic") p.input_ic = value;
    else if (arg == "--particles") p.particles = parse_size(value, "--particles");
    else if (arg == "--cells") p.cells = parse_size(value, "--cells");
    else if (arg == "--steps") p.steps = parse_size(value, "--steps");
    else if (arg == "--sample-every") p.sample_every = parse_size(value, "--sample-every");
    else if (arg == "--rho") p.rho = parse_double(value, "--rho");
    else if (arg == "--sigma-over-m") p.sigma_over_m = parse_double(value, "--sigma-over-m");
    else if (arg == "--temperature") p.temperature = parse_double(value, "--temperature");
    else if (arg == "--epsilon") p.epsilon = parse_double(value, "--epsilon");
    else if (arg == "--k-lambda") p.k_lambda = parse_double(value, "--k-lambda");
    else if (arg == "--dt") p.dt = parse_double(value, "--dt");
    else if (arg == "--bulk-boost") p.bulk_boost_x = parse_double(value, "--bulk-boost");
    else if (arg == "--seed") p.seed = parse_seed(value);
    else throw std::runtime_error("unknown option: " + arg);
  }
  if (p.output_dir.empty()) throw std::runtime_error("--output-dir is required");
  if (p.particles < 2) throw std::runtime_error("at least two particles are required");
  if (!(p.rho > 0.0) || !(p.sigma_over_m > 0.0) || !(p.temperature > 0.0))
    throw std::runtime_error("rho, sigma-over-m, and temperature must be positive");
  if (!(std::abs(p.epsilon) < 1.0)) throw std::runtime_error("|epsilon| must be below one");
  if (!(p.k_lambda > 0.0)) throw std::runtime_error("k-lambda must be positive");
  if (!std::isfinite(p.bulk_boost_x * p.bulk_boost_x))
    throw std::runtime_error("bulk-boost squared must be finite");
  if (std::isfinite(p.dt) && !(p.dt > 0.0)) throw std::runtime_error("dt must be positive");
  if (p.export_endpoint_states && p.input_ic.empty())
    throw std::runtime_error("endpoint state export requires --input-ic");
  return p;
}

double wrap(double x, double length) {
  x -= length * std::floor(x / length);
  if (x >= length) x = 0.0;
  if (x < 0.0) x += length;
  return x;
}

Complex phase(double x, double k) { return std::polar(1.0, -k * x); }

double sinc(double x) {
  if (std::abs(x) < 1.0e-8) {
    const double x2 = x * x;
    return 1.0 - x2 / 6.0 + x2 * x2 / 120.0;
  }
  return std::sin(x) / x;
}

std::vector<Particle> initialize(const Params& p, Rng& rng, double length,
                                 double wave_number) {
  const double requested_norm = 1.0 / std::sqrt(1.0 - p.epsilon * p.epsilon);
  const double initial_epsilon = p.energy_matched_uniform_control ? 0.0 : p.epsilon;
  const double initial_temperature = p.energy_matched_uniform_control
      ? p.temperature / requested_norm : p.temperature;
  const double density_norm = 1.0 / std::sqrt(1.0 - initial_epsilon * initial_epsilon);
  std::vector<Particle> particles(p.particles);
  std::uint64_t next_id = 1;
  for (Particle& particle : particles) {
    particle.id = next_id++;
    double x = 0.0;
    double local_temperature = initial_temperature;
    for (;;) {
      x = length * rng.uniform();
      const double cosine = std::cos(wave_number * x);
      const double density_factor = 1.0 / (density_norm * (1.0 + initial_epsilon * cosine));
      const double acceptance = density_factor * density_norm *
                                (1.0 - std::abs(initial_epsilon));
      if (rng.uniform() <= acceptance) {
        local_temperature = initial_temperature * (1.0 + initial_epsilon * cosine);
        break;
      }
    }
    const double sigma_v = std::sqrt(local_temperature);
    particle.x = x;
    particle.v = {p.bulk_boost_x + sigma_v * rng.normal(), sigma_v * rng.normal(),
                  sigma_v * rng.normal()};
  }
  return particles;
}

std::vector<Particle> read_common_ic(const Params& p, double length,
                                     double represented_mass) {
  std::ifstream input(p.input_ic);
  if (!input) throw std::runtime_error("cannot open common IC: " + p.input_ic.string());
  std::vector<Particle> particles;
  particles.reserve(p.particles);
  std::string line;
  std::size_t row_count = 0;
  const double mass_tolerance = 8.0 * std::numeric_limits<double>::epsilon() *
                                std::abs(represented_mass);
  while (std::getline(input, line)) {
    if (row_count >= p.particles)
      throw std::runtime_error("common IC contains too many rows");
    std::istringstream row_stream(line);
    std::array<double, 7> row{};
    if (!(row_stream >> row[0] >> row[1] >> row[2] >> row[3] >>
          row[4] >> row[5] >> row[6]))
      throw std::runtime_error("malformed common IC row " + std::to_string(row_count + 1));
    std::string extra;
    if (row_stream >> extra)
      throw std::runtime_error("unexpected field in common IC row " +
                               std::to_string(row_count + 1));
    for (double value : row) {
      if (!std::isfinite(value))
        throw std::runtime_error("non-finite common IC row " +
                                 std::to_string(row_count + 1));
    }
    if (std::abs(row[6] - represented_mass) > mass_tolerance)
      throw std::runtime_error("common IC mass mismatch at row " +
                               std::to_string(row_count + 1));
    Particle particle;
    particle.id = static_cast<std::uint64_t>(row_count + 1);
    particle.x = wrap(row[0] + 0.5 * length, length);
    particle.y = wrap(row[1] + 0.5 * length, length);
    particle.z = wrap(row[2] + 0.5 * length, length);
    particle.v = {row[3], row[4], row[5]};
    particles.push_back(particle);
    ++row_count;
  }
  if (input.bad()) throw std::runtime_error("I/O error while reading common IC");
  if (row_count != p.particles)
    throw std::runtime_error("common IC row count does not match --particles");
  return particles;
}

struct CollisionStats {
  std::uint64_t candidates = 0;
  std::uint64_t accepted = 0;
  double maximum_momentum_residual = 0.0;
  double maximum_energy_residual = 0.0;
  double maximum_speed_bound_ratio = 0.0;
  Complex source_k{};
  double source_zero = 0.0;
};

void scatter_isotropic(Vec3& vi, Vec3& vj, Rng& rng) {
  const Vec3 center = 0.5 * (vi + vj);
  const Vec3 relative = vi - vj;
  const double speed = norm(relative);
  const double mu = 2.0 * rng.uniform() - 1.0;
  const double phi = 2.0 * kPi * rng.uniform();
  const double radial = std::sqrt(std::max(0.0, 1.0 - mu * mu));
  const Vec3 direction{radial * std::cos(phi), radial * std::sin(phi), mu};
  vi = center + 0.5 * speed * direction;
  vj = center - 0.5 * speed * direction;
}

CollisionStats collide(std::vector<Particle>& particles,
                       const std::vector<std::vector<std::size_t>>& bins,
                       const Params& p, Rng& rng, double dt, double length,
                       double wave_number, double represented_mass) {
  CollisionStats stats;
  if (!p.collisions || p.sigma_over_m == 0.0) return stats;
  const double volume = length * length * length;
  const double cell_volume = volume / static_cast<double>(p.cells);
  const double machine_margin = 1.0 + 1024.0 * std::numeric_limits<double>::epsilon();

  for (const auto& bin : bins) {
    const std::size_t n = bin.size();
    if (n < 2) continue;
    Vec3 mean{};
    for (std::size_t index : bin) mean = mean + particles[index].v;
    mean = mean / static_cast<double>(n);
    double thermal_speed_square_sum = 0.0;
    for (std::size_t index : bin) {
      const Vec3 c = particles[index].v - mean;
      thermal_speed_square_sum += dot(c, c);
    }
    if (!std::isfinite(thermal_speed_square_sum) || thermal_speed_square_sum < 0.0)
      throw std::runtime_error("nonfinite or invalid bin thermal-energy sum");
    // For any pair, |vi-vj|^2 <= 2 sum_i |vi-uc|^2.  Elastic
    // collisions preserve the bin's momentum and kinetic energy, so this
    // majorant remains valid throughout the frozen-position collision stage.
    const double speed_bound = std::sqrt(2.0 * thermal_speed_square_sum) * machine_margin;
    if (!std::isfinite(speed_bound) || speed_bound < 0.0)
      throw std::runtime_error("nonfinite or invalid conservative DSMC speed bound");
    if (speed_bound == 0.0) continue;
    const double pairs = 0.5 * static_cast<double>(n) * static_cast<double>(n - 1);
    const double candidate_rate = pairs * p.sigma_over_m * represented_mass /
                                  cell_volume * speed_bound;
    if (!(candidate_rate > 0.0) || !std::isfinite(candidate_rate))
      throw std::runtime_error("invalid null-collision candidate rate");

    double event_time = 0.0;
    while (true) {
      event_time += -std::log(rng.uniform()) / candidate_rate;
      if (event_time > dt) break;
      ++stats.candidates;
      const auto pair = choose_unordered_pair(n, rng);
      const std::size_t i = bin[pair[0]];
      const std::size_t j = bin[pair[1]];
      Vec3 vi = particles[i].v;
      Vec3 vj = particles[j].v;
      const Vec3 relative = vi - vj;
      const double relative_speed = norm(relative);
      const double ratio = relative_speed / speed_bound;
      stats.maximum_speed_bound_ratio = std::max(stats.maximum_speed_bound_ratio, ratio);
      if (!std::isfinite(ratio) || ratio > 1.0)
        throw std::runtime_error("relative speed exceeded conservative DSMC bound");
      if (rng.uniform() >= ratio) continue;

      const Vec3 old_momentum = vi + vj;
      const double old_energy = 0.5 * (dot(vi, vi) + dot(vj, vj));
      scatter_isotropic(vi, vj, rng);
      const Vec3 momentum_error = (vi + vj) - old_momentum;
      const double energy_error = 0.5 * (dot(vi, vi) + dot(vj, vj)) - old_energy;
      const double momentum_residual = norm(momentum_error);
      const double energy_residual = std::abs(energy_error);
      stats.maximum_momentum_residual =
          std::max(stats.maximum_momentum_residual, momentum_residual);
      stats.maximum_energy_residual =
          std::max(stats.maximum_energy_residual, energy_residual);
      const double local_scale = std::max(1.0, old_energy);
      if (momentum_residual > 512.0 * std::numeric_limits<double>::epsilon() *
                                  std::max(1.0, norm(old_momentum)) ||
          energy_residual > 1024.0 * std::numeric_limits<double>::epsilon() * local_scale)
        throw std::runtime_error("elastic pair update failed its conservation check");

      const double delta_i = 0.5 * represented_mass * (dot(vi, vi) - dot(particles[i].v, particles[i].v));
      const double delta_j = 0.5 * represented_mass * (dot(vj, vj) - dot(particles[j].v, particles[j].v));
      stats.source_k += (delta_i * phase(particles[i].x, wave_number) +
                         delta_j * phase(particles[j].x, wave_number)) / volume;
      stats.source_zero += (delta_i + delta_j) / volume;
      particles[i].v = vi;
      particles[j].v = vj;
      ++stats.accepted;
    }
  }
  return stats;
}

Complex exact_drift_flux_integral(const std::vector<Particle>& particles,
                                 double h, double length, double k,
                                 double represented_mass) {
  const double volume = length * length * length;
  Complex sum{};
  for (const Particle& particle : particles) {
    const double v2 = dot(particle.v, particle.v);
    const double half_phase = 0.5 * k * particle.v.x * h;
    const Complex integral = h * phase(particle.x, k) * sinc(half_phase) *
        std::polar(1.0, -half_phase);
    sum += (0.5 * represented_mass * v2 * particle.v.x / volume) * integral;
  }
  return sum;
}

void drift(std::vector<Particle>& particles, double h, double length) {
  for (Particle& particle : particles)
    particle.x = wrap(particle.x + particle.v.x * h, length);
}

Complex particle_energy_mode(const std::vector<Particle>& particles,
                             double length, double k, double represented_mass) {
  const double volume = length * length * length;
  Complex result{};
  for (const Particle& particle : particles)
    result += (0.5 * represented_mass * dot(particle.v, particle.v) / volume) *
              phase(particle.x, k);
  return result;
}

Complex free_stream_energy_mode(const std::vector<Particle>& initial,
                                double time, double length, double k,
                                double represented_mass) {
  const double volume = length * length * length;
  Complex result{};
  for (const Particle& particle : initial) {
    const double x = wrap(particle.x + particle.v.x * time, length);
    result += (0.5 * represented_mass * dot(particle.v, particle.v) / volume) *
              phase(x, k);
  }
  return result;
}

Complex continuum_free_stream_energy_mode(const Params& p, double time,
                                          double length, double k) {
  constexpr std::size_t quadrature = 16384;
  const double normalization = 1.0 / std::sqrt(1.0 - p.epsilon * p.epsilon);
  const double epsilon = p.energy_matched_uniform_control ? 0.0 : p.epsilon;
  const double temperature_0 = p.energy_matched_uniform_control
      ? p.temperature / normalization : p.temperature;
  const double a = k * time;
  Complex result{};
  for (std::size_t j = 0; j < quadrature; ++j) {
    const double x = length * (static_cast<double>(j) + 0.5) /
                     static_cast<double>(quadrature);
    const double temperature = temperature_0 * (1.0 + epsilon * std::cos(k * x));
    const double density_norm = 1.0 / std::sqrt(1.0 - epsilon * epsilon);
    const double density = p.rho / (density_norm * (1.0 + epsilon * std::cos(k * x)));
    const double exponential = std::exp(-0.5 * a * a * temperature);
    const Complex energy_factor = 0.5 *
        Complex(p.bulk_boost_x * p.bulk_boost_x + 3.0 * temperature -
                a * a * temperature * temperature,
                -2.0 * p.bulk_boost_x * a * temperature) * exponential *
        std::polar(1.0, -k * p.bulk_boost_x * time);
    result += density * energy_factor * phase(x, k) /
              static_cast<double>(quadrature);
  }
  return result;
}

struct Observables {
  Complex rho_particle{}, px_particle{}, py_particle{}, pz_particle{};
  Complex energy_particle{}, flux_particle{};
  Complex rho_cell{}, px_cell{}, py_cell{}, pz_cell{}, energy_cell{};
  Complex temperature_cell{}, thermal_energy_cell{};
  Complex qx_cell{}, flux_cell{}, flux_advective_cell{}, flux_pressure_cell{};
  double total_mass = 0.0;
  Vec3 total_momentum{};
  double total_kinetic_energy = 0.0;
  double max_flux_decomposition_error = 0.0;
  std::size_t empty_cells = 0;
  bool temperature_mode_complete = true;
  std::size_t min_occupancy = std::numeric_limits<std::size_t>::max();
  std::size_t max_occupancy = 0;
};

Observables measure(const std::vector<Particle>& particles,
                    const std::vector<std::vector<std::size_t>>& bins,
                    const Params& p, double length, double k,
                    double represented_mass) {
  const double volume = length * length * length;
  const double cell_volume = volume / static_cast<double>(p.cells);
  const double dx_over_length = 1.0 / static_cast<double>(p.cells);
  Observables o;
  for (const Particle& particle : particles) {
    const double v2 = dot(particle.v, particle.v);
    const Complex wphase = represented_mass * phase(particle.x, k) / volume;
    o.rho_particle += wphase;
    o.px_particle += wphase * particle.v.x;
    o.py_particle += wphase * particle.v.y;
    o.pz_particle += wphase * particle.v.z;
    o.energy_particle += wphase * (0.5 * v2);
    o.flux_particle += wphase * (0.5 * v2 * particle.v.x);
    o.total_momentum = o.total_momentum + represented_mass * particle.v;
    o.total_kinetic_energy += 0.5 * represented_mass * v2;
  }
  o.total_mass = represented_mass * static_cast<double>(particles.size());

  for (std::size_t c = 0; c < bins.size(); ++c) {
    const auto& bin = bins[c];
    const std::size_t n = bin.size();
    o.min_occupancy = std::min(o.min_occupancy, n);
    o.max_occupancy = std::max(o.max_occupancy, n);
    if (n == 0) {
      ++o.empty_cells;
      o.temperature_mode_complete = false;
      continue;
    }
    Vec3 momentum{};
    double energy = 0.0;
    double direct_flux = 0.0;
    for (std::size_t i : bin) {
      const Vec3 v = particles[i].v;
      momentum = momentum + represented_mass * v;
      energy += 0.5 * represented_mass * dot(v, v);
      direct_flux += 0.5 * represented_mass * dot(v, v) * v.x;
    }
    const double mass = represented_mass * static_cast<double>(n);
    const Vec3 u = momentum / mass;
    std::array<double, 9> pressure{};
    double thermal_energy = 0.0;
    double qx = 0.0;
    for (std::size_t i : bin) {
      const Vec3 cv = particles[i].v - u;
      const double c2 = dot(cv, cv);
      thermal_energy += 0.5 * represented_mass * c2;
      qx += 0.5 * represented_mass * c2 * cv.x;
      pressure[0] += represented_mass * cv.x * cv.x;
      pressure[1] += represented_mass * cv.x * cv.y;
      pressure[2] += represented_mass * cv.x * cv.z;
      pressure[3] += represented_mass * cv.y * cv.x;
      pressure[4] += represented_mass * cv.y * cv.y;
      pressure[5] += represented_mass * cv.y * cv.z;
      pressure[6] += represented_mass * cv.z * cv.x;
      pressure[7] += represented_mass * cv.z * cv.y;
      pressure[8] += represented_mass * cv.z * cv.z;
    }
    const double rho = mass / cell_volume;
    const double total_e_density = energy / cell_volume;
    const double e_th_density = thermal_energy / cell_volume;
    const double temp = (2.0 / 3.0) * thermal_energy / mass;
    const double pxx = pressure[0] / cell_volume;
    const double pxy = pressure[1] / cell_volume;
    const double pxz = pressure[2] / cell_volume;
    const double f_adv = total_e_density * u.x;
    const double f_pressure = pxx * u.x + pxy * u.y + pxz * u.z;
    const double f_q = qx / cell_volume;
    const double f_direct = direct_flux / cell_volume;
    const double f_decomposed = f_adv + f_pressure + f_q;
    o.max_flux_decomposition_error = std::max(
        o.max_flux_decomposition_error, std::abs(f_direct - f_decomposed));
    const Complex cell_phase = phase((static_cast<double>(c) + 0.5) *
                                     length / static_cast<double>(p.cells), k) *
                               dx_over_length;
    o.rho_cell += rho * cell_phase;
    o.px_cell += momentum.x / cell_volume * cell_phase;
    o.py_cell += momentum.y / cell_volume * cell_phase;
    o.pz_cell += momentum.z / cell_volume * cell_phase;
    o.energy_cell += total_e_density * cell_phase;
    o.temperature_cell += temp * cell_phase;
    o.thermal_energy_cell += e_th_density * cell_phase;
    o.qx_cell += f_q * cell_phase;
    o.flux_advective_cell += f_adv * cell_phase;
    o.flux_pressure_cell += f_pressure * cell_phase;
    o.flux_cell += f_decomposed * cell_phase;
  }
  if (o.min_occupancy == std::numeric_limits<std::size_t>::max()) o.min_occupancy = 0;
  if (!o.temperature_mode_complete) {
    const double undefined = std::numeric_limits<double>::quiet_NaN();
    o.temperature_cell = Complex(undefined, undefined);
  }
  return o;
}

void write_complex_header(std::ostream& out, const char* name) {
  out << ',' << name << "_re," << name << "_im";
}

void write_complex(std::ostream& out, Complex value) {
  out << ',' << value.real() << ',' << value.imag();
}

void write_csv_header(std::ostream& out) {
  out << "step,time,empty_cells,temperature_cell_mode_complete,min_occupancy,max_occupancy,total_mass,"
         "total_px,total_py,total_pz,total_kinetic_energy,"
         "collision_candidates,collision_accepted,acceptance_fraction,"
         "max_flux_decomposition_error,collision_pair_max_dp,collision_pair_max_dE,"
         "max_g_over_bound";
  for (const char* name : {"rho_particle_k", "px_particle_k", "py_particle_k",
                           "pz_particle_k", "energy_particle_k", "flux_particle_k",
                           "rho_cell_k", "px_cell_k", "py_cell_k", "pz_cell_k",
                           "energy_cell_k", "temperature_cell_k", "thermal_energy_cell_k",
                           "qx_cell_k", "flux_advective_cell_k", "flux_pressure_cell_k",
                           "flux_cell_k", "energy_free_particle_k",
                           "energy_free_continuum_k", "collision_source_k_cumulative",
                           "flux_integral_k_cumulative", "closure_residual_k"})
    write_complex_header(out, name);
  out << ",collision_source_zero_cumulative\n";
}

void write_row(std::ostream& out, std::size_t step, double time,
               const Observables& o, std::uint64_t candidates,
               std::uint64_t accepted, double max_dp, double max_dE,
               double max_g_ratio, Complex free_particle, Complex free_continuum,
               Complex source, double source_zero, Complex flux_integral,
               Complex initial_energy, double k) {
  const double acceptance = candidates == 0 ? 0.0 :
      static_cast<double>(accepted) / static_cast<double>(candidates);
  const Complex residual = o.energy_particle - initial_energy +
                           Complex(0.0, k) * flux_integral - source;
  out << step << ',' << time << ',' << o.empty_cells << ','
      << (o.temperature_mode_complete ? 1 : 0) << ',' << o.min_occupancy << ','
      << o.max_occupancy << ',' << o.total_mass << ',' << o.total_momentum.x << ','
      << o.total_momentum.y << ',' << o.total_momentum.z << ','
      << o.total_kinetic_energy << ',' << candidates << ',' << accepted << ','
      << acceptance << ',' << o.max_flux_decomposition_error << ',' << max_dp << ','
      << max_dE << ',' << max_g_ratio;
  for (Complex z : {o.rho_particle, o.px_particle, o.py_particle, o.pz_particle,
                    o.energy_particle, o.flux_particle, o.rho_cell, o.px_cell,
                    o.py_cell, o.pz_cell, o.energy_cell, o.temperature_cell,
                    o.thermal_energy_cell, o.qx_cell, o.flux_advective_cell,
                    o.flux_pressure_cell, o.flux_cell, free_particle, free_continuum,
                    source, flux_integral, residual})
    write_complex(out, z);
  out << ',' << source_zero << '\n';
}

void write_summary(const std::filesystem::path& path, const Params& p,
                   double length, double lambda, double k, double dt,
                   double collision_time, double represented_mass,
                   std::uint64_t candidates, std::uint64_t accepted,
                   const Observables& final, Complex final_residual,
                   double source_zero, double max_dp, double max_dE,
                   double max_g_ratio) {
  std::ofstream out(path);
  if (!out) throw std::runtime_error("cannot create summary.json");
  out << std::setprecision(17)
      << "{\n"
      << "  \"status\": \"REFERENCE_RUN_COMPLETE_NOT_PHYSICAL_ADMISSION\",\n"
      << "  \"model\": \"independent 1D3V equal-weight isotropic elastic DSMC\",\n"
      << "  \"rng\": \"SplitMix64 plus explicit Box-Muller\",\n"
      << "  \"sampler\": \"continuous-time null collision; uniform unordered pair; accept g/g_max\",\n"
      << "  \"particles\": " << p.particles << ",\n"
      << "  \"cells\": " << p.cells << ",\n"
      << "  \"steps\": " << p.steps << ",\n"
      << "  \"seed\": " << p.seed << ",\n"
      << "  \"collision_enabled\": " << (p.collisions ? "true" : "false") << ",\n"
      << "  \"rho_bar\": " << p.rho << ",\n"
      << "  \"sigma_over_m\": " << p.sigma_over_m << ",\n"
      << "  \"temperature_0\": " << p.temperature << ",\n"
      << "  \"epsilon\": " << p.epsilon << ",\n"
      << "  \"initial_condition\": \""
      << (p.energy_matched_uniform_control ? "energy_matched_uniform_control" :
          "isobaric_temperature_mode") << "\",\n"
      << "  \"effective_initial_epsilon\": "
      << (p.energy_matched_uniform_control ? 0.0 : p.epsilon) << ",\n"
      << "  \"effective_initial_temperature_0\": "
      << (p.energy_matched_uniform_control
              ? p.temperature * std::sqrt(1.0 - p.epsilon * p.epsilon)
              : p.temperature) << ",\n"
      << "  \"bulk_boost_x\": " << p.bulk_boost_x << ",\n"
      << "  \"k_lambda\": " << p.k_lambda << ",\n"
      << "  \"wave_number\": " << k << ",\n"
      << "  \"mean_free_path_convention\": \"lambda=1/(rho_bar*sigma_over_m)\",\n"
      << "  \"mean_free_path\": " << lambda << ",\n"
      << "  \"box_length\": " << length << ",\n"
      << "  \"effective_transverse_area\": " << length * length << ",\n"
      << "  \"represented_particle_mass\": " << represented_mass << ",\n"
      << "  \"cell_width_over_lambda\": " << length / static_cast<double>(p.cells) / lambda << ",\n"
      << "  \"collision_time_convention\": \"1/(rho*sigma*mean_relative_speed)\",\n"
      << "  \"collision_time\": " << collision_time << ",\n"
      << "  \"time_step\": " << dt << ",\n"
      << "  \"candidate_events\": " << candidates << ",\n"
      << "  \"accepted_events\": " << accepted << ",\n"
      << "  \"candidate_acceptance_fraction\": "
      << (candidates ? static_cast<double>(accepted) / static_cast<double>(candidates) : 0.0) << ",\n"
      << "  \"empty_cells_final\": " << final.empty_cells << ",\n"
      << "  \"temperature_cell_mode_complete_final\": "
      << (final.temperature_mode_complete ? "true" : "false") << ",\n"
      << "  \"total_mass_final\": " << final.total_mass << ",\n"
      << "  \"total_energy_final\": " << final.total_kinetic_energy << ",\n"
      << "  \"collision_source_zero_mode\": " << source_zero << ",\n"
      << "  \"closure_residual_k_abs\": " << std::abs(final_residual) << ",\n"
      << "  \"max_pair_momentum_residual\": " << max_dp << ",\n"
      << "  \"max_pair_energy_residual\": " << max_dE << ",\n"
      << "  \"max_observed_g_over_bound\": " << max_g_ratio << ",\n"
#if defined(__linux__)
      << "  \"maximum_rss_kib_linux\": " << []() {
           struct rusage usage {};
           if (getrusage(RUSAGE_SELF, &usage) != 0) return -1L;
           return usage.ru_maxrss;
         }() << ",\n"
#else
      << "  \"maximum_rss_kib_linux\": null,\n"
#endif
      << "  \"compiler\": \"" << __VERSION__ << "\",\n"
      << "  \"scientific_admission\": \"NOT_EVALUATED\"\n"
      << "}\n";
  out.flush();
  if (!out) throw std::runtime_error("failed flushing summary.json");
  out.close();
  if (!out) throw std::runtime_error("failed closing summary.json");
}

void write_endpoint_state(const Params& p, const std::vector<Particle>& particles,
                          int target, std::size_t step, double time,
                          double collision_time, double represented_mass) {
  std::ostringstream stem;
  stem << "sidm_kinetic_state_t" << std::setw(2) << std::setfill('0') << target
       << "_rank" << std::setw(4) << 1;
  const auto binary_path = p.output_dir / (stem.str() + ".bin");
  std::ofstream binary(binary_path, std::ios::binary | std::ios::out);
  if (!binary) throw std::runtime_error("cannot create endpoint state: " + binary_path.string());
  auto write_value = [&](const auto& value) {
    binary.write(reinterpret_cast<const char*>(&value),
                 static_cast<std::streamsize>(sizeof(value)));
    if (!binary) throw std::runtime_error("failed writing endpoint state");
  };
  for (const Particle& particle : particles) {
    const std::int64_t id = static_cast<std::int64_t>(particle.id);
    write_value(id);
    write_value(represented_mass);
    write_value(particle.x);
    write_value(particle.y);
    write_value(particle.z);
    write_value(particle.v.x);
    write_value(particle.v.y);
    write_value(particle.v.z);
  }
  binary.flush();
  if (!binary) throw std::runtime_error("failed flushing endpoint state");
  binary.close();
  if (!binary) throw std::runtime_error("failed closing endpoint state");

  const auto meta_path = p.output_dir / (stem.str() + ".meta");
  std::ofstream meta(meta_path);
  if (!meta) throw std::runtime_error("cannot create endpoint state metadata");
  const std::uint32_t endian_word = 1;
  const bool little_endian = *reinterpret_cast<const unsigned char*>(&endian_word) == 1;
  meta << std::setprecision(17)
       << "format_version=1\n"
       << "byte_order=" << (little_endian ? "little-endian" : "big-endian") << "\n"
       << "rank=1\n"
       << "target=" << target << "\n"
       << "step=" << step << "\n"
       << "time=" << time << "\n"
       << "tau0_fraction=" << time / collision_time << "\n"
       << "interval_exposure=0\n"
       << "local_count=" << particles.size() << "\n"
       << "id_bytes=8\nreal_bytes=8\nrecord_bytes=64\n"
       << "fields=int64_id,float64_mass,float64_x,float64_y,float64_z,"
       << "float64_vx,float64_vy,float64_vz\n";
  meta.flush();
  if (!meta) throw std::runtime_error("failed writing endpoint state metadata");
  meta.close();
  if (!meta) throw std::runtime_error("failed closing endpoint state metadata");

  const auto targets_path = p.output_dir / "sidm_kinetic_state_targets.csv";
  std::ofstream targets(targets_path, target == 0
      ? (std::ios::out | std::ios::trunc)
      : (std::ios::out | std::ios::app));
  if (!targets) throw std::runtime_error("cannot create state-target metadata");
  if (target == 0) targets << "target,step,time,tau0_fraction,N\n";
  targets << std::setprecision(17) << target << "," << step << "," << time
          << "," << time / collision_time << "," << particles.size() << "\n";
  targets.flush();
  if (!targets) throw std::runtime_error("failed writing state-target metadata");
}

int run(const Params& p) {
  const double lambda = 1.0 / (p.rho * p.sigma_over_m);
  const double mean_relative_speed = 4.0 * std::sqrt(p.temperature / kPi);
  const double collision_time = 1.0 /
      (p.rho * p.sigma_over_m * mean_relative_speed);
  const double wave_number = p.k_lambda / lambda;
  const double length = 2.0 * kPi / wave_number;
  const double dt = std::isfinite(p.dt) ? p.dt : 0.005 * collision_time;
  const double volume = length * length * length;
  const double cell_volume = volume / static_cast<double>(p.cells);
  const double represented_mass = p.rho * length * length * length /
                                  static_cast<double>(p.particles);
  const double final_time = static_cast<double>(p.steps) * dt;
  for (double value : {lambda, mean_relative_speed, collision_time, wave_number,
                       length, dt, volume, cell_volume, represented_mass, final_time}) {
    if (!std::isfinite(value) || value <= 0.0)
      throw std::runtime_error("derived geometry, timestep, or particle weight is nonfinite or nonpositive");
  }

  const std::filesystem::path parent = p.output_dir.parent_path();
  std::error_code ec;
  if (!parent.empty()) {
    std::filesystem::create_directories(parent, ec);
    if (ec) throw std::runtime_error("cannot create output parent: " + ec.message());
  }
  ec.clear();
  if (!std::filesystem::create_directory(p.output_dir, ec)) {
    if (ec) throw std::runtime_error("cannot atomically acquire output directory: " + ec.message());
    throw std::runtime_error("refusing to overwrite existing output directory: " +
                             p.output_dir.string());
  }

  Rng rng(p.seed);
  std::vector<Particle> particles = p.input_ic.empty()
      ? initialize(p, rng, length, wave_number)
      : read_common_ic(p, length, represented_mass);
  if (p.export_endpoint_states)
    write_endpoint_state(p, particles, 0, 0, 0.0, collision_time, represented_mass);
  const std::vector<Particle> initial = particles;
  std::vector<std::vector<std::size_t>> bins(p.cells);
  auto assign_bins = [&]() {
    for (auto& bin : bins) bin.clear();
    const double inv_dx = static_cast<double>(p.cells) / length;
    for (std::size_t i = 0; i < particles.size(); ++i) {
      std::size_t cell = static_cast<std::size_t>(particles[i].x * inv_dx);
      if (cell >= p.cells) cell = p.cells - 1;
      bins[cell].push_back(i);
    }
  };
  assign_bins();

  const auto csv_path = p.output_dir / "modes.csv";
  std::ofstream csv(csv_path);
  if (!csv) throw std::runtime_error("cannot create modes.csv");
  csv << std::setprecision(17);
  write_csv_header(csv);

  Complex source_k{};
  double source_zero = 0.0;
  Complex flux_integral{};
  std::uint64_t total_candidates = 0;
  std::uint64_t total_accepted = 0;
  double max_dp = 0.0;
  double max_dE = 0.0;
  double max_g_ratio = 0.0;
  double max_flux_decomposition_error = 0.0;

  auto write_sample = [&](std::size_t step) {
    assign_bins();
    const double time = static_cast<double>(step) * dt;
    const Observables o = measure(particles, bins, p, length, wave_number,
                                  represented_mass);
    max_flux_decomposition_error = std::max(max_flux_decomposition_error,
                                               o.max_flux_decomposition_error);
    const Complex initial_energy = particle_energy_mode(initial, length,
                                                         wave_number, represented_mass);
    const Complex free_particle = free_stream_energy_mode(initial, time, length,
                                                            wave_number, represented_mass);
    const Complex free_continuum = continuum_free_stream_energy_mode(
        p, time, length, wave_number);
    const Complex residual = o.energy_particle - initial_energy +
                             Complex(0.0, wave_number) * flux_integral - source_k;
    write_row(csv, step, time, o, total_candidates, total_accepted, max_dp,
              max_dE, max_g_ratio, free_particle, free_continuum, source_k,
              source_zero, flux_integral, initial_energy, wave_number);
    csv.flush();
    if (!csv) throw std::runtime_error("failed writing modes.csv");
    if (step == p.steps) {
      csv.close();
      if (!csv) throw std::runtime_error("failed closing modes.csv");
      write_summary(p.output_dir / "summary.json", p, length, lambda,
                    wave_number, dt, collision_time, represented_mass,
                    total_candidates, total_accepted, o, residual, source_zero,
                    max_dp, max_dE, max_g_ratio);
    }
  };

  write_sample(0);
  for (std::size_t step = 1; step <= p.steps; ++step) {
    flux_integral += exact_drift_flux_integral(particles, 0.5 * dt, length,
                                               wave_number, represented_mass);
    drift(particles, 0.5 * dt, length);
    assign_bins();
    CollisionStats stats = collide(particles, bins, p, rng, dt, length,
                                   wave_number, represented_mass);
    source_k += stats.source_k;
    source_zero += stats.source_zero;
    total_candidates += stats.candidates;
    total_accepted += stats.accepted;
    max_dp = std::max(max_dp, stats.maximum_momentum_residual);
    max_dE = std::max(max_dE, stats.maximum_energy_residual);
    max_g_ratio = std::max(max_g_ratio, stats.maximum_speed_bound_ratio);
    flux_integral += exact_drift_flux_integral(particles, 0.5 * dt, length,
                                               wave_number, represented_mass);
    drift(particles, 0.5 * dt, length);
    if (step % p.sample_every == 0 || step == p.steps) write_sample(step);
  }
  if (p.export_endpoint_states)
    write_endpoint_state(p, particles, 20, p.steps, final_time, collision_time,
                         represented_mass);
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Params params = parse_args(argc, argv);
    return run(params);
  } catch (const std::exception& error) {
    std::cerr << "weak_temperature_dsmc: " << error.what() << '\n';
    return 2;
  }
}
