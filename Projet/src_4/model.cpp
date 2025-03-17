#include <stdexcept>
#include <cmath>
#include <iostream>
#include "model.hpp"
#include <stdexcept>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <omp.h>
#include <chrono>
#include "model.hpp"
#include "display.hpp"


namespace
{
    double pseudo_random( std::size_t index, std::size_t time_step )
    {
        std::uint_fast32_t xi = std::uint_fast32_t(index*(time_step+1));
        std::uint_fast32_t r  = (48271*xi)%2147483647;
        return r/2147483646.;
    }

    double log_factor( std::uint8_t value )
    {
        return std::log(1.+value)/std::log(256);
    }
}
Model::Model(double t_length, unsigned t_discretization, std::array<double,2> t_wind,
             LexicoIndices t_start_fire_position,
             int local_rows, int local_offset, int cols,
             std::vector<std::uint8_t>& vegetation_vec, std::vector<std::uint8_t>& fire_vec,
             double t_max_wind)
    : m_local_rows(local_rows),
      m_local_offset(local_offset),
      m_length(t_length),
      m_distance(-1),
      m_geometry(t_discretization),
      m_wind(t_wind),
      m_wind_speed(std::sqrt(t_wind[0] * t_wind[0] + t_wind[1] * t_wind[1])),
      m_max_wind(t_max_wind),
      m_vegetation_map(vegetation_vec),
      m_fire_map(fire_vec)
{
    if (t_discretization == 0) {
        throw std::range_error("Le nombre de cases par direction doit être plus grand que zéro.");
    }
    m_distance = m_length / double(m_geometry);
    // Initialiser m_fire_front avec les cellules en feu locales
    for (size_t idx = 0; idx < m_fire_map.size(); ++idx) {
        if (m_fire_map[idx] == 255u) {
            m_fire_front[idx] = 255u;
        }
    }

    // Initialisation des probabilités et des coefficients de vent
    constexpr double alpha0 = 4.52790762e-01;
    constexpr double alpha1 = 9.58264437e-04;
    constexpr double alpha2 = 3.61499382e-05;

    if (m_wind_speed < t_max_wind)
        p1 = alpha0 + alpha1 * m_wind_speed + alpha2 * (m_wind_speed * m_wind_speed);
    else
        p1 = alpha0 + alpha1 * t_max_wind + alpha2 * (t_max_wind * t_max_wind);
    p2 = 0.3;

    if (m_wind[0] > 0) {
        alphaEastWest = std::abs(m_wind[0] / t_max_wind) + 1;
        alphaWestEast = 1. - std::abs(m_wind[0] / t_max_wind);
    } else {
        alphaWestEast = std::abs(m_wind[0] / t_max_wind) + 1;
        alphaEastWest = 1. - std::abs(m_wind[0] / t_max_wind);
    }

    if (m_wind[1] > 0) {
        alphaSouthNorth = std::abs(m_wind[1] / t_max_wind) + 1;
        alphaNorthSouth = 1. - std::abs(m_wind[1] / t_max_wind);
    } else {
        alphaNorthSouth = std::abs(m_wind[1] / t_max_wind) + 1;
        alphaSouthNorth = 1. - std::abs(m_wind[1] / t_max_wind);
    }
}

void Model::update_fire_front() {
    m_fire_front.clear();
    for (size_t idx = 0; idx < m_fire_map.size(); ++idx) {
        if (m_fire_map[idx] > 0u) {  
            m_fire_front[idx] = m_fire_map[idx]; 
        }
    }
}
// --------------------------------------------------------------------------------------------------------------------
bool Model::update()
{

    auto next_front = m_fire_front;
    
    for (auto& f : m_fire_front)
    {
        LexicoIndices coord = get_lexicographic_from_index(f.first);
        int global_row = (coord.row -1) + m_local_offset;


        double power = log_factor(f.second);

        // Propagation Nord
        if (global_row > 0) {
            size_t neighbor_idx = f.first - m_geometry;
            if (neighbor_idx < m_fire_map.size()) {
                double tirage = pseudo_random(neighbor_idx * 13427 + m_time_step, m_time_step);
                double green_power = m_vegetation_map[neighbor_idx];
                if (tirage < alphaNorthSouth * p1 * power * log_factor(green_power)) {
                    m_fire_map[neighbor_idx] = 255;
                    next_front[neighbor_idx] = 255;
                }
            }
        }

        // Propagation Sud
        if (global_row < m_geometry - 1) {
            size_t neighbor_idx = f.first + m_geometry;
            if (neighbor_idx < m_fire_map.size()) {
                double tirage = pseudo_random(neighbor_idx * 75329 + m_time_step, m_time_step);
                double green_power = m_vegetation_map[neighbor_idx];
                if (tirage < alphaSouthNorth * p1 * power * log_factor(green_power)) {
                    m_fire_map[neighbor_idx] = 255;
                    next_front[neighbor_idx] = 255;
                }
            }
        }

        // Propagation Ouest
        if (coord.column > 0) {
            size_t neighbor_idx = f.first - 1;
            if (neighbor_idx < m_fire_map.size()) {
                double tirage = pseudo_random(neighbor_idx * 45673 + m_time_step, m_time_step);
                double green_power = m_vegetation_map[neighbor_idx];
                if (tirage < alphaWestEast * p1 * power * log_factor(green_power)) {
                    m_fire_map[neighbor_idx] = 255;
                    next_front[neighbor_idx] = 255;
                }
            }
        }

        // Propagation Est
        if (coord.column < m_geometry - 1) {
            size_t neighbor_idx = f.first + 1;
            if (neighbor_idx < m_fire_map.size()) {
                double tirage = pseudo_random(neighbor_idx * 29123 + m_time_step, m_time_step);
                double green_power = m_vegetation_map[neighbor_idx];
                if (tirage < alphaEastWest * p1 * power * log_factor(green_power)) {
                    m_fire_map[neighbor_idx] = 255;
                    next_front[neighbor_idx] = 255;
                }
            }
        }
        // Si le feu est à son max,
        if (f.second == 255)
        {   // On regarde si il commence à faiblir pour s'éteindre au bout d'un moment :
            double tirage = pseudo_random( f.first * 52513 + m_time_step, m_time_step);
            if (tirage < p2)
            {
                m_fire_map[f.first] >>= 1;
                next_front[f.first] >>= 1;
            }
        }
        else
        {
            // Foyer en train de s'éteindre.
            m_fire_map[f.first] >>= 1;
            next_front[f.first] >>= 1;
            if (next_front[f.first] == 0)
            {
                next_front.erase(f.first);
            }
        }

    }    
    for (auto f : m_fire_front)
    {
        if (m_vegetation_map[f.first] > 0)
            m_vegetation_map[f.first] -= 1;
    }

    m_time_step += 1;
    m_fire_front = next_front;  

    return !m_fire_front.empty();
}

std::size_t Model::get_index_from_lexicographic_indices(LexicoIndices t_lexico) const {
    // Convertit les indices locaux en indices dans le tableau local incluant les fantômes
    return (t_lexico.row + 1) * m_geometry + t_lexico.column;
}

Model::LexicoIndices Model::get_lexicographic_from_index(std::size_t t_local_index) const {
    LexicoIndices ind_coords;
    ind_coords.row = (t_local_index / m_geometry); 
    ind_coords.column = t_local_index % m_geometry;
    return ind_coords;
}
void Model::log_grids(std::size_t step) const {
    std::ofstream file("simulation_log.txt", std::ios::app);  // Mode append
    file << "Step " << step << " - Fire Map:\n";
    for (std::size_t i = 0; i < m_geometry; ++i) {
        for (std::size_t j = 0; j < m_geometry; ++j) {
            file << std::setw(4) << static_cast<int>(m_fire_map[i * m_geometry + j]);
        }
        file << "\n";
    }
    file << "Vegetation Map:\n";
    for (std::size_t i = 0; i < m_geometry; ++i) {
        for (std::size_t j = 0; j < m_geometry; ++j) {
            file << std::setw(4) << static_cast<int>(m_vegetation_map[i * m_geometry + j]);
        }
        file << "\n";
    }
    file << "\n";
}
