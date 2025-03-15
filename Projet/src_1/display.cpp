#include <cassert>
#include <stdexcept>
#include <string> 
#include "display.hpp"

using namespace std::string_literals;

std::shared_ptr<Displayer> Displayer::unique_instance{nullptr};

Displayer::Displayer( std::uint32_t t_width, std::uint32_t t_height )
{
    if ( SDL_Init( SDL_INIT_EVERYTHING) < 0 )
    {
        std::string err_msg = "Erreur lors de l'initialisation de SDL : "s + std::string(SDL_GetError());
        throw std::runtime_error(err_msg);
    }
    SDL_CreateWindowAndRenderer(t_width, t_height, 0, &m_pt_window, &m_pt_renderer);
    if (m_pt_window == nullptr)
    {
        std::string err_msg = "Erreur lors de la création de la fenêtre : "s + std::string(SDL_GetError());
        throw std::runtime_error(err_msg);
    }
    if (m_pt_renderer == nullptr)
    {
        std::string err_msg = "Erreur lors de la création du moteur de rendu : "s + std::string(SDL_GetError());
        throw std::runtime_error(err_msg);
    }
    m_pt_surface = SDL_GetWindowSurface( m_pt_window );
    if (m_pt_surface == nullptr)
    {
        std::string err_msg = "Erreur lors de la récupération de la surface : "s + std::string(SDL_GetError());
        throw std::runtime_error(err_msg);
    }
    SDL_GetWindowSize(m_pt_window, &m_window_width, &m_window_height);
}

Displayer::~Displayer()
{
    SDL_DestroyRenderer(m_pt_renderer);
    SDL_DestroyWindow( m_pt_window );
    SDL_Quit();
}

void
Displayer::update( std::vector<std::uint8_t> const & vegetation_global_map,
                   std::vector<std::uint8_t> const & fire_global_map )
{
    int w, h;
    SDL_GetWindowSize(m_pt_window, &w, &h );
    SDL_SetRenderDrawColor(m_pt_renderer, 0, 0, 0, 255);
    SDL_RenderClear(m_pt_renderer);
    for (int i = 0; i < h; ++i )
      for (int j = 0; j < w; ++j )
      {
        SDL_SetRenderDrawColor(m_pt_renderer, fire_global_map[j + w*i],
                               vegetation_global_map[j + w*i], 0, 255);
        SDL_RenderDrawPoint(m_pt_renderer, j, h-i-1);
      }
    SDL_RenderPresent(m_pt_renderer);
}

void Displayer::update_region(int x, int y, int width, int height,
                              std::vector<std::uint8_t> const & vegetation_sub_map,
                              std::vector<std::uint8_t> const & fire_sub_map)
{
    // Accès protégé au renderer
    std::lock_guard<std::mutex> lock(m_display_mutex);
    // Pour chaque pixel de la région, calculer sa position globale
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            // La position globale en x reste j (puisque x vaut toujours 0 ici)
            // La position globale en y est (y + i).
            // On convertit en coordonnées SDL avec l'origine en haut à gauche
            int global_y = y + i;
            SDL_SetRenderDrawColor(m_pt_renderer,
                                   fire_sub_map[j + width * i],
                                   vegetation_sub_map[j + width * i],
                                   0, 255);
            SDL_RenderDrawPoint(m_pt_renderer, x + j, m_window_height - global_y - 1);
        }
    }
    // Ne pas appeler SDL_RenderPresent ici
}


void Displayer::clear_renderer()
{
    std::lock_guard<std::mutex> lock(m_display_mutex);
    SDL_SetRenderDrawColor(m_pt_renderer, 0, 0, 0, 255);
    SDL_RenderClear(m_pt_renderer);
}

void Displayer::present_renderer()
{
    std::lock_guard<std::mutex> lock(m_display_mutex);
    SDL_RenderPresent(m_pt_renderer);
}

std::shared_ptr<Displayer> 
Displayer::init_instance( std::uint32_t t_width, std::uint32_t t_height )
{
    assert( ( "L'initialisation de l'instance ne doit etre appele qu'une seule fois !" && (unique_instance == nullptr) ) );
    unique_instance = std::make_shared<Displayer>(t_width, t_height);
    return unique_instance;
}

std::shared_ptr<Displayer> 
Displayer::instance()
{
    assert( ( "Il faut initialiser l'instance avant tout !" && (unique_instance != nullptr) ) );
    return unique_instance;
}
