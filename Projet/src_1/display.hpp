#pragma once
#include <vector>
#include <cstdint>
#include <memory>
#include <mutex>
#if defined(__linux__)
#  include <SDL2/SDL.h>
#endif

class Displayer
{
public:
    Displayer( std::uint32_t width, std::uint32_t height );
    Displayer( Displayer const & ) = delete;
    Displayer( Displayer      && ) = delete;
    ~Displayer();

    Displayer& operator = ( Displayer const & ) = delete;
    Displayer& operator = ( Displayer      && ) = delete;

    // Méthode existante : affichage complet
    void update( std::vector<std::uint8_t> const & vegetation_global_map,
                 std::vector<std::uint8_t> const & fire_global_map );

    // Nouvelle méthode : mise à jour d'une région (x, y, width, height)
    // Cette méthode dessine la région sans appeler SDL_RenderPresent
    void update_region( int x, int y, int width, int height,
                        std::vector<std::uint8_t> const & vegetation_sub_map,
                        std::vector<std::uint8_t> const & fire_sub_map );

    // Nouvelles méthodes pour encapsuler l'accès au renderer
    void clear_renderer();
    void present_renderer();

    static std::shared_ptr<Displayer> init_instance( std::uint32_t t_width, std::uint32_t t_height );
    static std::shared_ptr<Displayer> instance();

private:
    static std::shared_ptr<Displayer> unique_instance;

    SDL_Renderer *m_pt_renderer{nullptr};
    SDL_Surface  *m_pt_surface{nullptr};
    SDL_Window   *m_pt_window {nullptr};

    // Taille de la fenêtre pour le calcul des coordonnées
    int m_window_width{0};
    int m_window_height{0};

    // Mutex pour protéger l'accès au renderer (SDL n'étant pas thread-safe)
    std::mutex m_display_mutex;
};
