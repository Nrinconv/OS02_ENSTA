#include <string>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <thread>
#include <chrono>
#include <mpi.h>
#include "model.hpp"
#include <stdexcept>
#include <fstream>
#include "display.hpp"
#include <stdexcept>
#include <cmath>
#include <iostream>
#include "model.hpp"
#include <stdexcept>
#include <cmath>
#include <iostream>
#include <fstream>
#include <condition_variable>
#include <iomanip>
#include <omp.h>
#include <chrono>
#include <mutex>
#include "model.hpp"
#include "display.hpp"

using namespace std::string_literals;
using namespace std::chrono_literals;

struct ParamsType
{
    double length{1.};
    unsigned discretization{20u};
    std::array<double,2> wind{0.,0.};
    Model::LexicoIndices start{10u,10u};
};

void analyze_arg( int nargs, char* args[], ParamsType& params )
{
    if (nargs ==0) return;
    std::string key(args[0]);
    if (key == "-l"s)
    {
        if (nargs < 2)
        {
            std::cerr << "Manque une valeur pour la longueur du terrain !" << std::endl;
            exit(EXIT_FAILURE);
        }
        params.length = std::stoul(args[1]);
        analyze_arg(nargs-2, &args[2], params);
        return;
    }
    auto pos = key.find("--longueur=");
    if (pos < key.size())
    {
        auto subkey = std::string(key,pos+11);
        params.length = std::stoul(subkey);
        analyze_arg(nargs-1, &args[1], params);
        return;
    }

    if (key == "-n"s)
    {
        if (nargs < 2)
        {
            std::cerr << "Manque une valeur pour le nombre de cases par direction pour la discrétisation du terrain !" << std::endl;
            exit(EXIT_FAILURE);
        }
        params.discretization = std::stoul(args[1]);
        analyze_arg(nargs-2, &args[2], params);
        return;
    }
    pos = key.find("--number_of_cases=");
    if (pos < key.size())
    {
        auto subkey = std::string(key, pos+18);
        params.discretization = std::stoul(subkey);
        analyze_arg(nargs-1, &args[1], params);
        return;
    }

    if (key == "-w"s)
    {
        if (nargs < 2)
        {
            std::cerr << "Manque une paire de valeurs pour la direction du vent !" << std::endl;
            exit(EXIT_FAILURE);
        }
        std::string values =std::string(args[1]);
        params.wind[0] = std::stod(values);
        auto pos = values.find(",");
        if (pos == values.size())
        {
            std::cerr << "Doit fournir deux valeurs séparées par une virgule pour définir la vitesse" << std::endl;
            exit(EXIT_FAILURE);
        }
        auto second_value = std::string(values, pos+1);
        params.wind[1] = std::stod(second_value);
        analyze_arg(nargs-2, &args[2], params);
        return;
    }
    pos = key.find("--wind=");
    if (pos < key.size())
    {
        auto subkey = std::string(key, pos+7);
        params.wind[0] = std::stoul(subkey);
        auto pos = subkey.find(",");
        if (pos == subkey.size())
        {
            std::cerr << "Doit fournir deux valeurs séparées par une virgule pour définir la vitesse" << std::endl;
            exit(EXIT_FAILURE);
        }
        auto second_value = std::string(subkey, pos+1);
        params.wind[1] = std::stod(second_value);
        analyze_arg(nargs-1, &args[1], params);
        return;
    }

    if (key == "-s"s)
    {
        if (nargs < 2)
        {
            std::cerr << "Manque une paire de valeurs pour la position du foyer initial !" << std::endl;
            exit(EXIT_FAILURE);
        }
        std::string values =std::string(args[1]);
        params.start.column = std::stod(values);
        auto pos = values.find(",");
        if (pos == values.size())
        {
            std::cerr << "Doit fournir deux valeurs séparées par une virgule pour définir la position du foyer initial" << std::endl;
            exit(EXIT_FAILURE);
        }
        auto second_value = std::string(values, pos+1);
        params.start.row = std::stod(second_value);
        analyze_arg(nargs-2, &args[2], params);
        return;
    }
    pos = key.find("--start=");
    if (pos < key.size())
    {
        auto subkey = std::string(key, pos+8);
        params.start.column = std::stoul(subkey);
        auto pos = subkey.find(",");
        if (pos == subkey.size())
        {
            std::cerr << "Doit fournir deux valeurs séparées par une virgule pour définir la vitesse" << std::endl;
            exit(EXIT_FAILURE);
        }
        auto second_value = std::string(subkey, pos+1);
        params.start.row = std::stod(second_value);
        analyze_arg(nargs-1, &args[1], params);
        return;
    }
}

ParamsType parse_arguments( int nargs, char* args[] )
{
    if (nargs == 0) return {};
    if ( (std::string(args[0]) == "--help"s) || (std::string(args[0]) == "-h") )
    {
        std::cout << 
R"RAW(Usage : simulation [option(s)]
  Lance la simulation d'incendie en prenant en compte les [option(s)].
  Les options sont :
    -l, --longueur=LONGUEUR     Définit la taille LONGUEUR (réel en km) du carré représentant la carte de la végétation.
    -n, --number_of_cases=N     Nombre n de cases par direction pour la discrétisation
    -w, --wind=VX,VY            Définit le vecteur vitesse du vent (pas de vent par défaut).
    -s, --start=COL,ROW         Définit les indices I,J de la case où commence l'incendie (milieu de la carte par défaut)
)RAW";
        exit(EXIT_SUCCESS);
    }
    ParamsType params;
    analyze_arg(nargs, args, params);
    return params;
}

bool check_params(ParamsType& params)
{
    bool flag = true;
    if (params.length <= 0)
    {
        std::cerr << "[ERREUR FATALE] La longueur du terrain doit être positive et non nulle !" << std::endl;
        flag = false;
    }

    if (params.discretization <= 0)
    {
        std::cerr << "[ERREUR FATALE] Le nombre de cellules par direction doit être positive et non nulle !" << std::endl;
        flag = false;
    }

    if ( (params.start.row >= params.discretization) || (params.start.column >= params.discretization) )
    {
        std::cerr << "[ERREUR FATALE] Mauvais indices pour la position initiale du foyer" << std::endl;
        flag = false;
    }
    
    return flag;
}

void display_params(ParamsType const& params)
{
    std::cout << "Parametres définis pour la simulation : \n"
              << "\tTaille du terrain : " << params.length << std::endl 
              << "\tNombre de cellules par direction : " << params.discretization << std::endl 
              << "\tVecteur vitesse : [" << params.wind[0] << ", " << params.wind[1] << "]" << std::endl
              << "\tPosition initiale du foyer (col, ligne) : " << params.start.column << ", " << params.start.row << std::endl;
}
// Fonction mise à jour pour les cellules fantômes avec réordonnement des appels
// Fonction mise à jour pour les cellules fantômes avec réordonnement des appels
void update_ghost_cells(std::vector<std::uint8_t>& local_fire,
                        std::vector<std::uint8_t>& local_vegetation,
                        int local_rows, int cols,
                        MPI_Comm newComm, int comp_rank, int comp_size) 
{
    // Échange pour le FEU ---------------------------------------------
    MPI_Request fire_requests[4];
    int fire_req_count = 0;
    
    // Réception Nord / Envoi Nord
    if (comp_rank > 0) {
        MPI_Irecv(local_fire.data(), cols, MPI_UINT8_T, comp_rank-1, 0, newComm, &fire_requests[fire_req_count++]);
        MPI_Isend(local_fire.data() + cols, cols, MPI_UINT8_T, comp_rank-1, 1, newComm, &fire_requests[fire_req_count++]);
    }

    // Réception Sud / Envoi Sud
    if (comp_rank < comp_size-1) {
        MPI_Irecv(local_fire.data() + (local_rows+1)*cols, cols, MPI_UINT8_T, comp_rank+1, 1, newComm, &fire_requests[fire_req_count++]);
        MPI_Isend(local_fire.data() + local_rows*cols, cols, MPI_UINT8_T, comp_rank+1, 0, newComm, &fire_requests[fire_req_count++]);
    }

    MPI_Waitall(fire_req_count, fire_requests, MPI_STATUSES_IGNORE);

    // Échange pour la VÉGÉTATION ---------------------------------------
    MPI_Request veg_requests[4];
    int veg_req_count = 0;
    
    // Réception Nord / Envoi Nord (tags différents: 2/3)
    if (comp_rank > 0) {
        MPI_Irecv(local_vegetation.data(), cols, MPI_UINT8_T, comp_rank-1, 2, newComm, &veg_requests[veg_req_count++]);
        MPI_Isend(local_vegetation.data() + cols, cols, MPI_UINT8_T, comp_rank-1, 3, newComm, &veg_requests[veg_req_count++]);
    }

    // Réception Sud / Envoi Sud
    if (comp_rank < comp_size-1) {
        MPI_Irecv(local_vegetation.data() + (local_rows+1)*cols, cols, MPI_UINT8_T, comp_rank+1, 3, newComm, &veg_requests[veg_req_count++]);
        MPI_Isend(local_vegetation.data() + local_rows*cols, cols, MPI_UINT8_T, comp_rank+1, 2, newComm, &veg_requests[veg_req_count++]);
    }

    MPI_Waitall(veg_req_count, veg_requests, MPI_STATUSES_IGNORE);

    // Conditions aux limites NEUMANN -----------------------------------
    // Pour le FEU
    if (comp_rank == 0) {
        std::copy(local_fire.begin() + cols, local_fire.begin() + 2*cols, local_fire.begin());
    }
    if (comp_rank == comp_size-1) {
        std::copy(
            local_fire.begin() + (local_rows) * cols,
            local_fire.begin() + (local_rows + 1) * cols,
            local_fire.begin() + (local_rows + 1) * cols
        );
    }

    // Pour la VÉGÉTATION
    if (comp_rank == 0) {
        std::copy(local_vegetation.begin() + cols, 
                 local_vegetation.begin() + 2*cols, 
                 local_vegetation.begin());
    }
    if (comp_rank == comp_size-1) {
        std::copy(local_vegetation.begin() + local_rows*cols,
                 local_vegetation.begin() + (local_rows+1)*cols,
                 local_vegetation.begin() + (local_rows+1)*cols);
    }
    MPI_Barrier(newComm);
}

void print_grid_state(const std::vector<std::uint8_t>& fire_map,
                      const std::vector<std::uint8_t>& veg_map,
                      unsigned discretization, int steps) 
{
    std::ofstream file("simulation_log.txt", std::ios::app);
    
    // Vérification cohérence taille des grilles
    const auto size = discretization * discretization;
    if (fire_map.size() != size || veg_map.size() != size) {
        file << "Erreur: Taille de grille incompatible avec la discrétisation\n";
        return;
    }

    // Affichage carte de feu
    file << "\n Step " << steps << " - Fire Map:\n";
    for (unsigned row = 0; row < discretization; ++row) {
        for (unsigned col = 0; col < discretization; ++col) {
            file << std::setw(4) << static_cast<int>(fire_map[row * discretization + col]);
        }
        file << "\n";
    }

    // Affichage carte végétation
    file << "\nVegetation Map:\n";
    for (unsigned row = 0; row < discretization; ++row) {
        for (unsigned col = 0; col < discretization; ++col) {
            file << std::setw(4) << static_cast<int>(veg_map[row * discretization + col]);
        }
        file << "\n";
    }
}
// Structure pour partager les données d'affichage entre threads
struct DisplayData {
    std::vector<std::uint8_t> vegetation;
    std::vector<std::uint8_t> fire;
    bool data_ready = false;
    bool keep_running = true;
    bool rendering_finished = false;
    double display_time = 0.0;
    int steps = 0;
};

// Fonction d'affichage dans un thread séparé
void rendering_thread_func(DisplayData& data, 
                          std::mutex& mtx, 
                          std::condition_variable& cv,
                          int discretization) {
    auto displayer = Displayer::init_instance(discretization, discretization);
    bool local_keep_running = true;

    while (local_keep_running) {
        std::unique_lock<std::mutex> lock(mtx);
        cv.wait(lock, [&data]{ return data.data_ready || !data.keep_running; });

        if (!data.keep_running) {
            local_keep_running = false;
            break;
        }

        // Copie locale des données pour minimiser le temps de verrouillage
        auto veg = data.vegetation;
        auto fire = data.fire;
        data.data_ready = false;
        lock.unlock();

        // Mesure du temps d'affichage
        auto start_display = std::chrono::high_resolution_clock::now();
        displayer->update(veg, fire);
        auto end_display = std::chrono::high_resolution_clock::now();
        double display_duration = std::chrono::duration<double, std::nano>(end_display - start_display).count() / 1e9; // Convertit en secondes

        // Mise à jour du temps d'affichage cumulé
        {
            std::lock_guard<std::mutex> lock(mtx);
            data.display_time += display_duration;
        }

        // Gestion des événements SDL
        SDL_Event event;
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) {
                std::lock_guard<std::mutex> lock(mtx);
                data.keep_running = false;
                local_keep_running = false;
            }
        }
    }

    std::lock_guard<std::mutex> lock(mtx);
    data.rendering_finished = true;
    cv.notify_one();
}
int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Comm globComm;
    MPI_Comm_dup(MPI_COMM_WORLD, &globComm);

    int world_rank, world_size;
    MPI_Comm_rank(globComm, &world_rank);
    MPI_Comm_size(globComm, &world_size);

    // Définition des tags MPI
    const int tag_signal = 0, tag_veg = 1, tag_fire = 2;
    // Nouveaux paramètres de configuration
    const int MAX_ITERATIONS = 2000;
    const int OUTPUT_INTERVAL = 100;

    // Variables de mesure des performances
    double total_ghost_time = 0.0;
    double total_compute_time = 0.0;
    double total_gather_time = 0.0;
    double total_display_time = 0.0;
    double total_comm_time = 0.0;
    int actual_iterations = 0;

    // Split du communicateur
    MPI_Comm newComm;
    int color = (world_rank == 0) ? 0 : 1;
    MPI_Comm_split(globComm, color, world_rank, &newComm);

    // Parsing des arguments
    ParamsType params = parse_arguments(argc - 1, argv + 1);
    if (!check_params(params)) {
        std::cerr << "[Global " << world_rank << "] Erreur dans les paramètres. Abandon." << std::endl;
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (world_rank == 0) {
        DisplayData display_data;
        std::mutex mtx;
        std::condition_variable cv;

        // Lancement du thread d'affichage
        std::thread rendering_thread(rendering_thread_func, 
                                    std::ref(display_data), 
                                    std::ref(mtx), 
                                    std::ref(cv),
                                    params.discretization);

        bool keep_running = true;
        int steps = 0;

        while (keep_running && steps < MAX_ITERATIONS) {
            std::vector<std::uint8_t> global_vegetation(params.discretization * params.discretization);
            std::vector<std::uint8_t> global_fire(params.discretization * params.discretization);

            // Réception non-bloquante des données
            MPI_Status status;
            int flag;
            MPI_Iprobe(MPI_ANY_SOURCE, tag_veg, globComm, &flag, &status);
            
            if (flag) {
                MPI_Recv(global_vegetation.data(), global_vegetation.size(), 
                        MPI_UINT8_T, MPI_ANY_SOURCE, tag_veg, globComm, MPI_STATUS_IGNORE);
                MPI_Recv(global_fire.data(), global_fire.size(), 
                        MPI_UINT8_T, MPI_ANY_SOURCE, tag_fire, globComm, MPI_STATUS_IGNORE);

                // Vérification des données reçues
                if (global_vegetation.empty() || global_fire.empty()) {
                    std::cerr << "Erreur : données reçues invalides." << std::endl;
                    continue;
                }

                // Mise à jour des données d'affichage
                {
                    std::lock_guard<std::mutex> lock(mtx);
                    display_data.vegetation = global_vegetation;
                    display_data.fire = global_fire;
                    display_data.data_ready = true;
                    display_data.steps = steps;
                }
                cv.notify_one();
                steps++;
            }

            // Vérification de l'arrêt
            {
                std::lock_guard<std::mutex> lock(mtx);
                if (!display_data.keep_running) {
                    keep_running = false;
                }
            }
        }

        // Signal d'arrêt au thread de rendu
        {
            std::lock_guard<std::mutex> lock(mtx);
            display_data.keep_running = false;
        }
        cv.notify_one();

        // Attente de la fin du thread
        rendering_thread.join();

        // Récupération du temps d'affichage cumulé
        double total_display_time = 0.0;
        int total_steps = 0;
        {
            std::lock_guard<std::mutex> lock(mtx);
            total_display_time = display_data.display_time;
            total_steps = display_data.steps;
        }

        // Affichage des statistiques
        if (total_steps > 0) {
            std::cout << "\n=== Statistiques d'affichage ==="
                      << "\nTemps total d'affichage: " << total_display_time << "s"
                      << "\nTemps moyen par frame: " << total_display_time / total_steps << "s"
                      << std::endl;
        } else {
            std::cerr << "Aucune itération n'a été effectuée." << std::endl;
        }

        // Envoi du signal d'arrêt aux autres processus
        int signal = -1;
        for (int i = 1; i < world_size; ++i) {
            MPI_Send(&signal, 1, MPI_INT, i, tag_signal, globComm);
        }
    } 
    else {
        // Processus de calcul
        int num_compute_procs = world_size - 1;  // Tous sauf le rang 0 (affichage)
        int total_rows = params.discretization;
        int cols = params.discretization;
        
        // Décomposition du domaine pour les processus de calcul (workers)
        int rows_per_proc = total_rows / num_compute_procs;
        int remainder = total_rows % num_compute_procs;
        
        // Calcul de la portion locale : comp_rank = world_rank - 1
        int comp_rank = world_rank - 1; 
        int local_rows = rows_per_proc + (comp_rank < remainder ? 1 : 0);
        int local_offset = 0;
        for (int i = 0; i < comp_rank; ++i) {
            local_offset += (i < remainder) ? rows_per_proc + 1 : rows_per_proc;
        }
        
        // Allocation de la grille locale avec 2 lignes fantômes
        std::vector<std::uint8_t> local_vegetation((local_rows + 2) * cols, 255);
        std::vector<std::uint8_t> local_fire((local_rows + 2) * cols, 0);
        
        // Le maître de calcul (comp_rank == 0, global rank 1) initialise la grille globale et distribue les tranches
        if (comp_rank == 0) {
            std::vector<std::uint8_t> global_fire(total_rows * cols, 0);
            std::vector<std::uint8_t> global_veg(total_rows * cols, 255);
            size_t fire_idx = params.start.row * cols + params.start.column;
            global_fire[fire_idx] = 255;  // Activation du feu initial

            // Traitement de la tranche locale pour le maître de calcul
            int proc_rows = rows_per_proc + (0 < remainder ? 1 : 0);
            for (int r = 0; r < proc_rows; ++r) {
                std::copy(
                    global_veg.begin() + r * cols,
                    global_veg.begin() + (r + 1) * cols,
                    local_vegetation.begin() + (r + 1) * cols
                );
                std::copy(
                    global_fire.begin() + r * cols,
                    global_fire.begin() + (r + 1) * cols,
                    local_fire.begin() + (r + 1) * cols
                );
            }
            std::copy(local_vegetation.begin() + cols, 
                    local_vegetation.begin() + 2*cols, 
                    local_vegetation.begin());
            std::copy(local_vegetation.begin() + (proc_rows)*cols, 
                    local_vegetation.begin() + (proc_rows + 1)*cols, 
                    local_vegetation.begin() + (proc_rows + 1)*cols);
            
            std::copy(local_fire.begin() + cols, 
                    local_fire.begin() + 2*cols, 
                    local_fire.begin());
            std::copy(local_fire.begin() + (proc_rows)*cols, 
                    local_fire.begin() + (proc_rows + 1)*cols, 
                    local_fire.begin() + (proc_rows + 1)*cols);

            // Envoi aux autres workers (processus de calcul dont comp_rank != 0)
            for (int i = 1; i < num_compute_procs; ++i) {
                int proc_rows = rows_per_proc + (i < remainder ? 1 : 0);
                int proc_offset = 0;
                for (int j = 0; j < i; ++j)
                    proc_offset += rows_per_proc + (j < remainder ? 1 : 0);

                // Buffer fire avec 2 lignes fantômes
                std::vector<std::uint8_t> buffer_fire((proc_rows + 2) * cols, 0);
                // Buffer végétation avec 2 lignes fantômes
                std::vector<std::uint8_t> buffer_veg((proc_rows + 2) * cols, 255);

                // Remplir les lignes réelles (entre les fantômes)
                for (int r = 0; r < proc_rows; ++r) {
                    std::copy(
                        global_fire.begin() + (proc_offset + r) * cols,
                        global_fire.begin() + (proc_offset + r + 1) * cols,
                        buffer_fire.begin() + (r + 1) * cols
                    );
                }
                // Remplir les lignes réelles (entre les fantômes)
                for (int r = 0; r < proc_rows; ++r) {
                    std::copy(
                        global_veg.begin() + (proc_offset + r) * cols,
                        global_veg.begin() + (proc_offset + r + 1) * cols,
                        buffer_veg.begin() + (r + 1) * cols
                    );
                }

                // Gestion des fantômes HAUT (sauf pour premier worker)
                if (i > 0) { 
                    std::copy(
                        global_fire.begin() + (proc_offset - 1) * cols,
                        global_fire.begin() + proc_offset * cols,
                        buffer_fire.begin()
                    );
                }

                // Gestion des fantômes BAS (sauf pour dernier worker)
                if (i < num_compute_procs - 1) {
                    std::copy(
                        global_fire.begin() + (proc_offset + proc_rows) * cols,
                        global_fire.begin() + (proc_offset + proc_rows + 1) * cols,
                        buffer_fire.end() - cols
                    );
                }
                // Gestion des fantômes HAUT (sauf pour premier worker)
                if (i > 0) { 
                    std::copy(
                        global_veg.begin() + (proc_offset - 1) * cols,
                        global_veg.begin() + proc_offset * cols,
                        buffer_veg.begin()
                    );
                }

                // Gestion des fantômes BAS (sauf pour dernier worker)
                if (i < num_compute_procs - 1) {
                    std::copy(
                        global_veg.begin() + (proc_offset + proc_rows) * cols,
                        global_veg.begin() + (proc_offset + proc_rows + 1) * cols,
                        buffer_veg.end() - cols
                    );
                }

                // Envoi du buffer fire
                MPI_Send(buffer_fire.data(), buffer_fire.size(), MPI_UINT8_T, i + 1, tag_fire, globComm);
                // Envoi du buffer végétation
                MPI_Send(buffer_veg.data(), buffer_veg.size(), MPI_UINT8_T, i + 1, tag_veg, globComm);
            }
        }
        
        // Les workers (comp_rank != 0) reçoivent la tranche initiale depuis le maître de calcul (global rank 1)
        if (comp_rank != 0) {
            MPI_Recv(local_fire.data(), local_fire.size(), MPI_UINT8_T, 1, tag_fire, globComm, MPI_STATUS_IGNORE);
            MPI_Recv(local_vegetation.data(), local_vegetation.size(), MPI_UINT8_T, 1, tag_veg, globComm, MPI_STATUS_IGNORE);
        }

        // Initialisation du modèle avec les données locales
        Model model(
            params.length, params.discretization, params.wind,
            params.start, local_rows, local_offset, cols,
            local_vegetation,  // Référence directe
            local_fire         // Référence directe
        );

        // Variables locales de timing
        std::chrono::duration<double> ghost_time{0};
        std::chrono::duration<double> compute_time{0};
        std::chrono::duration<double> gather_time{0};

        bool simulation_running = true;
        int iteration = 0;

        while (simulation_running && iteration < MAX_ITERATIONS) {
            auto iter_start = std::chrono::high_resolution_clock::now();

            // Mesure des communications fantômes
            auto start_ghost = std::chrono::high_resolution_clock::now();
            update_ghost_cells(local_fire, local_vegetation, local_rows, cols, 
                              newComm, comp_rank, num_compute_procs);
            auto end_ghost = std::chrono::high_resolution_clock::now();
            ghost_time += (end_ghost - start_ghost);

            // Mesure du calcul
            auto start_compute = std::chrono::high_resolution_clock::now();
            model.update_fire_front();
            bool local_continue = model.update();
            auto end_compute = std::chrono::high_resolution_clock::now();
            compute_time += (end_compute - start_compute);

            // Mesure du rassemblement des données
            auto start_gather = std::chrono::high_resolution_clock::now();
            std::vector<std::uint8_t> local_fire_data(local_fire.begin() + cols, local_fire.end() - cols);
            std::vector<std::uint8_t> local_veg_data(local_vegetation.begin() + cols, local_vegetation.end() - cols);

            std::vector<int> counts(num_compute_procs);
            std::vector<int> displs(num_compute_procs);
            int offset = 0;
            for (int i = 0; i < num_compute_procs; ++i) {
                counts[i] = ((i < remainder) ? rows_per_proc + 1 : rows_per_proc) * cols;
                displs[i] = offset;
                offset += counts[i];
            }

            std::vector<std::uint8_t> global_fire(total_rows * cols), global_veg(total_rows * cols);
            MPI_Gatherv(local_fire_data.data(), local_fire_data.size(), MPI_UINT8_T,
                       global_fire.data(), counts.data(), displs.data(), MPI_UINT8_T,
                       0, newComm);
            MPI_Gatherv(local_veg_data.data(), local_veg_data.size(), MPI_UINT8_T,
                       global_veg.data(), counts.data(), displs.data(), MPI_UINT8_T,
                       0, newComm);
            auto end_gather = std::chrono::high_resolution_clock::now();
            gather_time += (end_gather - start_gather);

            // Vérification de l'activité globale
            int local_active = local_continue ? 1 : 0;
            int global_active = 0;
            MPI_Allreduce(&local_active, &global_active, 1, MPI_INT, MPI_SUM, newComm);
            simulation_running = (global_active > 0);

            // Envoi au processus d'affichage
            if (comp_rank == 0) {
                MPI_Send(global_veg.data(), global_veg.size(), MPI_UINT8_T, 0, tag_veg, globComm);
                MPI_Send(global_fire.data(), global_fire.size(), MPI_UINT8_T, 0, tag_fire, globComm);
            }

            // Vérification arrêt
            int flag;
            MPI_Iprobe(0, tag_signal, globComm, &flag, MPI_STATUS_IGNORE);
            if (flag) {
                int signal;
                MPI_Recv(&signal, 1, MPI_INT, 0, tag_signal, globComm, MPI_STATUS_IGNORE);
                simulation_running = (signal != -1);
            }

            iteration++;
            actual_iterations = iteration;
        }

        // Collecte des statistiques locales
        double local_ghost = ghost_time.count();
        double local_compute = compute_time.count();
        double local_gather = gather_time.count();

        // Réduction MPI pour obtenir les totaux globaux
        double global_ghost, global_compute, global_gather;
        MPI_Reduce(&local_ghost, &global_ghost, 1, MPI_DOUBLE, MPI_MAX, 0, newComm);
        MPI_Reduce(&local_compute, &global_compute, 1, MPI_DOUBLE, MPI_MAX, 0, newComm);
        MPI_Reduce(&local_gather, &global_gather, 1, MPI_DOUBLE, MPI_MAX, 0, newComm);

        // Affichage des résultats sur le processus racine du calcul
        if (comp_rank == 0) {
            int num_procs = world_size - 1;
            std::cout << "\n=== Statistiques de calcul ==="
                      << "\nTemps moyen par iteration:"
                      << "\n- Communications fantômes: " << global_ghost/(num_procs*actual_iterations) << "s"
                      << "\n- Calcul modèle: " << global_compute/(num_procs*actual_iterations) << "s"
                      << "\n- Rassemblement données: " << global_gather/(num_procs*actual_iterations) << "s"
                      << std::endl;
        }
    }

    // Affichage final des stats de communication
    // if (world_rank == 0) {
    //     double avg_comm = total_comm_time / actual_iterations;
    //     std::cout << "\n=== Statistiques de communication ==="
    //               << "\nTemps moyen de communication globale: " << avg_comm << "s"
    //               << std::endl;
    // }

    MPI_Finalize();
    return EXIT_SUCCESS;
}