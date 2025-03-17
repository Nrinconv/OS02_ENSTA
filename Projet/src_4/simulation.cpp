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
#include <iomanip>
#include <omp.h>
#include <chrono>
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

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Comm globComm;
    MPI_Comm_dup(MPI_COMM_WORLD, &globComm);

    int world_rank, world_size;
    MPI_Comm_rank(globComm, &world_rank);
    MPI_Comm_size(globComm, &world_size);

    std::cout << "[Global " << world_rank << "] Lancement de la simulation sur " << world_size << " processus." << std::endl;

    // Définition des tags MPI
    const int tag_signal = 0, tag_veg = 1, tag_fire = 2;

    // Split du communicateur : le processus de rang 0 sera dédié à l'affichage, les autres au calcul.
    MPI_Comm newComm;
    int color = (world_rank == 0) ? 0 : 1;
    MPI_Comm_split(globComm, color, world_rank, &newComm);
    std::cout << "[Global " << world_rank << "] Communicateur splitté (color " << color << ")." << std::endl;

    // Parsing des arguments
    ParamsType params = parse_arguments(argc - 1, argv + 1);
    if (!check_params(params)) {
        std::cerr << "[Global " << world_rank << "] Erreur dans les paramètres. Abandon." << std::endl;
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (world_rank == 0) {
        // Processus d'affichage
        std::cout << "[Global 0] Processus d'affichage lancé." << std::endl;
        display_params(params);
        auto displayer = Displayer::init_instance(params.discretization, params.discretization);
        bool keep_running = true;
        int steps = 0;
        while (keep_running) {
            std::vector<std::uint8_t> global_vegetation(params.discretization * params.discretization);
            std::vector<std::uint8_t> global_fire(params.discretization * params.discretization);
            
            std::cout << "[Global 0] Attente de la réception des grilles depuis un processus calcul." << std::endl;
            MPI_Recv(global_vegetation.data(), global_vegetation.size(), MPI_UINT8_T, MPI_ANY_SOURCE, tag_veg, globComm, MPI_STATUS_IGNORE);
            MPI_Recv(global_fire.data(), global_fire.size(), MPI_UINT8_T, MPI_ANY_SOURCE, tag_fire, globComm, MPI_STATUS_IGNORE);
            std::cout << "[Global 0] Grilles reçues. Mise à jour de l'affichage." << std::endl;
            
            displayer->update(global_vegetation, global_fire);
            if(steps == 100)
                print_grid_state(global_fire, global_vegetation, params.discretization, steps);
            steps++;

            SDL_Event event;
            if (SDL_PollEvent(&event) && event.type == SDL_QUIT) {
                std::cout << "[Global 0] Signal d'arrêt (SDL_QUIT) détecté." << std::endl;
                int signal = -1;
                for (int i = 1; i < world_size; ++i) {
                    MPI_Send(&signal, 1, MPI_INT, i, tag_signal, globComm);
                }
                keep_running = false;
            }
        }
        std::cout << "[Global 0] Processus d'affichage terminé." << std::endl;
    } 
    else {
        // Processus de calcul
        std::cout << "[Global " << world_rank << "] Processus de calcul lancé." << std::endl;
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
        std::cout << "[Global " << world_rank << " / Compute " << comp_rank 
                  << "] Portion locale : " << local_rows << " lignes, offset = " << local_offset << "." << std::endl;
        
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
            std::cout << "[Global " << world_rank << " / Compute 0] Tranche locale initialisée." << std::endl;

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
            std::cout << "[Global " << world_rank << " / Compute " << comp_rank << "] Réception des données initiales." << std::endl;
            MPI_Recv(local_fire.data(), local_fire.size(), MPI_UINT8_T, 1, tag_fire, globComm, MPI_STATUS_IGNORE);
            MPI_Recv(local_vegetation.data(), local_vegetation.size(), MPI_UINT8_T, 1, tag_veg, globComm, MPI_STATUS_IGNORE);
            std::cout << "[Global " << world_rank << " / Compute " << comp_rank << "] Données initiales reçues." << std::endl;
        } else {
            std::cout << "[Global " << world_rank << " / Compute 0] Utilisation des données initiales déjà disponibles localement." << std::endl;
        }

        // Initialisation du modèle avec les données locales
        Model model(
            params.length, params.discretization, params.wind,
            params.start, local_rows, local_offset, cols,
            local_vegetation,  // Référence directe
            local_fire         // Référence directe
        );
        std::cout << "[Global " << world_rank << " / Compute " << comp_rank << "] Modèle initialisé." << std::endl;

        // Boucle de simulation
        bool simulation_running = true;
        int iteration = 0;

        while (simulation_running) {
            // 1. Échange des cellules fantômes
            update_ghost_cells(local_fire, local_vegetation, local_rows, cols, newComm, comp_rank, num_compute_procs);
            std::cout << "[Global " << world_rank << " / Compute " << comp_rank << "] Échange des cellules fantômes terminé." << std::endl;

            //  Mettre à jour le front de feu avec les nouvelles cellules fantômes
            model.update_fire_front();
            // 2. Mise à jour du modèle
            bool local_continue = model.update();
            std::cout << "[Global " << world_rank << " / Compute " << comp_rank << "] Calcul de la nouvelle génération terminé." << std::endl;

            // 3. Extraction des données réelles (sans fantômes)
            std::vector<std::uint8_t> local_fire_data(local_fire.begin() + cols, local_fire.end() - cols);
            std::vector<std::uint8_t> local_veg_data(local_vegetation.begin() + cols, local_vegetation.end() - cols);
            
            // Vérifications de taille
            assert(local_fire.size() == (local_rows + 2) * cols);
            assert(local_vegetation.size() == (local_rows + 2) * cols);
            assert(local_fire_data.size() == local_rows * cols);
            assert(local_veg_data.size() == local_rows * cols);
            std::cout << "[Global " << world_rank << " / Compute " << comp_rank << "] Itération " << iteration << " début." << std::endl;

            // 4. Vérification de l'activité globale
            int local_active = local_continue ? 1 : 0;
            int global_active = 0;
            MPI_Allreduce(&local_active, &global_active, 1, MPI_INT, MPI_SUM, newComm);
            simulation_running = (global_active > 0);
            std::cout << "[Global " << world_rank << " / Compute " << comp_rank << "] AllReduce terminé. Feu actif: " << global_active << std::endl;

            // 5. Préparation données pour Gatherv
            // Calcul des counts/displs exacts
            std::vector<int> counts(num_compute_procs);
            std::vector<int> displs(num_compute_procs);
            int offset = 0;
            for (int i = 0; i < num_compute_procs; ++i) {
                counts[i] = ((i < remainder) ? rows_per_proc + 1 : rows_per_proc) * cols;
                displs[i] = offset;
                offset += counts[i];
            }

            // Vérification finale
            assert(offset == total_rows * cols);

            // Log des paramètres Gatherv
            std::cout << "[Global " << world_rank << " / Compute " << comp_rank << "] counts: ";
            for (auto c : counts) std::cout << c << " ";
            std::cout << "\nDispls: ";
            for (auto d : displs) std::cout << d << " ";
            std::cout << std::endl;

            // 6. Rassemblement des données
            std::vector<std::uint8_t> global_fire(total_rows * cols), global_veg(total_rows * cols);
            
            MPI_Gatherv(local_fire_data.data(), local_fire_data.size(), MPI_UINT8_T,
                        global_fire.data(), counts.data(), displs.data(), MPI_UINT8_T,
                        0, newComm);
            
            MPI_Gatherv(local_veg_data.data(), local_veg_data.size(), MPI_UINT8_T,
                        global_veg.data(), counts.data(), displs.data(), MPI_UINT8_T,
                        0, newComm);
            std::cout << "[Global " << world_rank << " / Compute " << comp_rank << "] Gatherv terminé." << std::endl;

            // 7. Envoi au processus d'affichage
            if (comp_rank == 0) {
                std::cout << "[Global " << world_rank << " / Compute 0] Envoi des données (" 
                        << global_veg.size() << " éléments) à l'affichage." << std::endl;
                MPI_Send(global_veg.data(), global_veg.size(), MPI_UINT8_T, 0, tag_veg, globComm);
                MPI_Send(global_fire.data(), global_fire.size(), MPI_UINT8_T, 0, tag_fire, globComm);
            }

            // 8. Vérification arrêt
            int flag;
            MPI_Iprobe(0, tag_signal, globComm, &flag, MPI_STATUS_IGNORE);
            if (flag) {
                int signal;
                MPI_Recv(&signal, 1, MPI_INT, 0, tag_signal, globComm, MPI_STATUS_IGNORE);
                simulation_running = (signal != -1);
                std::cout << "[Global " << world_rank << " / Compute " << comp_rank << "] Signal d'arrêt: " << signal << std::endl;
            }

            iteration++;
        }
        std::cout << "[Global " << world_rank << " / Compute " << comp_rank << "] Simulation terminée." << std::endl;
    }
    MPI_Finalize();
    std::cout << "[Global " << world_rank << "] Fin de la simulation." << std::endl;
    return EXIT_SUCCESS;
}