#include <string>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <thread>
#include <chrono>
#include <mpi.h>
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

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // On exige au moins 2 processus : maître (rank 0) et esclave (rank 1)
    if (size < 2)
    {
        if (rank == 0)
            std::cerr << "Le programme doit être lancé avec au moins 2 processus." << std::endl;
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    // Récupération commune des paramètres
    auto params = parse_arguments(argc - 1, &argv[1]);

    // Variables pour accumuler le temps global de chaque itération
    double total_global_iter_time = 0.0;
    int nb_iter = 0;

    if (rank == 0)
    {
        // ------------------ Processus Maître : Affichage ------------------
        display_params(params);
        auto displayer = Displayer::init_instance(params.discretization, params.discretization);
        bool keep_running = true;

        while (keep_running)
        {
            // Vérification d'un message de terminaison envoyé par l'esclave (tag = 2)
            int flag = 0;
            MPI_Status status;
            MPI_Iprobe(1, 2, MPI_COMM_WORLD, &flag, &status);
            if (flag)
            {
                int term;
                MPI_Recv(&term, 1, MPI_INT, 1, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                keep_running = false;
                break;
            }

            double iter_start = MPI_Wtime();

            std::vector<std::uint8_t> vegetation(params.discretization * params.discretization);
            std::vector<std::uint8_t> fire(params.discretization * params.discretization);
            MPI_Request reqs[2];

            // Réception asynchrone des tableaux depuis l'esclave
            MPI_Irecv(vegetation.data(), vegetation.size(), MPI_UINT8_T, 1, 0, MPI_COMM_WORLD, &reqs[0]);
            MPI_Irecv(fire.data(), fire.size(), MPI_UINT8_T, 1, 1, MPI_COMM_WORLD, &reqs[1]);

            // Attente de la réception complète
            MPI_Waitall(2, reqs, MPI_STATUSES_IGNORE);

            // Mise à jour de l'affichage via SDL
            displayer->update(vegetation, fire);

            // Gestion d'événement pour permettre une fermeture manuelle
            SDL_Event event;
            if (SDL_PollEvent(&event) && event.type == SDL_QUIT)
            {
                // Envoi d'un message de terminaison à l'esclave si l'utilisateur ferme la fenêtre
                int term = 1;
                MPI_Send(&term, 1, MPI_INT, 1, 2, MPI_COMM_WORLD);
                keep_running = false;
            }

            double iter_end = MPI_Wtime();
            double local_iter_time = iter_end - iter_start;

            // Récupération du temps maximal de l'itération parmi les processus
            double global_iter_time = 0.0;
            MPI_Reduce(&local_iter_time, &global_iter_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

            total_global_iter_time += global_iter_time;
            nb_iter++;
        }

        if (nb_iter > 0)
        {
            std::cout << "Temps moyen global par itération (affichage, calcul et communication) : "
                      << (total_global_iter_time / nb_iter) << " s" << std::endl;
        }
    }
    else if (rank == 1)
    {
        // ------------------ Processus Esclave : Calcul ------------------
        if (!check_params(params))
        {
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        Model simu(params.length, params.discretization, params.wind, params.start);
        bool simulation_running = true;

        while (simulation_running)
        {
            double iter_start = MPI_Wtime();

            // Mise à jour de la simulation ; update retourne false si on a atteint la limite d'itérations
            simulation_running = simu.update();

            // Récupération des tableaux mis à jour
            auto vegetation = simu.vegetal_map();
            auto fire = simu.fire_map();

            MPI_Request reqs[2];
            // Envoi asynchrone des tableaux vers le maître
            MPI_Isend(vegetation.data(), vegetation.size(), MPI_UINT8_T, 0, 0, MPI_COMM_WORLD, &reqs[0]);
            MPI_Isend(fire.data(), fire.size(), MPI_UINT8_T, 0, 1, MPI_COMM_WORLD, &reqs[1]);

            // Attente de la complétion des envois regroupés
            MPI_Waitall(2, reqs, MPI_STATUSES_IGNORE);

            double iter_end = MPI_Wtime();
            double local_iter_time = iter_end - iter_start;

            // Réduction pour récupérer le maximum des temps locaux 
            double global_iter_time = 0.0;
            MPI_Reduce(&local_iter_time, &global_iter_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

            total_global_iter_time += local_iter_time;
            nb_iter++;
        }

        // Envoi d'un message de terminaison au maître (tag = 2) pour indiquer la fin de la simulation
        int term = 1;
        MPI_Send(&term, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);

    }

    MPI_Finalize();
    return EXIT_SUCCESS;
}
