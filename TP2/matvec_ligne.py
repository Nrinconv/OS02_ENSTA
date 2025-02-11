# Produit matrice-vecteur v = A.u
from mpi4py import MPI
import numpy as np
from time import time


# MPI init
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# Dimension du problème (peut-être changé)
dim = 120

# Calcul de nombre de lignes local
lign_loc = dim // size

# Initialisation de la matrice local par blocs de colonnes
A_loc = np.array([[(i + j) % dim + 1 for i in range(dim)] for j in range(rank * lign_loc, (rank + 1) * lign_loc)], dtype=np.float64)
#print(f"A_loc = {A_loc}")
# Initialisation du vecteur u (partagé par tous les processus)
u = np.array([i + 1 for i in range(dim)], dtype=np.float64)
#print(f"u = {u}")

deb = time()

# Produit matrice-vecteur partiel
v_loc = np.dot(A_loc, u)
#print(f"v_loc = {v_loc}")
# Réduction pour assembler le vecteur résultat complet
v = None
if rank == 0:
    v = np.zeros(dim, dtype=np.float64)
    recvcounts = np.array([lign_loc] * size, dtype=np.int64)
else:
    recvcounts = None

comm.Gatherv(v_loc, [v, recvcounts], root=0)

fin = time()

# Affichage du résultat final
if rank == 0:
    print(f"Temps total de calcul : {fin - deb:.4f}")
    print(f"Vecteur résultat v = {v}")
