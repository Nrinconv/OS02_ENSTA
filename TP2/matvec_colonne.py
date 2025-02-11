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

# Calcul de nombre de colonnes local
col_loc = dim // size

# Initialisation de la matrice local par blocs de colonnes
A_loc = np.array([[(i + j) % dim + 1 for i in range(rank * col_loc, (rank + 1) * col_loc)] for j in range(dim)], dtype=np.float64)
#print(f"A_loc = {A_loc}")
# Initialisation du vecteur u (partagé par tous les processus)
u = np.array([i + 1 for i in range(dim)], dtype=np.float64)
#print(f"u = {u}")

deb = time()

# Produit matrice-vecteur partiel
v_loc = np.dot(A_loc, u[rank * col_loc:(rank + 1) * col_loc])
#print(f"v_loc = {v_loc}")
# Réduction pour assembler le vecteur résultat complet
v = None
if rank == 0:
    v = np.zeros(dim, dtype=np.float64)
comm.Reduce(v_loc, v, op=MPI.SUM, root=0)

fin = time()

# Affichage du résultat final
if rank == 0:
    print(f"Temps total de calcul : {fin - deb:.4f}")
    print(f"Vecteur résultat v = {v}")
