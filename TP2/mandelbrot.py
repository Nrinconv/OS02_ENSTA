# Calcul de l'ensemble de Mandelbrot en python
from mpi4py import MPI
import numpy as np
from dataclasses import dataclass
from PIL import Image
from math import log
from time import time
import matplotlib.cm


@dataclass
class MandelbrotSet:
    max_iterations: int
    escape_radius:  float = 2.0

    def __contains__(self, c: complex) -> bool:
        return self.stability(c) == 1

    def convergence(self, c: complex, smooth=False, clamp=True) -> float:
        value = self.count_iterations(c, smooth)/self.max_iterations
        return max(0.0, min(value, 1.0)) if clamp else value

    def count_iterations(self, c: complex,  smooth=False) -> int | float:
        z:    complex
        iter: int

        # On vérifie dans un premier temps si le complexe
        # n'appartient pas à une zone de convergence connue :
        #   1. Appartenance aux disques  C0{(0,0),1/4} et C1{(-1,0),1/4}
        if c.real*c.real+c.imag*c.imag < 0.0625:
            return self.max_iterations
        if (c.real+1)*(c.real+1)+c.imag*c.imag < 0.0625:
            return self.max_iterations
        #  2.  Appartenance à la cardioïde {(1/4,0),1/2(1-cos(theta))}
        if (c.real > -0.75) and (c.real < 0.5):
            ct = c.real-0.25 + 1.j * c.imag
            ctnrm2 = abs(ct)
            if ctnrm2 < 0.5*(1-ct.real/max(ctnrm2, 1.E-14)):
                return self.max_iterations
        # Sinon on itère
        z = 0
        for iter in range(self.max_iterations):
            z = z*z + c
            if abs(z) > self.escape_radius:
                if smooth:
                    return iter + 1 - log(log(abs(z)))/log(2)
                return iter
        return self.max_iterations

#MPI init
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

# On peut changer les paramètres des deux prochaines lignes
mandelbrot_set = MandelbrotSet(max_iterations=50, escape_radius=10)
width, height = 1024, 1024
scaleX = 3./width
scaleY = 2.25/height

# Répartition du travail entre les processus
lignes_par_proc = height // size
#premier_ligne = rank*lignes_par_proc
#dernier_ligne = (rank+1)*lignes_par_proc if rank != size-1 else height

local_convergence = np.empty((width, lignes_par_proc), dtype=np.double)

# Calcul de l'ensemble de mandelbrot :
deb = time()
""" Implementation question 1
for y in range(premier_ligne, dernier_ligne):
    for x in range(width):
        c = complex(-2. + scaleX*x, -1.125 + scaleY * y)
        local_convergence[x, y-premier_ligne] = mandelbrot_set.convergence(c, smooth=True)"""

#question2
for i, y in enumerate(range(rank, height, size)):  # Répartition intercalée
    for x in range(width):
        c = complex(-2.0 + scaleX * x, -1.125 + scaleY * y)
        local_convergence[x, i] = mandelbrot_set.convergence(c, smooth=True)


# Collecte des résultats sur le rank 0
if rank == 0:
    convergence = np.empty((width, height), dtype=np.double)
else:
    convergence = None

recvcounts = np.array([width * (height // size) for _ in range(size)])

# MPI.Gatherv (pour tailles différentes)
comm.Gatherv(local_convergence, [convergence, recvcounts], root=0)

# Calcul du temps total
total_time = comm.reduce(fin - deb, op=MPI.MAX, root=0)

if rank == 0:
    print(f"Temps du calcul de l'ensemble de Mandelbrot : {total_time:.4f}")
    deb = time()
    image = Image.fromarray(np.uint8(matplotlib.cm.plasma(convergence.T)*255))
    fin = time()
    print(f"Temps de constitution de l'image : {fin-deb}")
    image.show()

