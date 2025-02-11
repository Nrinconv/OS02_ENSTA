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
    escape_radius: float = 2.0

    def __contains__(self, c: complex) -> bool:
        return self.stability(c) == 1

    def convergence(self, c: complex, smooth=False, clamp=True) -> float:
        value = self.count_iterations(c, smooth)/self.max_iterations
        return max(0.0, min(value, 1.0)) if clamp else value

    def count_iterations(self, c: complex, smooth=False) -> int | float:
        if c.real*c.real + c.imag*c.imag < 0.0625:
            return self.max_iterations
        if (c.real+1)*(c.real+1) + c.imag*c.imag < 0.0625:
            return self.max_iterations
        if (c.real > -0.75) and (c.real < 0.5):
            ct = c.real-0.25 + 1.j * c.imag
            ctnrm2 = abs(ct)
            if ctnrm2 < 0.5*(1-ct.real/max(ctnrm2, 1.E-14)):
                return self.max_iterations
        
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

# quand il y a un seul processus
if size == 1:
    convergence = np.empty((width, height), dtype=np.double)
    deb = time()
    for y in range(height):
        for x in range(width):
            c = complex(-2.0 + scaleX * x, -1.125 + scaleY * y)
            convergence[x, y] = mandelbrot_set.convergence(c, smooth=True)
    fin = time()
    print(f"Temps total de calcul de Mandelbrot (séquentiel) : {fin - deb:.4f} secondes")
    deb = time()
    image = Image.fromarray(np.uint8(matplotlib.cm.plasma(convergence.T) * 255))
    fin = time()
    print(f"Temps de génération de l'image : {fin - deb:.4f} secondes")
    image.show()
    exit()

# Stratégie maître-esclave
if rank == 0:
    convergence = np.empty((width, height), dtype=np.double)
    lignes_restantes = list(range(height))
    processus_actifs = size - 1
    workers = list(range(1, size))
    deb = time()
    
    # Distribuer une première ligne à chaque processus
    for worker in workers:
        if lignes_restantes:
            ligne = lignes_restantes.pop(0)
            comm.send(ligne, dest=worker, tag=1)
        else:
            comm.send(None, dest=worker, tag=0)
            processus_actifs -= 1
    
    # Gerer la récupération des résultats et l'attribution des nouvelles tâches
    while processus_actifs > 0:
        worker, ligne, resultats = comm.recv(source=MPI.ANY_SOURCE, tag=2)
        convergence[:, ligne] = resultats
        
        if lignes_restantes:
            ligne = lignes_restantes.pop(0)
            comm.send(ligne, dest=worker, tag=1)
        else:
            comm.send(None, dest=worker, tag=0)
            processus_actifs -= 1
    
    fin = time()
    print(f"Temps total de calcul de Mandelbrot (maître-esclave) : {fin - deb:.4f} secondes")
    
    deb = time()
    image = Image.fromarray(np.uint8(matplotlib.cm.plasma(convergence.T) * 255))
    fin = time()
    print(f"Temps de génération de l'image : {fin - deb:.4f} secondes")
    image.show()

else:
    while True:
        ligne = comm.recv(source=0, tag=MPI.ANY_TAG)
        if ligne is None:
            break
        local_convergence = np.empty(width, dtype=np.double)
        for x in range(width):
            c = complex(-2.0 + scaleX * x, -1.125 + scaleY * ligne)
            local_convergence[x] = mandelbrot_set.convergence(c, smooth=True)
        comm.send((rank, ligne, local_convergence), dest=0, tag=2)
