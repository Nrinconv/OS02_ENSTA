""" from PIL import Image
import os
import numpy as np
from scipy import signal
import time

# fonction pour doubler la taille d'une image sans trop la pixeliser
def double_size(image):
    # démarrer le chronomètre
    start_time = time.time()
    
    # on charge l'image
    img = Image.open(image)
    print(f"original size: {img.size}")
    
    # conversion en hsv
    img = img.convert('HSV')
    
    # on convertit l'image en tableau numpy
    img = np.array(img, dtype=np.double)
    
    # on double sa taille et on la normalise
    img = np.repeat(np.repeat(img, 2, axis=0), 2, axis=1) / 255.
    print(f"new size: {img.shape}")
    
    # on crée un masque de flou gaussien
    mask = np.array([[1., 2., 1.], [2., 4., 2.], [1., 2., 1.]]) / 16.
    
    # on applique le filtre de flou
    blur_image = np.zeros_like(img, dtype=np.double)
    for i in range(3):
        blur_image[:, :, i] = signal.convolve2d(img[:, :, i], mask, mode='same')
    
    # on crée un masque de netteté
    mask = np.array([[0, -1, 0], [-1, 5, -1], [0, -1, 0]])
    
    # on applique le filtre de netteté uniquement sur la luminance
    sharpen_image = np.zeros_like(img, dtype=np.double)
    sharpen_image[:, :, :2] = blur_image[:, :, :2]
    sharpen_image[:, :, 2] = np.clip(signal.convolve2d(blur_image[:, :, 2], mask, mode='same'), 0., 1.)
    
    # on retourne l'image modifiée
    sharpen_image = (255. * sharpen_image).astype(np.uint8)
    modified_image = Image.fromarray(sharpen_image, 'HSV').convert('RGB')
    # arrêter le chronomètre après la sauvegarde
    end_time = time.time()
    total_time = end_time - start_time
    print(f"total execution time: {total_time:.2f} seconds")
    return modified_image

path = "datas/"
image = path + "paysage2.jpg"
doubled_image = double_size(image)

# on sauvegarde l'image modifiée
doubled_image.save("sorties/paysage_double.jpg")
print("image saved") """

from PIL import Image
import os
import numpy as np
from scipy import signal
import time
from mpi4py import MPI

# fonction pour doubler la taille d'une image sans trop la pixeliser
def double_size(image_section):
    # on convertit l'image en tableau numpy
    img = np.array(image_section, dtype=np.double)
    # on double sa taille et on la normalise
    img = np.repeat(img, 2, axis=0)  # étendre verticalement
    img = np.repeat(img, 2, axis=1)  # étendre horizontalement
    img = img / 255.0
    # on crée un masque de flou gaussien
    mask = np.array([[1., 2., 1.], [2., 4., 2.], [1., 2., 1.]]) / 16.
    # on applique le filtre de flou
    blur_image = np.zeros_like(img, dtype=np.double)
    for i in range(3):
        blur_image[:, :, i] = signal.convolve2d(img[:, :, i], mask, mode='same')
    # on crée un masque de netteté
    mask = np.array([[0, -1, 0], [-1, 5, -1], [0, -1, 0]])
    # on applique le filtre de netteté uniquement sur la luminance
    sharpen_image = np.zeros_like(img, dtype=np.double)
    sharpen_image[:, :, :2] = blur_image[:, :, :2]
    sharpen_image[:, :, 2] = np.clip(signal.convolve2d(blur_image[:, :, 2], mask, mode='same'), 0., 1.)
    # on retourne l'image modifiée
    return (255. * sharpen_image).astype(np.uint8)

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

path = "datas/"
image_path = path + "paysage2.jpg"

if rank == 0:
    # on charge l'image complète sur le processus principal
    img = Image.open(image_path)
    print(f"original size: {img.size}")
    img = img.convert('HSV')
    img_array = np.array(img, dtype=np.uint8)
    height, width, channels = img_array.shape
    # diviser l'image en sections pour chaque processus avant l'agrandissement
    sections = np.array_split(img_array, size, axis=0)
else:
    sections = None

# distribuer les sections d'image aux processus
local_section = comm.scatter(sections, root=0)

# démarrer le chronomètre
start_time = time.time()

# doubler la taille de la section locale
processed_section = double_size(local_section)

# arrêter le chronomètre
comm.Barrier()
end_time = time.time()
total_time = end_time - start_time
print(f"process {rank}: execution time: {total_time:.2f} seconds")

# récupérer toutes les sections agrandies avec leurs indices pour maintenir l'ordre
gathered_sections = comm.gather((rank, processed_section), root=0)

if rank == 0:
    # trier les sections pour garantir l'ordre correct
    gathered_sections.sort()
    sorted_sections = [section for _, section in gathered_sections]
    
    # reconstruire l'image complète
    final_image = np.vstack(sorted_sections)
    final_image = Image.fromarray(final_image, 'HSV').convert('RGB')
    final_image.save("sorties/paysage_double.jpg")
    print("image saved")

# calculer le temps total d'exécution
max_time = comm.reduce(total_time, op=MPI.MAX, root=0)
if rank == 0:
    print(f"total execution time: {max_time:.2f} seconds")