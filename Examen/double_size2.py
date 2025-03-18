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
    print(f"taille originale {img.size}")
    img = img.convert('HSV')
    # on convertit l'image en tableau numpy
    img = np.array(img, dtype=np.double)
    # on double sa taille et on la normalise
    img = np.repeat(np.repeat(img, 2, axis=0), 2, axis=1) / 255.
    print(f"nouvelle taille : {img.shape}")
    # on crée un masque de flou gaussien pour la teinte et la saturation (H et S)
    mask = np.array([[1., 2., 1.], [2., 4., 2.], [1., 2., 1.]]) / 16.
    # on applique le filtre de flou
    blur_image = np.zeros_like(img, dtype=np.double)
    for i in range(2):
        blur_image[:, :, i] = signal.convolve2d(img[:, :, i], mask, mode='same')
    blur_image[:, :, 2] = img[:, :, 2]
    # on crée un masque de flou 5x5
    mask = -np.array([[1., 4., 6., 4., 1.], [4., 16., 24., 16., 4.], [6., 24., -476., 24., 6.], [4., 16., 24., 16., 4.], [1., 4., 6., 4., 1.]]) / 256
    # on applique le filtre sur la luminance
    blur_image[:, :, 2] = np.clip(signal.convolve2d(blur_image[:, :, 2], mask, mode='same'), 0., 1.)
    blur_image = (255. * blur_image).astype(np.uint8)
    # arrêter le chronomètre
    end_time = time.time()
    total_time = end_time - start_time
    print(f"temps total d'exécution : {total_time:.2f} secondes")
    # on retourne l'image modifiée
    return Image.fromarray(blur_image, 'HSV').convert('RGB')

path = "datas/"
image = path + "paysage2.jpg"
doubled_image = double_size(image)

# on sauvegarde l'image modifiée
doubled_image.save("sorties/paysage_double_2.jpg")
print("image sauvegardée") """

from mpi4py import MPI
from PIL import Image
import os
import numpy as np
from scipy import signal

def compute_local_indices(total_rows, size, rank):
    """
    Répartit équitablement les lignes de l'image entre les processus.
    Renvoie les indices valides (début, fin) pour la tranche du processus.
    """
    rows_per_proc = total_rows // size
    remainder = total_rows % size
    if rank < remainder:
        start = rank * (rows_per_proc + 1)
        end = start + rows_per_proc + 1
    else:
        start = remainder * (rows_per_proc + 1) + (rank - remainder) * rows_per_proc
        end = start + rows_per_proc
    return start, end

def apply_gaussian_filter_local(slice_array):
    """
    Applique le filtre gaussien FG (3x3) sur une tranche locale.
    """
    mask = np.array([[1., 2., 1.],[2., 4., 2.],[1., 2., 1.]], dtype=np.float32) / 16.
    L, W, C = slice_array.shape
    result = np.empty((L, W, C), dtype=np.float32)
    for ch in range(C):
        result[:, :, ch] = signal.convolve2d(slice_array[:, :, ch], mask, mode='same')
    return result

def apply_sharpen_filter_local(channel_array):
    """
    Applique le filtre de netteté FS (3x3) sur un canal 2D (la luminance V).
    """
    mask = np.array([[0, -1, 0],[-1, 5, -1],[0, -1, 0]], dtype=np.float32)
    conv_result = signal.convolve2d(channel_array, mask, mode='same')
    return np.clip(conv_result, 0., 1.)

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

if rank == 0:
    image_path = "datas/paysage2.jpg"
    img = Image.open(image_path)
    print(f"[Processus 0] Taille originale: {img.size}")
    img = img.convert('HSV')
    # Doubler la taille de l'image pour limiter la pixellisation
    img = np.repeat(np.repeat(np.array(img), 2, axis=0), 2, axis=1)
    # Utiliser np.float32 pour réduire l'utilisation de mémoire
    img = img.astype(np.float32) / 255.
    
    total_rows = img.shape[0]
    slices = []
    # Découper l'image en tranches horizontales pour chaque processus
    for r in range(size):
        valid_start, valid_end = compute_local_indices(total_rows, size, r)
        ghost_top = 1 if valid_start > 0 else 0
        ghost_bottom = 1 if valid_end < total_rows else 0
        # Extraire la tranche avec lignes fantômes
        slice_data = img[valid_start - ghost_top : valid_end + ghost_bottom, :, :].copy()
        valid_rows = valid_end - valid_start
        # On enregistre aussi l'indice de départ pour reconstituer l'image dans l'ordre
        slices.append((slice_data, ghost_top, ghost_bottom, valid_rows, valid_start))
else:
    slices = None

# Le processus 0 envoie les sous-images aux autres processus via send/recv
if rank == 0:
    # Envoyer aux processus 1 à size-1
    for r in range(1, size):
        comm.send(slices[r], dest=r, tag=77)
    local_data = slices[0]
else:
    local_data = comm.recv(source=0, tag=77)

(local_slice, ghost_top, ghost_bottom, valid_rows, valid_start) = local_data

if rank == 0:
    start_time = MPI.Wtime()

# Chaque processus applique le filtre gaussien sur sa tranche locale
local_gauss = apply_gaussian_filter_local(local_slice)
# Extraction de la partie valide (sans les lignes fantômes)
local_valid = local_gauss[ghost_top : ghost_top + valid_rows, :, :]

# Application du filtre de netteté sur le canal V (indice 2) de la partie valide
local_channel_V = local_valid[:, :, 2]
local_sharp = apply_sharpen_filter_local(local_channel_V)
local_valid[:, :, 2] = local_sharp

# Rassemblement des tranches traitées via gather : chaque tuple contient (valid_start, tranche valide)
gathered = comm.gather((valid_start, local_valid), root=0)

if rank == 0:
    # Trier les tranches par ordre croissant de valid_start
    gathered.sort(key=lambda x: x[0])
    processed_slices = [slice_valid for (start, slice_valid) in gathered]
    full_img = np.vstack(processed_slices)
    final_array = (full_img * 255).astype(np.uint8)
    final_img = Image.fromarray(final_array, 'HSV').convert('RGB')
    out_path = "sorties/"
    if not os.path.exists(out_path):
        os.makedirs(out_path)
    final_img.save(os.path.join(out_path, "paysage_double.jpg"))
    end_time = MPI.Wtime()
    total_time = end_time - start_time
    print("Image sauvegardée")
    print(f"Traitement terminé en {total_time:.4f} secondes avec {size} processus.")