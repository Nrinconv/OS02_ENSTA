# Ce programme va charger n images et y appliquer un filtre de netteté
# puis les sauvegarder dans un dossier de sortie

from PIL import Image
import os
import numpy as np
from scipy import signal
import time
from mpi4py import MPI


# Fonction pour appliquer un filtre de netteté à une image
def apply_filter(image):
    # On charge l'image
    img = Image.open(image)
    print(f"Taille originale {img.size}")
    # Conversion en HSV :
    img = img.convert('HSV')
    # On convertit l'image en tableau numpy et on normalise
    img = np.repeat(np.repeat(np.array(img), 2, axis=0), 2, axis=1)
    img = np.array(img, dtype=np.double)/255.
    print(f"Nouvelle taille : {img.shape}")
    # Tout d'abord, on crée un masque de flou gaussien
    mask = np.array([[1., 2., 1.], [2., 4., 2.], [1., 2., 1.]]) / 16.
    # On applique le filtre de flou
    blur_image = np.zeros_like(img, dtype=np.double)
    for i in range(3):
        blur_image[:,:,i] = signal.convolve2d(img[:,:,i], mask, mode='same')
    # On crée un masque de netteté
    mask = np.array([[0., -1., 0.], [-1., 5., -1.], [0., -1., 0.]])
    # On applique le filtre de netteté
    sharpen_image = np.zeros_like(img)
    sharpen_image[:,:,:2] = blur_image[:,:,:2]
    sharpen_image[:,:,2] = np.clip(signal.convolve2d(blur_image[:,:,2], mask, mode='same'), 0., 1.)

    sharpen_image *= 255.
    sharpen_image = sharpen_image.astype(np.uint8)
    # On retourne l'image modifiée
    return Image.fromarray(sharpen_image, 'HSV').convert('RGB')


""" path = "datas/perroquets/"
# On crée un dossier de sortie
if not os.path.exists("sorties/perroquets"):
    os.makedirs("sorties/perroquets")
out_path = "sorties/perroquets/"

output_images = []

# démarrer le chronomètre
start_time = time.time()

for i in range(37):
    image = path + "Perroquet{:04d}.jpg".format(i+1)
    sharpen_image = apply_filter(image)
    output_images.append(sharpen_image)
    print(f"image {i+1} processed")

# arrêter le chronomètre
end_time = time.time()
total_time = end_time - start_time
print(f"total execution time: {total_time:.2f} seconds")

# on sauvegarde les images modifiées
for i, img in enumerate(output_images):
    img.save(out_path + "Perroquet{:04d}.jpg".format(i+1))
print("images saved") """

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

path = "datas/perroquets/"
out_path = "sorties/perroquets/"
if rank == 0:
    # on crée un dossier de sortie
    if not os.path.exists(out_path):
        os.makedirs(out_path)
    image_list = [f"Perroquet{i+1:04d}.jpg" for i in range(37)]
    # diviser la liste des images pour chaque processus
    chunks = [image_list[i::size] for i in range(size)]
else:
    chunks = None

# répartir les chunks entre les processus
local_images = comm.scatter(chunks, root=0)

start_time = time.time()

# traiter les images de chque processus
processed_images = []
for image in local_images:
    sharpen_image = apply_filter(path + image)
    print(f"image {image} processed")
    processed_images.append((image, sharpen_image))

# s'assurer que tous les processus ont terminé pour calculer le temp total
comm.Barrier()
end_time = time.time()
total_time = end_time - start_time
max_time = comm.reduce(total_time, op=MPI.MAX, root=0)
if rank == 0:
    print(f"total execution time: {max_time:.2f} seconds")

# chaque processus envoie ses images modifiées au processus 0
gathered_images = comm.gather(processed_images, root=0)

if rank == 0:
    all_images = [img for sublist in gathered_images for img in sublist]
    all_images.sort()
    for img_name, img in all_images:
        img.save(out_path + img_name)
    print("images saved")