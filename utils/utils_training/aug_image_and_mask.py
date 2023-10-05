import os
import cv2
import numpy as np
from random import randrange
import random
from PIL import Image, ImageOps
import glob


if __name__ == '__main__':
    folder_path = r"D:\BioLab\img\training_sets\microtubules_comets\img"
    mask_folder_path =r"D:\BioLab\img\training_sets\microtubules_comets\mask"

    for i, img_path in enumerate(glob.glob(os.path.join(folder_path, "*.png"))):
        mask_path = glob.glob(os.path.join(mask_folder_path, "*.bmp"))[i]

        img_file = os.path.basename(img_path)
        mask_file = os.path.basename(mask_path)
        nucleus_img = Image.open(img_path)
        nucleus_mask = Image.open(mask_path)
        width, height = nucleus_img.size
        width_shift = 1 / 2 * (width - ((height / 2) ** 2 + (width / 2) ** 2) ** (1 / 2))
        height_shift = 1 / 2 * (height - ((height / 2) ** 2 + (width / 2) ** 2) ** (1 / 2))
        gamma = [0.80, 0.90, 1, 1, 1, 1.10, 1.20]

        for j in np.arange(3):
            #angle_to_rotate = randrange(360)
            gamma_corrected = Image.fromarray(np.array(255 * (np.array(nucleus_img)/ 255) ** random.choice(gamma), dtype='uint8'))
            angle_to_rotate = random.choice([0, 90, 180, 270])
            aug_img = gamma_corrected.rotate(angle_to_rotate)
            aug_mask = nucleus_mask.rotate(angle_to_rotate)

            #rotated_img = nucleus_img.rotate(angle_to_rotate)
            #cropped_img = rotated_img.crop((height_shift, width_shift, height - height_shift, width - width_shift))
            aug_mask = nucleus_mask.rotate(angle_to_rotate)
            #cropped_mask = rotated_mask.crop((height_shift, width_shift, height - height_shift, width - width_shift))
            if (random.choice([1, 0])):
                aug_img = ImageOps.flip(aug_img)
                aug_mask = ImageOps.flip(aug_mask)
            if (random.choice([1, 0])):
                aug_img = ImageOps.mirror(aug_img)
                aug_mask = ImageOps.mirror(aug_mask)

            img_file_name = '_aug_' + str(j) + '_' + img_file
            mask_file_name = '_aug_' + str(j) + '_' + mask_file
            img_path_to_save = os.path.join(r"D:\BioLab\img\training_sets\microtubules_comets\img_aug", img_file_name)
            mask_path_to_save = os.path.join(r"D:\BioLab\img\training_sets\microtubules_comets\mask_aug", mask_file_name)
            aug_img.save(img_path_to_save)
            aug_mask.save(mask_path_to_save)
            print("Processing " + img_path)

        # for k in np.arange(10):
        #     random_left_shift = randrange((int)(width_shift * 2))
        #     random_up_shift = randrange((int)(height_shift * 2))
        #     right_shift = width - (2 * width_shift - random_left_shift)
        #     bottom_shift = height - (2 * height_shift - random_up_shift)
        #     cropped_shifted_img = nucleus_img.crop((random_up_shift, random_left_shift, bottom_shift, right_shift))
        #     cropped_shifted_mask = nucleus_mask.crop((random_up_shift, random_left_shift, bottom_shift, right_shift))
        #     img_file_name = '_translated_' + str(k) + '_' + img_file
        #     mask_file_name = '_translated_' + str(k) + '_' + mask_file
        #     img_path_to_save = os.path.join(r"C:\Users\nnina\Desktop\Aug_all_labled_nuclei_IMG", img_file_name)
        #     mask_path_to_save = os.path.join(r"C:\Users\nnina\Desktop\Aug_all_labled_nuclei_MASK", mask_file_name)
        #     cropped_shifted_img.save(img_path_to_save)
        #     cropped_shifted_mask.save(mask_path_to_save)
 








