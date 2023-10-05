import os
import cv2.cv2 as cv2
from pathlib import Path


input_folder = r"D:\BioLab\img\training_sets\Nucleus_training_img_and_masks\from_the_top_max_projection\top_nuc-4additional_40x_zeiss\mask_bmp_color_converted"
imges_path = os.listdir(input_folder)
output_folder = r"D:\BioLab\img\training_sets\Nucleus_training_img_and_masks\from_the_top_max_projection\top_nuc-4additional_40x_zeiss\mask"
for img_path in imges_path:
    img = cv2.imread(input_folder + "\\" + img_path, cv2.IMREAD_UNCHANGED)
    img_name = Path(img_path).stem
    b = os.path.join(output_folder, img_name, ".png")
    cv2.imwrite(os.path.join(output_folder, img_name + ".png"), img)