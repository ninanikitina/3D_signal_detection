import glob
import cv2.cv2 as cv2
import os

if __name__ == '__main__':
    input_folder = r"C:\Users\nnina\Desktop\mask_max_pr"
    output_folder = r"C:\Users\nnina\Desktop\512_mask_max_pr"

    for img_path in glob.glob(os.path.join(input_folder, "*.bmp")):
        img = cv2.imread(img_path, cv2.IMREAD_GRAYSCALE)
        name_1 = os.path.split(img_path)[1].split("_")[0]
        name_2 = os.path.split(img_path)[1].split(" ")[1].split("_")[-1]
        resized_img = cv2.resize(img.astype('uint8'), (512, 512), interpolation=cv2.INTER_AREA)

        output_img_path = os.path.join(output_folder, name_1 + "_" + name_2 + ".bmp")
        cv2.imwrite(output_img_path, resized_img)
        print(output_img_path)