import glob
import os
import math
import cv2
import numpy as np
from scipy import ndimage as ndi
from skimage.segmentation import watershed
from skimage.feature import peak_local_max
import pyvista as pv
import matplotlib.pyplot as plt
from objects import Contour

from unet.predict import run_predict_unet, run_predict_unet_one_img


def prepare_folder(folder):
    """
    Create folder if it has not been created before
    or clean the folder
    ---
    Parameters:
    -   folder (string): folder's path
    """
    if not os.path.exists(folder):
        os.makedirs(folder)
    for f in glob.glob(folder + "/*"):
        os.remove(f)


def find_max_projection(input_folder, identifier, show_img=False):
    object_layers = []
    for img_path in glob.glob(os.path.join(input_folder, "*_" + identifier + "_*.png")):
        layer = int(img_path.rsplit(".", 1)[0].rsplit("_", 1)[1])
        object_img = cv2.imread(img_path, cv2.IMREAD_ANYDEPTH)

        object_layers.append([object_img, layer])

    object_layers = sorted(object_layers, key=lambda x: x[1], reverse=True)
    image_3d = np.asarray([img for img, layer in object_layers])
    image_3d = cv2.normalize(image_3d, None, 0, 255, cv2.NORM_MINMAX, dtype=cv2.CV_8U) #This line was included to solve problem with weak nucleus signal cells
    max_projection_origin_size = image_3d[:, :, :].max(axis=0, out=None, keepdims=False, where=True)
    max_progection_unet_size = cv2.resize(max_projection_origin_size, (512, 512))
    if show_img:
        cv2.imshow("Max projection", max_progection_unet_size)
        cv2.waitKey()
    return max_projection_origin_size, max_progection_unet_size


def find_rotation_angle(input_folder, norm_th, rotation_trh=50):
    """
    Find rotation angle based on Hough lines of edges of maximum actin projection.
    ---
    Parameters:
        - input_folder (string): folder's path
        - rotation_trh (int): threshold for canny edges
    ___
    Returns:
        - rot_angle (int): rotaion angle in degrees
        - max_progection_img (img): image for verification
        - hough_lines_img (img): image for verification
    """
    identifier = "actin"
    max_projection, max_progection_img = find_max_projection(input_folder, identifier, norm_th)

    # Find the edges in the image using canny detector
    edges = cv2.Canny(max_projection, rotation_trh, 100)

    # Detect points that form a line
    lines = cv2.HoughLinesP(edges, 1, np.pi / 180, threshold=15, minLineLength=50, maxLineGap=250)

    # Draw lines on the image
    angles = []
    if lines is None:
        rot_angle = 0
        hough_lines_img = np.zeros((1000, 1000), dtype=np.uint8)

    else:
        for line in lines:
            x1, y1, x2, y2 = line[0]
            cv2.line(max_projection, (x1, y1), (x2, y2), (255, 0, 0), 3)
            if x1 == x2:
                angles.append(-90)
            else:
                angles.append(np.arctan((y2 - y1) / (x2 - x1)) * 180 / np.pi)
        # Create window with freedom of dimensions
        rot_angle = (np.median(angles))
        hough_lines_img = cv2.resize(max_projection, (1000, 1000))
        # cv2.imshow("output", cv2.resize(max_progection,(1000, 1000))) #keep it for debugging
        # cv2.waitKey()

    return -rot_angle, max_progection_img, hough_lines_img


def rotate_bound(image, angle):
    """
    Rotate provided image to specified angle
        ---
    Parameters:
        - image (img): image to rotate
        - angle (int): angle to rotate
    ___
    Returns:
        - rotated image
    """
    # grab the dimensions of the image and then determine the
    # center
    (h, w) = image.shape[:2]
    (cX, cY) = (w // 2, h // 2)
    # grab the rotation matrix (applying the negative of the
    # angle to rotate clockwise), then grab the sine and cosine
    # (i.e., the rotation components of the matrix)
    M = cv2.getRotationMatrix2D((cX, cY), -angle, 1.0)
    cos = np.abs(M[0, 0])
    sin = np.abs(M[0, 1])
    # compute the new bounding dimensions of the image
    nW = int((h * sin) + (w * cos))
    nH = int((h * cos) + (w * sin))
    # adjust the rotation matrix to take into account translation
    M[0, 2] += (nW / 2) - cX
    M[1, 2] += (nH / 2) - cY
    # perform the actual rotation and return the image
    return cv2.warpAffine(image, M, (nW, nH))


def get_3D(input_folder, identifier):
    """

    ---
    """
    object_layers = []
    for img_path in glob.glob(os.path.join(input_folder, "*_" + identifier + "_*.png")):
        layer = int(img_path.rsplit(".", 1)[0].rsplit("_", 1)[1])
        object_img = cv2.imread(img_path, cv2.IMREAD_ANYDEPTH)
        object_layers.append([object_img, layer])

    object_layers = sorted(object_layers, key=lambda x: x[1], reverse=True)
    image_3d = np.asarray([img for img, layer in object_layers])
    max_progection = image_3d[:, :, :].max(axis=0, out=None, keepdims=False, where=True)
    max_progection_img = cv2.resize(max_progection, (512, 512))
    image_3d = np.moveaxis(image_3d, 0, -1)  # Before moving axis (z, x, y), after moving axis (x, y, z)

    return image_3d, max_progection_img


def save_as_8bit(input_folder, output_folder, norm_th):
    """
    Convert images in the input folder to 8-bit using a common normalization threshold and save them in the output folder.

    Args:
        input_folder (str): Path to the input folder containing PNG images.
        output_folder (str): Path to the output folder where the 8-bit images will be saved.
        norm_th (int): Normalization threshold value.

    Returns:
        None.

    Raises:
        None.

    """

    for img_path in glob.glob(os.path.join(input_folder, "*.png")):
        object_img = cv2.imread(img_path, cv2.IMREAD_ANYDEPTH)

        #It is important to normalize using the same threshold for each image of the current image set, as using different
        # normalization thresholds might lead to inconsistent results. Therefore, the OpenCV normalization function with the
        # MinMax flag is not suitable in this case: img = cv2.normalize(img, None, alpha=0, beta=255, norm_type=cv2.NORM_INF, dtype=cv2.CV_8UC1)

        object_img = np.clip(object_img, 0, norm_th)
        img_8bit = np.uint8(object_img / norm_th * 255)

        output_img_path = os.path.join(output_folder, os.path.basename(img_path))
        cv2.imwrite(output_img_path, img_8bit)


def get_yz_xsection(img_3d, output_folder, identifier, cnt_extremes, unet_img_size=512):
    """
    Save jpg cross-section of img_3d with padding in output_folder.
    ---
        Parameters:
        - img_3d (np. array): the three-dimensional array that represents 3D image of the identifier (nucleus/actin)
        - output_folder (string): path to the folder to save processed jpg images
        - identifier (string) "actin" or "nucleus"
        - cnt_extremes (CntExtremes object) where left, right, top, bottom attributes are coordinates
                                    of the corresponding extreme points of the biggest nucleus contour
    ---
        Returns:
        - mid_cut (img): for verification
        - x_xsection_start (img): for verification
        - x_section_end (img): for verification
    """
    top, bottom = cnt_extremes.top[1], cnt_extremes.bottom[1]
    x_start, x_end, step = cnt_extremes.left[0], cnt_extremes.right[0], 1
    x_middle = x_start + (x_end - x_start) // 2
    mid_cut = None  # saved later for verification
    x_xsection_start, x_section_end = 0, 0
    for x_slice in range(x_start, x_end, step):
        xsection = img_3d[top: bottom, x_slice, :]
        img = xsection
        padded_img, x_xsection_start, x_section_end = make_padding(img, unet_img_size)
        img_path = os.path.join(output_folder, "xsection_" + identifier + "_" + str(x_slice) + ".png")
        cv2.imwrite(img_path, padded_img)
        if x_slice == x_middle:
            mid_cut = padded_img
    return mid_cut, x_xsection_start, x_section_end


def get_yz_xsection_3d(img_3d, cnt_extremes, unet_img_size=512):
    """
    Get yz cross-sections of img_3d with padding.
    ---
    Parameters:
    - img_3d (np. array): the three-dimensional array that represents a 3D image
    - cnt_extremes (CntExtremes object) where left, right, top, bottom attributes are coordinates
                                    of the corresponding extreme points of the biggest nucleus contour
    - unet_img_size (int): size for padding images
    ---
    Returns:
    - cross_sections (np. array): a 3D array containing the cross-sections
    """
    top, bottom = cnt_extremes.top[1], cnt_extremes.bottom[1]
    x_start, x_end, step = cnt_extremes.left[0], cnt_extremes.right[0], 1

    # Create a list to store the cross-sections
    cross_sections = []

    for x_slice in range(x_start, x_end, step):
        xsection = img_3d[top: bottom, x_slice, :]
        img = xsection
        padded_img, _, _ = make_padding(img, unet_img_size)

        # Append the cross-section to the list
        cross_sections.append(padded_img)

    # Convert the list of cross-sections to a 3D numpy array
    cross_sections = np.array(cross_sections)

    return cross_sections


def make_padding(img, final_img_size):

    h, w = img.shape[:2]

    if max([h, w]) > final_img_size:
        final_img_size = max([h, w])

    h_out = w_out = final_img_size

    top = (h_out - h) // 2
    bottom = h_out - h - top
    left = (w_out - w) // 2
    right = w_out - w - left

    padded_img = cv2.copyMakeBorder(img, top, bottom, left, right, cv2.BORDER_CONSTANT, value=0)
    return padded_img, left, (padded_img.shape[0] - right)


def get_3d_img(input_folder):
    """
    Reads and combines images from input_folder into image_3d according to layer number.
    ---
        Parameters:
        - input_folder (string): path to the input folder with jpg images
    ---
        Returns:
        - image_3d (np.array): three-dimensional array of images combined layer by layer together
    """
    object_layers = []
    i = 0

    for img_path in glob.glob(os.path.join(input_folder, "*.png")):
    #for img_path in glob.glob(input_folder + r"\*"): #FOR HISTORY this line caused issues on BSU cluster
        img_name, img_ext = os.path.splitext(os.path.basename(img_path))
        layer = int(img_name.rsplit("_", 1)[1])  # layer number is part of the image name

        img = cv2.imread(img_path, cv2.IMREAD_ANYDEPTH)
        print(f"Debugging image reading:\n"
              f"img # {i} size: {img.shape}\n"
              f"img path: {img_path}\n")
        i+=1
        object_layers.append([img, layer])

    object_layers = sorted(object_layers, key=lambda x: x[1], reverse=True)
    img_3d = np.asarray([mask for mask, layer in object_layers])

    return img_3d


def save_rotation_verification(cell, max_progection_img, hough_lines_img, rotated_max_projection,
                               mid_cut_img, part, folders, unet_parm):
    """
    Save images in specify folder
    ---
    """
    max_projection_init_file_path = os.path.join(folders["rotatation_verification"],
                                                 "img_num_" + str(cell.img_number) + "__cell_num_" + str(cell.number)
                                                 + "_" + part + "_max_projection_initial.png")
    max_projection_rotated_file_path = os.path.join(folders["rotatation_verification"],
                                                    "img_num_" + str(cell.img_number) + "__cell_num_" + str(cell.number)
                                                    + "_" + part + "_max_projection_rotated.png")
    hough_lines_file_path = os.path.join(folders["rotatation_verification"],
                                         "img_num_" + str(cell.img_number) + "__cell_num_" + str(cell.number)
                                         + "_" + part + "_hough_lines.png")
    cv2.imwrite(max_projection_init_file_path, max_progection_img)
    cv2.imwrite(max_projection_rotated_file_path, rotated_max_projection)
    cv2.imwrite(hough_lines_file_path, hough_lines_img)

    # save cross section for verification
    pixel_value = 65536 if mid_cut_img.dtype == "uint16" else 256
    mid_cut_img_8bit = np.uint8(mid_cut_img / (pixel_value / 256))
    middle_xsection_file_path = os.path.join(folders["middle_xsection"],
                                             "img_num_" + str(cell.img_number) + "__cell_num_" + str(cell.number)
                                             + "_" + part + "_middle_xsection.png")
    cv2.imwrite(middle_xsection_file_path, mid_cut_img_8bit)

    #save cross section unet output for verification

    mask_xsection_file_path = os.path.join(folders["middle_xsection"],
                                                 "img_num_" + str(cell.img_number) + "__cell_num_" + str(cell.number)
                                                 + "_" + part + "_mask_xsection.png")

    mid_cut_mask = run_predict_unet_one_img(middle_xsection_file_path, unet_parm.actin_unet_model, unet_parm.unet_model_scale,
                                            unet_parm.unet_model_thrh).astype(np.uint8)
    mid_cut_mask *= 255

    cv2.imwrite(mask_xsection_file_path, mid_cut_mask)



def find_biggest_nucleus_layer(temp_folders, treshold, find_biggest_mode, unet_parm=None):
    """
    Finds and analyzes image (layer) with the biggest area of the nucleus
    ---
        Parameters:
        - input_folder (string): path to the folder where all slices of the nucleus
                            in jpg format is located
    ---
        Returns:
        - biggest_nucleus_mask (np. array): array of 0 and 1 where 1 is white pixels
                                        which represent the shape of the biggest area
                                        of nucleus over all layers and 0 is a background
    """
    nucleus_area = 0
    mask, center = None, None

    if find_biggest_mode == "unet":
        run_predict_unet(temp_folders["raw"], temp_folders["nucleus_top_mask"], unet_parm.from_top_nucleus_unet_model,
                         unet_parm.unet_model_scale,
                         unet_parm.unet_model_thrh)
        folder = temp_folders["nucleus_top_mask"]
    elif find_biggest_mode == "trh":
        folder = temp_folders["raw"]

    for img_path in glob.glob(os.path.join(folder, "*_nucleus_*.png")):
        nucleus_img = cv2.imread(img_path, cv2.IMREAD_GRAYSCALE)

        current_nucleus_cnt = Contour.get_biggest_cnt(nucleus_img, treshold)
        if current_nucleus_cnt is None:
            continue
        current_nucleus_cnt_area = cv2.contourArea(current_nucleus_cnt)
        if current_nucleus_cnt_area > nucleus_area:
            nucleus_area = current_nucleus_cnt_area
            mask = Contour.draw_cnt(current_nucleus_cnt, nucleus_img.shape[:2])
    return mask


def сut_out_mask(mask, input_folder, output_folder, identifier):
    """
    Cuts out an area that corresponds to the mask on each image (layer) located in the input_folder,
    saves processed images in the output_folder, and returns processed images combined into image_3d
    ---
        Parameters:
        - input_folder (string): path to the input folder with jpg images
        - output_folder (string): path to the folder to save processed jpg images
        - identifier (string): "actin" or "nucleus"
        - mask (np. array): stencil to cut out from the images
    ---
        Returns:
        - image_3d (np. array): three-dimensional array of processed (cut out) images combined layer by layer together

    """
    object_layers = []
    for img_path in glob.glob(os.path.join(input_folder, "*_" + identifier + "_*.png")):
        layer = int(img_path.rsplit(".", 1)[0].rsplit("_", 1)[1])
        object_img = cv2.imread(img_path, cv2.IMREAD_ANYDEPTH)
        mask[mask == 255] = 1
        object_layer = np.multiply(object_img, mask)
        object_layers.append([object_layer, layer])
        cv2.imwrite(os.path.join(output_folder, os.path.basename(img_path)), object_layer)

    object_layers = sorted(object_layers, key=lambda x: x[1], reverse=True)
    image_3d = np.asarray([img for img, layer in object_layers], dtype=np.uint8)
    image_3d = np.moveaxis(image_3d, 0, -1)
    return image_3d


def combine_masked_images_into_3d(mask, input_folder, identifier):
    """
    Cuts out an area that corresponds to the mask on each image (layer) located in the input_folder,
    and returns processed images combined into image_3d
    ---
        Parameters:
        - input_folder (string): path to the input folder with jpg images
        - identifier (string): channel name that is part of the image file names
        - mask (np. array): stencil to cut out from the images
    ---
        Returns:
        - image_3d (np. array): three-dimensional array of processed (cut out) images combined layer by layer together
    """
    object_layers = []
    for img_path in glob.glob(os.path.join(input_folder, "*.png")):
        # Check if the identifier is anywhere in the image file name
        if identifier in os.path.basename(img_path):
            layer = int(img_path.rsplit(".", 1)[0].rsplit("_", 1)[1])
            object_img = cv2.imread(img_path, cv2.IMREAD_ANYDEPTH)
            mask[mask == 255] = 1
            object_layer = np.multiply(object_img, mask)
            object_layers.append([object_layer, layer])

    object_layers = sorted(object_layers, key=lambda x: x[1], reverse=True)
    image_3d = np.asarray([img for img, layer in object_layers])
    image_3d = np.moveaxis(image_3d, 0, -1)
    return image_3d


def plot_histogram(title, image):
    histogram, bin_edges = np.histogram(image, bins=256 * 256)
    plt.figure()
    plt.title(title)
    plt.xlabel("grayscale value")
    plt.ylabel("pixel count")

    plt.plot(histogram)  # <- or here
    plt.show()


def is_point_in_pyramid_old_version(x, y, z, x_test, y_test, z_test, con_angle, min_len, resolution):
    is_x_in_pyramid = x_test < x + min_len
    # Find y frame based on known angle and the height of the triangle (con_angle)
    sin = math.sin(math.radians(con_angle))
    y_frame = abs(sin * (x - x_test))
    is_y_in_pyramid = y - y_frame <= y_test <= y + y_frame
    z_frame = int(y_frame * (
                resolution.y / resolution.z))  # Since pixel size is different for xy and xz planes, adjustment should be done
    is_z_in_pyramid = z - z_frame <= z_test <= z + z_frame
    return is_x_in_pyramid and is_y_in_pyramid and is_z_in_pyramid


def is_point_in_pyramid(right_rotated_xy, right_rotated_xz, left_rotated_xy, left_rotated_xz, con_angle, con_angle_z,
                        min_len, resolution):
    x_y, y = right_rotated_xy
    x_y_test, y_test = left_rotated_xy
    x_z, z = right_rotated_xz
    x_z_test, z_test = left_rotated_xz
    is_x_y_in_pyramid = x_y_test < x_y + min_len
    # Find y frame based on known angle and the height of the triangle (con_angle)
    sin = math.sin(math.radians(con_angle))
    sin_z = math.sin(math.radians(con_angle_z))
    y_frame = abs(sin * (x_y - x_y_test))
    is_y_in_pyramid = y - y_frame <= y_test <= y + y_frame
    z_frame = int(abs(sin_z * (x_z - x_z_test)) * (
                resolution.y / resolution.z))  # Since pixel size is different for xy and xz planes, adjustment should be done
    is_x_z_in_pyramid = x_z_test < x_z + min_len
    is_z_in_pyramid = z - z_frame <= z_test <= z + z_frame
    return is_x_y_in_pyramid and is_y_in_pyramid and is_x_z_in_pyramid and is_z_in_pyramid


def rotate_point(rotated_point, main_point, rot_angle):
    """
    Based on https://math.stackexchange.com/questions/3244392/getting-transformed-rotated-coordinates-relative-to
    -the-center-of-a-rotated-ima example the new coordinates of the candidates' starting point (p1,p2) after a rotation by θ
    degrees around the ending point of the main line (a,b): (a + x cos θ − y sin θ,   b + x sin θ + y cos θ)
    where x = p1−a, y = p2−b
    """
    x = rotated_point[0] - main_point[0]
    y = rotated_point[1] - main_point[1]
    new_coordinates = (
    int(main_point[0] + x * math.cos(math.radians(rot_angle)) - y * math.sin(math.radians(rot_angle))),
    int(main_point[1] + x * math.sin(math.radians(rot_angle)) + y * math.cos(math.radians(rot_angle))))
    return new_coordinates


def get_nuclei_masks(temp_folders, output_analysis_folder, image_path, nuc_theshold,
                     nuc_area_min_pixels_num, find_biggest_mode, img_num, unet_parm=None):
    img_base_path = os.path.splitext(os.path.basename(image_path))[0]
    max_projection_origin_size, max_progection_unet_size = find_max_projection(temp_folders["nuc_raw"], "nucleus", show_img=False)
    max_projection_path = os.path.join(temp_folders["nucleus_top_img"], img_base_path + ".png")
    cv2.imwrite(max_projection_path, max_projection_origin_size)
    dim = 0
    cnts = []
    if find_biggest_mode == "unet":
        run_predict_unet(temp_folders["nucleus_top_img"], temp_folders["nucleus_top_mask"],
                         unet_parm.from_top_nucleus_unet_model,
                         unet_parm.unet_model_scale,
                         unet_parm.unet_model_thrh)
        mask_path = os.path.join(temp_folders["nucleus_top_mask"], img_base_path + ".png")
        mask_img = cv2.imread(mask_path, cv2.IMREAD_GRAYSCALE)
        cnts = Contour.get_img_cnts(mask_img, nuc_theshold)
        dim = mask_img.shape

    elif find_biggest_mode == "trh":
        img_path = os.path.join(temp_folders["nucleus_top_img"], img_base_path + ".png")
        nucleus_img = cv2.imread(img_path, cv2.IMREAD_GRAYSCALE)
        cnts = Contour.get_img_cnts(nucleus_img, nuc_theshold)
        dim = nucleus_img.shape

    # Remove nuclei that touch edges of the image
    # cnts = remove_edge_nuc(cnts, dim) #I commented it out since for current analysis it would be beneficial remove it manually
    cnts = [cnt for cnt in cnts if cv2.contourArea(cnt) > nuc_area_min_pixels_num]  # removes noise

    nuclei_masks = []

    for i, cnt in enumerate(cnts):
        one_nuc_mask = np.zeros(dim, dtype="uint8")
        cv2.drawContours(one_nuc_mask, [cnt], -1, color=(255, 255, 255), thickness=cv2.FILLED)
        kernel = np.ones((10, 10), np.uint8)
        one_nuc_mask_dilation = cv2.dilate(one_nuc_mask, kernel, iterations=1)

        one_nuc_mask_path = os.path.join(temp_folders["nuclei_top_masks"], img_base_path + "_" + str(i) + ".png")
        cv2.imwrite(one_nuc_mask_path, one_nuc_mask_dilation)
        nuclei_masks.append(one_nuc_mask_dilation)

    draw_and_save_cnts_verification(output_analysis_folder, image_path, cnts, max_projection_origin_size, img_num)
    nuclei_xy_centers = [Contour.get_cnt_center(cont) for cont in cnts]

    return nuclei_masks, nuclei_xy_centers


def draw_and_save_cnts_verification(output_analysis_folder, image_path, cnts, max_progection_img, img_num):
    base_img_name = os.path.splitext(os.path.basename(image_path))[0]
    ver_img_path = os.path.join(output_analysis_folder,
                                base_img_name + "_img-num_" + str(img_num) + "_max_projection.png")

    cv2.drawContours(max_progection_img, cnts, -1, (255, 255, 255), 5)

    for i, cnt in enumerate(cnts):
        org = Contour.get_cnt_center(cnt)
        cv2.putText(max_progection_img, str(i), org, fontFace=cv2.FONT_HERSHEY_PLAIN, fontScale=3, color=(255, 255, 0),
                    thickness=3)
    cv2.imwrite(ver_img_path, max_progection_img)


def remove_edge_nuc(cnts, img_dim):  # removes cells that touch the edges of the frame
    max_y, max_x = img_dim

    # Countours should not touch the edges of the image
    new_cnts = []
    for cnt in cnts:
        cnt_max_x = cnt[:, 0, 0].max()
        cnt_max_y = cnt[:, 0, 1].max()
        cnt_min_x = cnt[:, 0, 0].min()
        cnt_min_y = cnt[:, 0, 1].min()
        if (cnt_max_x < max_x - 2) and (cnt_max_y < max_y - 2) and (cnt_min_x > 1) and (cnt_min_y > 1):
            new_cnts.append(cnt)

    return new_cnts


def detect_circles(img):

    # Create a copy of the original image for displaying the result
    img_display = img.copy()

    # Apply a Gaussian blur to reduce noise.
    blurred_img = cv2.GaussianBlur(img, (9, 9), 2)

    # Apply simple thresholding.
    _, thresh = cv2.threshold(blurred_img, 100, 255, cv2.THRESH_BINARY)

    # Calculate the distance from the thresholded image to the nearest zero pixel
    distance = ndi.distance_transform_edt(thresh)

    # Find peaks in the distance map as markers for the foreground
    coords = peak_local_max(distance, min_distance=5, labels=thresh)
    mask = np.zeros(distance.shape, dtype=bool)
    mask[tuple(coords.T)] = True
    markers, _ = ndi.label(mask)

    # Use watershed to segment the image
    labels = watershed(-distance, markers, mask=thresh)
    count = 0
    for region in np.unique(labels):
        # skip background
        if region == 0:
            continue

        # Get the region properties
        region_mask = ((labels == region) * 255).astype(np.uint8)
        props = cv2.connectedComponentsWithStats(region_mask, connectivity=8)
        count = count + 1
        # Draw a circle for each connected component
        for i in range(1, props[0]):
            center = (int(props[3][i][0]), int(props[3][i][1]))
            radius = int((props[2][i][2] ** 2 + props[2][i][3] ** 2) ** 0.5 / 2)

            # Draw the circle on the original image.
            cv2.circle(img_display, center, radius, (0, 255, 255), 2)  # Yellow color

    return img_display, count


def save_nuclei_meshes(meshes, nuclei_xy_centers, xy_size, folder_name, file_name):
    """
    Saves all meshes of nuclei as a single STL file using pyvista, positioned according to their centers,
    with a black border of size xy_size.

    :param meshes: List of pyvista mesh objects to be combined and saved as a single STL file.
    :param nuclei_xy_centers: List of (x, y) centers for each mesh in pixel coordinates.
    :param xy_size: Tuple (x, y) size in pixels of the surface border.
    :param folder_name: Name of the folder to save the combined mesh.
    :param file_name: Name for the output STL file.
    """
    # Create the directory if it doesn't exist
    mesh_folder = os.path.join('mesh', folder_name)
    if not os.path.exists(mesh_folder):
        os.makedirs(mesh_folder)

    combined_mesh = None

    # Create a black border as a flat mesh
    border = pv.Plane(center=(xy_size[0] / 2, xy_size[1] / 2, 0), i_size=xy_size[0], j_size=xy_size[1])
    border.color = 'black'

    for mesh_obj, center in zip(meshes, nuclei_xy_centers):
        # Calculate translation vector
        current_center = np.mean(mesh_obj.points, axis=0)
        translated_center = (center[0], center[1], 0)  # Assuming z-coordinate is 0
        translation_vector = np.array(translated_center) - current_center[:3]

        # Apply translation
        mesh_obj.translate(translation_vector, inplace=True)

        # Combine the mesh
        if combined_mesh is None:
            combined_mesh = mesh_obj
        else:
            combined_mesh = combined_mesh.merge(mesh_obj)

    # Merge the border with the combined mesh
    combined_mesh_with_border = combined_mesh.merge(border)

    # Save the combined mesh with border as a single STL file
    combined_mesh_with_border.save(os.path.join(mesh_folder, f"{file_name}.stl"))



