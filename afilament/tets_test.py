import cv2
import numpy as np
from scipy import ndimage as ndi
from skimage.segmentation import watershed
from skimage.feature import peak_local_max


def detect_circles(input_path, output_path):
    # Read the image in grayscale mode.
    img = cv2.imread(input_path, cv2.IMREAD_GRAYSCALE)

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

    for region in np.unique(labels):
        # skip background
        if region == 0:
            continue

        # Get the region properties
        region_mask = ((labels == region) * 255).astype(np.uint8)
        props = cv2.connectedComponentsWithStats(region_mask, connectivity=8)

        # Draw a circle for each connected component
        for i in range(1, props[0]):
            center = (int(props[3][i][0]), int(props[3][i][1]))
            radius = int((props[2][i][2] ** 2 + props[2][i][3] ** 2) ** 0.5 / 2)

            # Draw the circle on the original image.
            cv2.circle(img_display, center, radius, (0, 255, 255), 2)  # Yellow color

    return img_display


# Test the function.
detect_circles(r"C:\Users\nnina\Desktop\23-5-19 Kash 12J -LIV PFA-Triton 0hr-01.png", r"C:\Users\nnina\Desktop\out.png")
