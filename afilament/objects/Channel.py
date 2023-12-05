import numpy as np
import cv2
import os
import random
import string


class Channel:
    def __init__(self, name: str):
        self.name = name
        self.av_signal_in_nuc_area_3D = None
        self.cut_3d_img = None
        self.sum_pix_in_nuc_cylinder = None
        self.has_ring = None
        self.ring_intensity_coef = None

    def convert_to_8bit(self, image):
        # Determine the minimum and maximum pixel values in the image
        min_val = np.min(image)
        max_val = np.max(image)

        # Scale the image to the 8-bit range
        scaled_image = ((image - min_val) / (max_val - min_val)) * 255

        # Convert the scaled image to the 8-bit unsigned integer format
        converted_image = np.uint8(scaled_image)

        return converted_image


    def detect_ring_on_nucleus(self, biggest_slice_index, cut_img_3d, nuc_3D_whole_img):
        # For 40x oil lens, average nucleus is 8 slices, so we can check 5 slices: 2 above and 2 below
        # the biggest slice if there is a ring on any of these slices
        self.ring_intensity_coef = -1
        print(f"Channel {self.name}\n")

        for i in range(-2, 3):

            current_slice = cut_img_3d[:, :, biggest_slice_index + i]
            nuc_slice = nuc_3D_whole_img[:, :, biggest_slice_index + i]
            nuc_slice_mask = (nuc_slice != 0).astype(np.uint8)

            # Convert the image to 8-bit unsigned integer format
            gray = self.convert_to_8bit(current_slice)


            # Erode the mask to create a ring (using a 20x20 kernel)
            kernel_size = 40
            eroded_mask = cv2.erode(nuc_slice_mask.astype(np.uint8), np.ones((kernel_size, kernel_size), np.uint8),
                                    iterations=1)

            # Subtract eroded mask from the original mask to get the ring area
            ring_mask = nuc_slice_mask - eroded_mask

            # Check if there are non-zero pixels in the ring mask
            if np.count_nonzero(ring_mask) == 0:
                continue  # Skip this iteration if there are no non-zero pixels in the ring mask

            # Create a mask where pixels greater than background_signal are set to 1, and others to 0
            # This method probably will help in identifying the ring
            # Calculate the threshold for the dimmest 5% of pixels
            background_signal = np.percentile(current_slice[current_slice > 0], 20)
            threshold_mask = (current_slice > background_signal).astype(np.uint8)

            # Calculate average signal in the ring area
            a = ring_mask * threshold_mask * gray
            av_signal_in_ring_area = np.sum(ring_mask * threshold_mask * gray)/np.sum(ring_mask*threshold_mask)

            # Calculate average signal in the non-ring area of the nucleus
            b = eroded_mask * threshold_mask * gray
            av_signal_in_nuc_non_ring_area = np.sum(eroded_mask * threshold_mask * gray)/np.sum(eroded_mask*threshold_mask)

            # Determine if the slice has a ring based on average signals
            current_ring_intensity_coef = av_signal_in_ring_area / av_signal_in_nuc_non_ring_area
            print(f"Ring intensity coefficient for slice {biggest_slice_index + i}: {current_ring_intensity_coef}\n")
            if current_ring_intensity_coef > self.ring_intensity_coef:
                self.ring_intensity_coef = current_ring_intensity_coef

        self.has_ring = self.ring_intensity_coef > 1


    def quantify_signals(self, nucleus, cut_img_3d, biggest_slice_index):
        """
        Quantifies signals in the 3D nucleus area, 2D area and ring area of the biggest slice.

        Parameters:
            nucleus (object): Nucleus object containing the 3D mask of the nucleus.
            nuc_max_projection_mask (ndarray): 2D binary mask of the biggest slice.
            cut_img_3d (ndarray): 3D image to be analyzed.
        """
        # Validate inputs
        if not hasattr(nucleus, 'nuc_3D_mask'):
            raise ValueError("Input 'nucleus' must have 'nuc_3D_mask' attribute")

        if not isinstance(cut_img_3d, np.ndarray) or cut_img_3d.ndim != 3:
            raise ValueError("'cut_img_3d' must be a 3-dimensional numpy array")

        # Calculate total signal in area that is corresponded to the biggest nucleus slice

        # Calculate total signal in nucleus area (3D)
        roi_3d = np.multiply(self.cut_3d_img, nucleus.nuc_3D_mask)
        self.sum_pix_in_nuc_cylinder = np.sum(cut_img_3d, dtype=np.uint64)
        total_signal_in_nucleus_3d = np.sum(roi_3d, dtype=np.int64)
        self.av_signal_in_nuc_area_3D = total_signal_in_nucleus_3d / np.count_nonzero(nucleus.nuc_3D_mask)

        self.detect_ring_on_nucleus(biggest_slice_index, cut_img_3d, nucleus.nuc_3D_whole_img)


    def quantify_intersection_with_nuc_channel(self, channel_biggest_slice, nuc_biggest_slice, intensity_persent=0.7):


        # Calculate average intensity
        nuc_slice_percentile_intensity = self.get_percentile_intensity(nuc_biggest_slice, intensity_persent)
        channel_slice_percentile_intensity = self.get_percentile_intensity(channel_biggest_slice, intensity_persent)

        # Create a mask of pixels with intensity greater than the threshold
        mask1 = cv2.inRange(nuc_biggest_slice, nuc_slice_percentile_intensity + 1, 255)
        mask2 = cv2.inRange(channel_biggest_slice, channel_slice_percentile_intensity + 1, 255)

        # Create empty color layers for each mask
        color_mask1 = np.zeros_like(mask1)
        color_mask2 = np.zeros_like(mask1)

        # Make the nucleus mask yellow (BGR format: [0,255,255])
        color_mask1[mask1 > 0] = [0, 255, 255]

        # Make the channel mask green (BGR format: [0,255,0])
        color_mask2[mask2 > 0] = [0, 255, 0]

        # Combine the masks
        combined = cv2.addWeighted(color_mask1, 1, color_mask2, 1, 0)

        # Create a mask for the intersection
        intersection = np.logical_and(mask1, mask2)

        # Make the intersection red (BGR format: [0,0,255])
        combined[intersection] = [0, 0, 255]


        # Calculate the percentage of intersection
        total_pixels = np.size(mask1)
        intersection_pixels = np.sum(intersection)
        intersection_percent = (intersection_pixels / total_pixels) * 100

        print(f"Percentage of intersection: {intersection_percent}%")
        print("Combined image saved as 'combined.jpg'")

    def get_percentile_intensity(self, img, intensity_persent):
        # Find all non-zero pixels
        non_zero_pixels = img[np.nonzero(img)]

        # Sort the pixel intensities in descending order
        sorted_pixels = np.sort(non_zero_pixels)[::-1]

        # Find the number of pixels in the top x%
        num_top_x_percent = int(len(sorted_pixels) * intensity_persent)

        # Get the intensity threshold for the top x%
        threshold_intensity = sorted_pixels[num_top_x_percent]

        return threshold_intensity




