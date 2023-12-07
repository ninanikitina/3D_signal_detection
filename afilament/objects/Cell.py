import cv2
import numpy as np
from objects.Nucleus import Nucleus
from objects.Channel import Channel
from objects import Utils


class Cell(object):
    def __init__(self, img_num, cell_num, channel_names):
        self.img_number = img_num
        self.number = cell_num
        self.nucleus = None
        self.foci_count = None
        self.channels = [Channel(channel_name) for channel_name in channel_names]

    def count_foci(self, foci_image):
        nuc_3D_normalized = cv2.normalize(self.nucleus.nuc_3D_whole_img, None, 0, 255, cv2.NORM_MINMAX, dtype=cv2.CV_8U)
        mip_8bit = np.amax(nuc_3D_normalized, axis=2)
        img, foci_count = Utils.detect_circles(mip_8bit)

        # Add the mip_8bit image to the foci_image.
        foci_image = cv2.add(foci_image, img)

        self.foci_count = foci_count

        return foci_image  # Returns the accumulated foci_image.


    def analyze_nucleus(self, folders, unet_parm, img_resolution, analysis_folder, norm_th, cnt_extremes, nuc_max_projection_mask):
        """
        Run nucleus analysis of the cell
        """
        nucleus = Nucleus(nuc_max_projection_mask)
        nucleus.reconstruct(folders, unet_parm, img_resolution, analysis_folder, norm_th, cnt_extremes)
        self.nucleus = nucleus

    def analyze_channels(self, raw_img_folder, cnt_extremes, nuc_max_projection_mask, nuc_theshold, nucleus_channel_name):
        cut_nuc_img_3d_full_size = Utils.combine_masked_images_into_3d(nuc_max_projection_mask, raw_img_folder,
                                                                   nucleus_channel_name)
        biggest_slice_index = self.find_biggest_slice(cut_nuc_img_3d_full_size)
        for channel in self.channels:
            cut_img_3d_full_size = Utils.combine_masked_images_into_3d(nuc_max_projection_mask, raw_img_folder, channel.name)
            channel.cut_3d_img = Utils.get_yz_xsection_3d(cut_img_3d_full_size, cnt_extremes)
            channel.quantify_signals(self.nucleus, cut_img_3d_full_size, biggest_slice_index)

            # #This method of finding correspondence between two channels is based on masks intersection as an outcome this
            # channel_biggest_slice = cut_img_3d_full_size[:, :, biggest_slice_index]
            # nuc_biggest_slice = cut_nuc_img_3d_full_size[:, :, biggest_slice_index]
            # channel.quantify_intersection_with_nuc_channel(channel_biggest_slice, nuc_biggest_slice)

            print(f"{channel.name}: has ring {channel.has_ring} \n")


    def find_biggest_slice(self, cut_img_3d):
        # Find the index of the biggest slice (one with the highest intensity)
        slice_intensities = np.sum(cut_img_3d, axis=(0, 1))
        biggest_slice_index = np.argmax(slice_intensities)

        return biggest_slice_index

    def get_aggregated_cell_stat(self):
        """
        [_, "Img_num", "Cell_num", "Nucleus_volume, cubic_micrometre",  "Nucleus_cylinder, pixels_number",
        "Nucleus_length, micrometre", "Nucleus_width, micrometre", "Nucleus_high, micrometre", "Nucleus_total_intensity",
        "Channel_signal_in_nuc_area", "Channel_signal_in_ring_area", "Does_this_channel_has_ring"]
        """

        basic_info = [
            self.img_number,
            self.number,
            self.nucleus.nuc_volume,
            self.nucleus.nuc_cylinder_pix_num,
            self.nucleus.nuc_length,
            self.nucleus.nuc_width,
            self.nucleus.nuc_high_alternative
        ]

        channel_info = [
            [channel.av_signal_in_nuc_area_3D, channel.sum_pix_in_nuc_cylinder, channel.has_ring, channel.ring_intensity_coef]
            for channel in self.channels
        ]

        # Flatten channel_info and concatenate with basic_info
        return basic_info + [item for sublist in channel_info for item in sublist]
