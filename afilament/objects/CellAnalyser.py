import os
import csv
import json
import cv2
from pathlib import Path
from datetime import datetime
import numpy as np

from objects import Utils
from objects.ConfocalImgReader import ConfocalImgReader
from objects import Contour
from objects.Cell import Cell
from objects.Parameters import UnetParam


class CellAnalyser(object):
    def __init__(self, config):
        self.initial_conf = config
        self.nucleus_channel_name = config.nucleus_channel_name
        self.confocal_path = config.confocal_img
        self.nuc_theshold = config.nuc_theshold
        self.unet_parm = UnetParam(config)
        self.norm_th = config.norm_th
        self.img_resolution = None
        self.nuc_area_min_pixels_num = config.nuc_area_min_pixels_num
        self.total_img_number = 0
        self.total_cells_number = 0
        self.find_biggest_mode = config.find_biggest_mode
        # Set output_data_folder to include confocal_dir_name
        self.output_data_folder = config.output_analysis_path
        self.channel_names = None

        confocal_dir_name = Path(self.confocal_path).stem
        self.output_data_folder = config.output_analysis_path + f"/{confocal_dir_name}"

        self.temp_folders = {
            "raw": f'../afilament/img_data/{confocal_dir_name}/temp/czi_layers',
            "nuc_raw": f'../afilament/img_data/{confocal_dir_name}/temp/czi_layers_nuc',
            "cut_out_nuc": f'../afilament/img_data/{confocal_dir_name}/temp/actin_and_nucleus_layers',
            "cut_out_channel": f'../afilament/img_data/{confocal_dir_name}/temp/actin_and_nucleus_layers',
            "nucleous_xsection": f'../afilament/img_data/{confocal_dir_name}/temp/nucleus_layers',
            "nucleous_xsection_unet": f'../afilament/img_data/{confocal_dir_name}/temp/nucleus_layers_unet',
            "nucleus_mask": f'../afilament/img_data/{confocal_dir_name}/temp/nucleus_mask',
            "nucleus_top_mask": f'../afilament/img_data/{confocal_dir_name}/temp/nucleus_top_mask',
            "nucleus_top_img": f'../afilament/img_data/{confocal_dir_name}/temp/nucleus_top_img',
            "nuclei_top_masks": f'../afilament/img_data/{confocal_dir_name}/temp/nuclei_top_masks',
        }


        Utils.prepare_folder(self.output_data_folder)


    def analyze_img(self, img_num):

        # For metadata statistics
        self.total_img_number += 1

        for folder in self.temp_folders.values():
            Utils.prepare_folder(folder)

        reader = ConfocalImgReader(self.confocal_path, self.nucleus_channel_name, img_num, self.norm_th)
        self.img_resolution = reader.img_resolution
        mask_size = reader.read_nucleus_layers(self.temp_folders["nuc_raw"])
        self.channel_names = reader.read_all_layers(self.temp_folders["raw"])
        nuclei_masks, nuclei_xy_centers = Utils.get_nuclei_masks(self.temp_folders, self.output_data_folder,
                                              reader.image_path, self.nuc_theshold, self.nuc_area_min_pixels_num,
                                              self.find_biggest_mode, img_num, self.unet_parm)

        cells = []
        # Prepare a grayscale image to show the results.
        foci_image = np.zeros(mask_size, dtype=np.uint8)

        meshes = []

        for i, nuc_mask in enumerate(nuclei_masks):
            cell = self.run_analysis(img_num, i, nuc_mask, reader)

            folder_name = os.path.basename(os.path.dirname(reader.image_path))
            file_name = f"img-{os.path.basename(reader.image_path)}_cell-{i}"


            # cell.nucleus.save_nucleus_mesh(self.img_resolution, folder_name, file_name)
            meshes.append(cell.nucleus.get_nucleus_mesh(self.img_resolution))
            cells.append(cell)
            # foci_image = cell.count_foci(foci_image)
            self.total_cells_number += 1

        # #Save foci image to show Anamaria if foci is accuratelly found
        # img_base_path = os.path.splitext(os.path.basename(reader.image_path))[0]
        # foci_path = os.path.join(self.output_data_folder, "foci_" + img_base_path + ".png")
        # cv2.imwrite(foci_path, foci_image)

        folder_name = os.path.basename(os.path.dirname(reader.image_path))
        file_name = os.path.basename(reader.image_path)

        Utils.save_nuclei_meshes(meshes, nuclei_xy_centers, mask_size, folder_name, file_name)

        return cells, Path(reader.image_path).name


    def run_analysis(self, img_num, cell_num, nucleus_mask, reader):
        cell = Cell(img_num, cell_num, self.channel_names)

        Utils.сut_out_mask(nucleus_mask, self.temp_folders["nuc_raw"], self.temp_folders["cut_out_nuc"], 'nucleus')
        cnt_extremes = Contour.get_cnt_extremes(Contour.get_mask_cnt(nucleus_mask))
        cell.analyze_nucleus(self.temp_folders, self.unet_parm, self.img_resolution,
                             self.output_data_folder, self.norm_th, cnt_extremes, nucleus_mask)

        print(f"Cell #{cell_num}\n")
        cell.analyze_channels(self.temp_folders["raw"], cnt_extremes, nucleus_mask, self.nuc_theshold, self.nucleus_channel_name)
        return cell

    def add_aggregated_cells_stat(self, cell_stat_list, cells, img_name):
        for cell in cells:
            cell_stat_list.append([str(img_name)] + cell.get_aggregated_cell_stat())
        return cell_stat_list

    def save_aggregated_cells_stat_list(self, stat_list, channels):
        basic_info = ["Image_name", "Img_num", "Cell_num", "Nucleus_volume, cubic_micrometre",
                      "Nucleus_cylinder, pixels_number", "Nucleus_length, micrometre",
                      "Nucleus_width, micrometre", "Nucleus_high, micrometre"]

        channel_info = [
            [f"{channel.name} av_signal_in_nuc_area_3D", f"{channel.name} sum_pix_in_nuc_cylinder",
             f"{channel.name} has ring", f"{channel.name} ring intensity coef"]
            for channel in channels
        ]
        header_row = basic_info + [item for sublist in channel_info for item in sublist]

        path = os.path.join(self.output_data_folder, 'cell_stat.csv')
        with open(path, mode='w') as stat_file:
            csv_writer = csv.writer(stat_file, delimiter=',')
            csv_writer.writerow(header_row)
            for raw in stat_list:
                csv_writer.writerow(raw)

        print("Aggregated stat created")

    def save_config(self, img_folder):
        # Add additional information
        analysis_date = datetime.today().strftime('%Y-%m-%d %H:%M:%S')
        self.initial_conf.analysis_date_time = analysis_date
        self.initial_conf.total_img_number = self.total_img_number
        self.initial_conf.total_cells_number = self.total_cells_number
        # Serializing json
        json_conf_str = json.dumps(self.initial_conf, indent=4, default=lambda o: o.__dict__,
                                   sort_keys=True)

        # Writing to sample.json in analysis output folder
        file_path = os.path.join(self.output_data_folder, "analysis_configurations.json")
        with open(file_path, "w") as outfile:
            outfile.write(json_conf_str)

    def save_nuc_verification(self, img_num, output_folder):
        """
        Save nucleus area verificatin imagies. This function is helpful to verify different settings
        """
        for folder in self.temp_folders.values():
            Utils.prepare_folder(folder)


        reader = ConfocalImgReader(self.confocal_path, self.nucleus_channel_name, img_num, self.norm_th)
        reader.read_nucleus_layers(self.temp_folders["nuc_raw"])
        Utils.get_nuclei_masks(self.temp_folders, output_folder,
                               reader.image_path, self.nuc_theshold, self.nuc_area_min_pixels_num,
                               self.find_biggest_mode, img_num, self.unet_parm)


    def save_nuc_verification_and_mask(self, img_num, output_folder_ver, output_folder_masks):
        """
        Save nucleus area verificatin imagies. This function is helpful to verify different settings
        """
        for folder in self.temp_folders.values():
            Utils.prepare_folder(folder)


        reader = ConfocalImgReader(self.confocal_path, self.nucleus_channel_name, img_num, self.norm_th)
        reader.read_nucleus_layers(self.temp_folders["nuc_raw"])
        nuclei_masks = Utils.get_nuclei_masks(self.temp_folders, output_folder_ver,
                               reader.image_path, self.nuc_theshold, self.nuc_area_min_pixels_num,
                               self.find_biggest_mode, img_num, self.unet_parm)

        stacked_masks = np.stack(nuclei_masks, axis=0)
        sum_projection = np.sum(stacked_masks, axis=0).astype(dtype=np.uint8)

        img_base_path = os.path.splitext(os.path.basename(reader.image_path))[0]

        mask_path = os.path.join(output_folder_masks, img_base_path + ".png")
        cv2.imwrite(mask_path, sum_projection)

