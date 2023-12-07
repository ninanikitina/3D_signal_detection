import os
import cv2
from tqdm import tqdm
from pathlib import Path
import bioformats

from objects.Parameters import ImgResolution


class ConfocalImgReader(object):
    """
    Creates an object that reads confocal microscopy images of two channels (actin and nucleus)
    """

    def __init__(self, path, nucleus_channel_name, cell_number, norm_th):
        """
            Parameters:
            img_path (string): path to the file to read
            nucleus_channel(int): channel of nucleus images at the provided microscopic image
            actin_channel(int): channel of actin images at the provided microscopic image
        """
        self.image_path, self.series = self.get_img_path_and_series(path, cell_number)
        self.norm_th = norm_th
        metadata = bioformats.get_omexml_metadata(str(self.image_path))
        self.metadata_obj = bioformats.OMEXML(metadata)
        self.channel_nums = self.metadata_obj.image(self.series).Pixels.get_channel_count()
        self.nuc_channel = self.find_channel(nucleus_channel_name)
        self.img_resolution = self.get_resolution()
        self.depth = self.metadata_obj.image(self.series).Pixels.PixelType  # Where depth (8 bit, etc.) is identified?
        # What does this actually return?
        self.t_num = self.metadata_obj.image(self.series).Pixels.SizeT
        self.magnification = self.metadata_obj.instrument(
            self.series).Objective.NominalMagnification  # "63.0" for 63x, "20.0" for 20x


    def find_channel(self, channel_name):
        channel_num = None
        for i in range(self.channel_nums):
            if channel_name in self.metadata_obj.image().Pixels.Channel(i).get_Name():
                channel_num = i

        return channel_num

    def get_resolution(self):
        scale_x = self.metadata_obj.image(self.series).Pixels.get_PhysicalSizeX()
        scale_y = self.metadata_obj.image(self.series).Pixels.get_PhysicalSizeY()
        scale_z = self.metadata_obj.image(self.series).Pixels.get_PhysicalSizeZ()
        img_resolution = ImgResolution(scale_x, scale_y, scale_z)
        return img_resolution

    def get_img_path_and_series(self, path, cell_number):
        """
        CZI and LIF files, in our case organized differently.
        LIF is a project file that has different images as a Series.
        CZI is a path to the folder that contains separate images.
        This method checks what is the case and finds the path-specific image and Series.
        Args:
            path: str, path to folder or project file

        Returns:
            img_path: path to file
            series: series to analyze
        """
        img_path = None
        if os.path.isfile(path):
            series = cell_number
            img_path = path

        else:
            series = 0
            folder_path = path
            for i, current_path in enumerate(Path(folder_path).rglob('*.czi')):
                if i == cell_number:
                    img_path = current_path
                    break
        print(img_path)

        return img_path, series




    def read_nucleus_layers(self, output_folder):
        """
        Converts confocal microscopic images into a set of png images specified in reader object

        """
        channel = self.nuc_channel
        image_stack = []
        z_layers_num = self.metadata_obj.image(self.series).Pixels.get_SizeZ()
        for i in range(z_layers_num):
            img = bioformats.load_image(str(self.image_path), c=channel, z=i, t=0, series=self.series, index=None, rescale=False,
                                        wants_max_intensity=False,
                                        channel_names=None)

            type = self.metadata_obj.image(self.series).Pixels.get_PixelType()

            image_stack.append(img)

        self._save_img(image_stack, "nucleus", output_folder)

        return image_stack[0].shape

    def read_all_layers(self, output_folder):
        """
        Converts confocal microscopic images into a set of png images specified in reader object

        """
        channel_names = []
        for channel in range(self.channel_nums):
            channel_name = self.metadata_obj.image().Pixels.Channel(channel).get_Name()
            channel_names.append(channel_name)
            image_stack = []
            z_layers_num = self.metadata_obj.image(self.series).Pixels.get_SizeZ()
            for i in range(z_layers_num):
                img = bioformats.load_image(str(self.image_path), c=channel, z=i, t=0, series=self.series, index=None,
                                            rescale=False,
                                            wants_max_intensity=False,
                                            channel_names=None)


                image_stack.append(img)

            self._save_img(image_stack, channel_name, output_folder)
        return channel_names




    def _save_img(self, img_stack, bio_structure, output_folder):
        for i in tqdm(range(len(img_stack))):
            img_name = os.path.splitext(os.path.basename(self.image_path))[0]
            img_path = os.path.join(output_folder,
                                    img_name + '_series_' + str(self.series) +'_' + bio_structure + '_layer_' + str(i) +'.png')
            cv2.imwrite(img_path, img_stack[i])

    def close(self):
        bioformats.clear_image_reader_cache()
