import time
import pickle
import javabridge
import bioformats
import logging
import json
import os
from types import SimpleNamespace
from pathlib import Path

from afilament.objects.CellAnalyser import CellAnalyser
from afilament.objects.Parameters import ImgResolution, CellsImg
from afilament.objects import Utils

def main():

    # Specify image numbers to be analyzed
    img_nums = range(0, 23)

    # Start Java virtual machine for Bioformats library
    javabridge.start_vm(class_path=bioformats.JARS)

    # Initialize CellAnalyser object with configuration settings
    # There are some outputs for example "nuc area verification" or "cell_img_objects", that should be saved when we
    start = time.time()

    # Load JSON configuration file.
    with open("config.json", "r") as f:
        config = json.load(f, object_hook=lambda d: SimpleNamespace(**d))

    # Set up logging to record errors
    logging.basicConfig(filename='myapp.log', level=logging.DEBUG,
                        format='%(asctime)s %(levelname)s %(name)s %(message)s')
    logger = logging.getLogger(__name__)


    analyser = CellAnalyser(config)

    Utils.prepare_folder(config.imgs_objects)
    # Analyze each specified image and store cell data in all_cells list

    aggregated_stat_list = []
    channes = None

    for img_num in img_nums:
        try:
            cells, img_name = analyser.analyze_img(img_num)

            # Save analyzed image to a pickle file
            cells_img = CellsImg(img_name, analyser.img_resolution, cells)
            aggregated_stat_list = analyser.add_aggregated_cells_stat(aggregated_stat_list, cells_img.cells,
                                                                      cells_img.name)
            channes = cells_img.cells[0].channels

        except Exception as e:
            # Log error message if analysis fails for an image
            logger.error(f"\n----------- \n Img #{img_num} from file {config.confocal_img} was not analysed. "
                         f"\n Error: {e} \n----------- \n")
            print("An exception occurred")

    analyser.save_aggregated_cells_stat_list(aggregated_stat_list, channes)
    analyser.save_config(config.imgs_objects)

    end = time.time()
    print("Total time is: ")
    print(end - start)

    # Kill Java virtual machine
    javabridge.kill_vm()

if __name__ == '__main__':
    main()