import os
import time
import javabridge
import bioformats
import logging
import json
from types import SimpleNamespace
import argparse
import sys
import glob
sys.path.insert(0, 'D:/BioLab/src_3D_signal_detection')

from objects.CellAnalyser import CellAnalyser
from objects.Parameters import CellsImg

def load_config(config_path, img_path=None):
    with open(config_path, "r") as f:
        config = json.load(f, object_hook=lambda d: SimpleNamespace(**d))
    if img_path:
        config.confocal_img = img_path
    return config


def setup_logging():
    log_directory = "err_logs"
    log_file_name = "myapp.log"

    # Create 'err_logs' directory if it doesn't exist
    if not os.path.exists(log_directory):
        os.makedirs(log_directory)

    # Set up logging with the correct file path
    log_file_path = os.path.join(log_directory, log_file_name)
    logging.basicConfig(filename=log_file_path, level=logging.DEBUG,
                        format='%(asctime)s %(levelname)s %(name)s %(message)s')
    return logging.getLogger(__name__)


class JavaVM:
    def __enter__(self):
        # Set the directory for JVM crash logs
        crash_log_directory = os.path.abspath("err_logs")
        os.environ['JAVA_TOOL_OPTIONS'] = f'-XX:ErrorFile={crash_log_directory}/hs_err_pid%p.log'

        javabridge.start_vm(class_path=bioformats.JARS)
        javabridge.static_call("loci/common/DebugTools", "setRootLevel", "(Ljava/lang/String;)V", "OFF")
        return self


    def __exit__(self, exc_type, exc_value, traceback):
        javabridge.kill_vm()


def parse_arguments():
    parser = argparse.ArgumentParser(description="Cell Analysis Script")
    parser.add_argument('img_path', type=str, nargs='?', default=None, help='Path to the confocal image file')
    args = parser.parse_args()
    return args.img_path

def main():

    img_path = parse_arguments()
    config = load_config("config.json", img_path)
    logger = setup_logging()

    with JavaVM():
        start = time.time()
        analyser = CellAnalyser(config)

        aggregated_stat_list = []
        channels = None

        o = glob.glob(os.path.join(config.confocal_img, '*.czi'))


        for img_num, czi_file in enumerate(glob.glob(os.path.join(config.confocal_img, '*.czi'))):
            try:
                cells, img_name = analyser.analyze_img(img_num)
                cells_img = CellsImg(img_name, analyser.img_resolution, cells)
                if len(cells) > 0:
                    aggregated_stat_list = analyser.add_aggregated_cells_stat(aggregated_stat_list, cells_img.cells,
                                                                              cells_img.name)
                channels = cells_img.cells[0].channels

                if len(aggregated_stat_list) > 0:
                    analyser.save_aggregated_cells_stat_list(aggregated_stat_list, channels)

                analyser.save_config("")

            except Exception as e:
                logger.error(
                    f"\n----------- \n Img #{img_num} from file {config.confocal_img} was not analysed. \n Error: {e} \n----------- \n")
                print("An exception occurred")

        end = time.time()
        print(f"Total time is: {end - start}")


if __name__ == '__main__':
    main()
