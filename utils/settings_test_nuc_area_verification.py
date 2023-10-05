import time
import javabridge
import bioformats
import logging
import json
from types import SimpleNamespace

from afilament.objects.CellAnalyser import CellAnalyser

def main():
    # Load JSON configuration file. This file can be produced by GUI in the future implementation
    with open("../afilament/config.json", "r") as f:
        config = json.load(f, object_hook=lambda d: SimpleNamespace(**d))

    img_nums = range(0, 152)
    print(img_nums)

    javabridge.start_vm(class_path=bioformats.JARS)
    analyser = CellAnalyser(config)
    start = time.time()
    logging.basicConfig(filename='../afilament/myapp.log', level=logging.DEBUG,
                        format='%(asctime)s %(levelname)s %(name)s %(message)s')
    logger = logging.getLogger(__name__)

    for img_num in img_nums:
        try:
            analyser.save_nuc_verification(img_num, output_folder=r"C:\Users\nnina\Desktop\analysis_data")
        except Exception as e:
            logger.error(f"\n----------- \n Img #{img_num} from file {config.confocal_img} was not analysed. "
                                 f"\n Error: {e} \n----------- \n")
            print("An exception occurred")

    end = time.time()
    print("Total time is: ")
    print(end - start)
    javabridge.kill_vm()


if __name__ == '__main__':
    main()