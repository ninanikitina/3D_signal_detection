import argparse
import logging
import os
import glob

import numpy as np
import torch
import torch.nn.functional as F
from torchvision import transforms
from PIL import Image
from tqdm import tqdm

from unet import UNet
from unet.data_vis import plot_img_and_mask
from unet.dataset import BasicDataset
os.environ["KMP_DUPLICATE_LIB_OK"] ="TRUE"


def predict_img(net,
                full_img,
                device,
                scale_factor=1,
                out_threshold=0.5):
    net.eval()

    img = torch.from_numpy(BasicDataset.preprocess(full_img, scale_factor))

    img = img.unsqueeze(0)
    img = img.to(device=device, dtype=torch.float32)

    with torch.no_grad():
        output = net(img)

        if net.n_classes > 1:
            probs = F.softmax(output, dim=1)
        else:
            probs = torch.sigmoid(output)

        probs = probs.squeeze(0)

        tf = transforms.Compose(
            [
                transforms.ToPILImage(),
                transforms.Resize(full_img.size[1]),
                transforms.ToTensor()
            ]
        )

        probs = tf(probs.cpu())
        full_mask = probs.squeeze().cpu().numpy()

    return full_mask > out_threshold


def get_args():
    parser = argparse.ArgumentParser(description='Predict masks from input images',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--model', '-m', default='MODEL.pth',
                        metavar='FILE',
                        help="Specify the file in which the model is stored")
    parser.add_argument('--input', '-i', metavar='INPUT',
                        help='folder name of input imagies', required=True)
    parser.add_argument('--output', '-o', metavar='INPUT',
                        help='folder name of ouput images')
    parser.add_argument('--viz', '-v', action='store_true',
                        help="Visualize the images as they are processed",
                        default=False)
    parser.add_argument('--decision-threshold', '-t', type=float,
                        help="Minimum probability value to consider a mask pixel white",
                        default=0.5)
    parser.add_argument('--scale', '-s', type=float,
                        help="Scale factor for the input images",
                        default=0.5)

    return parser.parse_args()


def get_output_filenames(args):
    in_files_folder = args.input
    in_files = glob.glob(os.path.join(in_files_folder, "*.png"))
    out_files = []

    if not args.output:
        for f in in_files:
            pathsplit = os.path.splitext(f)
            out_files.append("{}_OUT{}".format(pathsplit[0], pathsplit[1]))
    else:
        out_files_folder = args.output
        for f in in_files:
            input_file_basename = os.path.basename(f)
            out_files.append(os.path.join(out_files_folder, input_file_basename))

    return out_files


def mask_to_image(mask):
    return Image.fromarray((mask * 255).astype(np.uint8))


def run_predict_unet(folder_path, output_folder_path, model_path, scale, decision_threshold, viz=False):
    if folder_path is None:
        if not os.path.exists('temp/true_nucleus_imgs'):
            raise RuntimeError("There is no folder {}\nCan't process images".format("temp/true_nucleus_imgs"))

        if not os.path.exists('temp/true_nucleus_masks'):
            os.makedirs('temp/true_nucleus_masks')

        folder_path = 'temp/true_nucleus_imgs'
        output_folder_path = 'temp/true_nucleus_masks'

    net = UNet(n_channels=1, n_classes=1)

    logging.info("Loading model {}".format(model_path))

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    logging.info(f'Using device {device}')
    net.to(device=device)
    net.load_state_dict(torch.load(model_path, map_location=device))

    logging.info("Model loaded !")

    for img_path in tqdm(glob.glob(os.path.join(folder_path, "*.png"))):
        logging.info("\nPredicting image {} ...".format(img_path))

        img = Image.open(img_path).convert('L')  # .convert('RGB')

        mask = predict_img(net=net,
                           full_img=img,
                           scale_factor=scale,
                           out_threshold=decision_threshold,
                           device=device)

        img_name = os.path.basename(img_path)
        img_path_to_save = os.path.join(output_folder_path, img_name)
        mask_to_image(mask).save(img_path_to_save)

        if viz:
            logging.info("Visualizing results for image {}, close to continue ...".format(img_path))
            plot_img_and_mask(img, mask)

def run_predict_unet_one_img(img_path, model_path, scale, decision_threshold):

    net = UNet(n_channels=1, n_classes=1)

    logging.info("Loading model {}".format(model_path))

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    logging.info(f'Using device {device}')
    net.to(device=device)
    net.load_state_dict(torch.load(model_path, map_location=device))

    logging.info("Model loaded !")
    img = Image.open(img_path).convert('L')
    mask = predict_img(net=net,
                      full_img=img,
                      scale_factor=scale,
                      out_threshold=decision_threshold,
                      device=device)

    return mask



if __name__ == "__main__":
    args = get_args()
    in_files_folder = args.input
    out_files = get_output_filenames(args)

    run_predict_unet(args.input, args.output, args.model, args.scale, args.decision_threshold, args.viz)
