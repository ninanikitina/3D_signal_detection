import os
import cv2
import numpy as np
import math
from skimage import measure
import trimesh
import pyvista as pv
import pymeshfix as mf

from objects import Utils
from objects import Contour
from unet.predict import run_predict_unet


class Nucleus(object):
    """
    Creates an object that track all nucleus information
    """

    def __init__(self, nuc_max_projection_mask):
        self.nuc_max_projection_mask = nuc_max_projection_mask
        self.nuc_volume = None
        self.nuc_length = None
        self.nuc_width = None
        self.nuc_high = None
        self.nuc_high_alternative = None
        self.nuc_3D_mask = None
        self.point_cloud = None
        self.nucleus_3d_img = None
        self.nuc_3D_whole_img = None #This image depict whole picture
        self.nuc_intensity = None
        self.nuc_cylinder_pix_num = None


    def reconstruct(self, temp_folders, unet_parm, resolution, analysis_folder, norm_th, cnt_extremes):
        """
        Reconstruct nucleus and saves all stat info into this Nicleus objects
        ---
        Parameters:
        - rot_angle (int): rotation angle to rotate images of the nucleus before making cross-section images for segmentation
        - rotated_cnt_extremes (CntExtremes object):
        - folders: (list): list of folders to save intermediate images during analysis
        - unet_parm (UnetParam object):
        img_resolution (ImgResolution object):
        """
        print("\nFinding nuclei...")
        Utils.prepare_folder(temp_folders["nucleous_xsection"])
        Utils.prepare_folder(temp_folders["nucleous_xsection_unet"])
        Utils.prepare_folder(temp_folders["nucleus_mask"])

        if len(os.listdir(temp_folders["cut_out_nuc"])) == 0:
            raise ValueError(f"Can not reconstruct nucleus, the directory {temp_folders['cut_out_nuc']} that should "
                             f"include preprocessed images is empty")
        nucleus_3d_img, _ = Utils.get_3D(temp_folders["cut_out_nuc"], 'nucleus')


        Utils.get_yz_xsection(nucleus_3d_img, temp_folders["nucleous_xsection"], 'nucleus', cnt_extremes)
        Utils.save_as_8bit(temp_folders["nucleous_xsection"], temp_folders["nucleous_xsection_unet"], norm_th)
        run_predict_unet(temp_folders["nucleous_xsection_unet"], temp_folders["nucleus_mask"], unet_parm.nucleus_unet_model,
                         unet_parm.unet_model_scale,
                         unet_parm.unet_model_thrh)
        self.nuc_3D_whole_img = nucleus_3d_img
        self.nucleus_3d_img = Utils.get_3d_img(temp_folders["nucleous_xsection"])
        self.nuc_3D_mask = Utils.get_3d_img(temp_folders["nucleus_mask"])
        self.nucleus_reco_3d(resolution, analysis_folder)
        self.nuc_length = (cnt_extremes.right[0] - cnt_extremes.left[0]) * resolution.x
        self.nuc_width = (cnt_extremes.bottom[1] - cnt_extremes.top[1]) * resolution.y


        point_cloud = np.array(self.point_cloud)
        print("Debug Info:")
        print("point_cloud shape:", point_cloud.shape)
        print("point_cloud size:", point_cloud.size)
        print("point_cloud dtype:", point_cloud.dtype)
        self.nuc_high_alternative = (max(point_cloud[:, 2]) - min(point_cloud[:, 2])) * resolution.z
        self.nuc_high = 2 * self.nuc_volume * 3/4 / (math.pi * self.nuc_length/2 * self.nuc_width/2)
        self.nuc_cylinder_pix_num = np.count_nonzero(self.nuc_max_projection_mask) * nucleus_3d_img.shape[2]


    def  get_cut_off_z(self, cap_bottom_ratio):
        nuc_point_cloud = np.array(self.point_cloud)
        nuc_z_min = min(nuc_point_cloud[:, 2])
        nuc_z_max = max(nuc_point_cloud[:, 2])
        cut_off_z = cap_bottom_ratio * (nuc_z_max - nuc_z_min) + nuc_z_min
        return cut_off_z

    def nucleus_reco_3d(self, resolution, analysis_folder):
        points = []

        xdata, ydata, zdata = [], [], []
        volume, intensity = 0, 0
        print(f"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
              f"Debugging Nucleus 3 D mask:\n"
              f"nuc_3D_mask.shape: {self.nuc_3D_mask.shape}")
        i = 0
        for slice in range(self.nuc_3D_mask.shape[0]):
            xsection_mask = self.nuc_3D_mask[slice, :, :]
            xsection_img = self.nucleus_3d_img[slice, :, :]

            xsection_mask[xsection_mask == 255] = 1
            xsection_img = np.multiply(xsection_img, xsection_mask)
            intensity += np.sum(xsection_img, dtype=np.int64)

            slice_cnts = cv2.findContours(xsection_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)[0]
            print(f"Debuging slice_cnts {i} info: {slice_cnts} \n")
            i += 1
            if len(slice_cnts) != 0:
                slice_cnt = slice_cnts[np.argmax([len(cnt) for cnt in slice_cnts])]
                volume += cv2.contourArea(
                    slice_cnt) * resolution.x * resolution.y * resolution.z  # calculate volume based on counter area
                # volume += np.count_nonzero(utils.draw_cnts((512, 512), [slice_cnt])) * scale_x * scale_y * scale_z #calculate volume based on pixel number

                if slice % 15 == 0:
                    ys = [pt[0, 0] for idx, pt in enumerate(slice_cnt) if idx % 4 == 0 and pt[0, 0] < 720]  # 720
                    zs = [pt[0, 1] for idx, pt in enumerate(slice_cnt) if idx % 4 == 0 and pt[0, 0] < 720]

                    xdata.extend([slice] * len(ys))
                    ydata.extend(ys)
                    zdata.extend(zs)

                cnt_ys = slice_cnt[:, 0, 1]
                cnt_zs = slice_cnt[:, 0, 0]

                points.extend([[x, y, z] for x, y, z in
                               zip([slice] * len(cnt_ys), cnt_ys, cnt_zs)])

        print("Nucleus volume: {}".format(volume))

        # TODO: to print out nucleus data commented this out
        # ax = plt.axes(projection='3d')
        # ax.scatter3D(ydata, zdata, xdata, cmap='Greens', alpha=0.5)
        # plt.show()
        self.nuc_volume = volume
        self.nuc_intensity = intensity
        self.point_cloud = points
        # print(f"###########################POINTS###########################################\n\n\n"
        #       f"{points}\n\n\n")

    def get_nucleus_origin(self):
        """
        Finds coordinate of nucleus anchor position which is:
        x is x of the center of the nucleus
        y is y of the center of the nucleus
        z is z the bottom of the nucleus
        ---
            Returns:
            - center_x, center_y, center_z (int, int, int) coordinates of the nucleus anchor position
        """
        center_x = self.nuc_3D_mask.shape[0] // 2
        slice_cnts = cv2.findContours(self.nuc_3D_mask[center_x, :, :], cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)[0]
        slice_cnt = slice_cnts[np.argmax([len(cnt) for cnt in slice_cnts])]
        cnt_extremes = Contour.get_cnt_extremes(slice_cnt)
        center_y = self.nuc_3D_mask.shape[1] // 2
        center_z = cnt_extremes.right[0]

        return center_x, center_y, center_z

    def save_nucleus_mesh(self, resolution, folder_name, file_name):

        image_3d = self.nuc_3D_mask

        # Add padding to the image
        p = 1
        z0 = np.zeros((p, image_3d.shape[1], image_3d.shape[2]), dtype=image_3d.dtype)
        z1 = np.zeros((image_3d.shape[0] + p * 2, p, image_3d.shape[2]), dtype=image_3d.dtype)
        image_3d = np.concatenate((z1, np.concatenate((z0, image_3d, z0), axis=0), z1), axis=1)

        # Generate mesh using marching cubes
        verts, faces, normals, values = measure.marching_cubes(image_3d, 0, step_size=3, allow_degenerate=False)

        # Flip the z-axis to orient the mesh correctly
        for v in verts:
            v[2] *= -1
        surf_mesh = trimesh.Trimesh(verts, faces, validate=True)

        meshfix = mf.MeshFix(pv.wrap(surf_mesh))
        meshfix.repair(verbose=True)
        fixed_mesh = meshfix.mesh

        # Scale the mesh to match the image resolution
        scaled_mesh = fixed_mesh.scale([1.0, 1.0, resolution.z / resolution.y], inplace=False)

        # Smooth the mesh and subdivide it
        smoothed_mesh = scaled_mesh.smooth(n_iter=400)
        final_mesh = smoothed_mesh.subdivide(1, subfilter="loop")


        if not os.path.exists(f'mesh/{folder_name}'):
            os.makedirs(f'mesh/{folder_name}')

        final_mesh.save(f'mesh/{folder_name}/{file_name}.stl')


    def get_nucleus_mesh(self, resolution):
        image_3d = self.nuc_3D_mask

        # Add padding to the image
        p = 1
        z0 = np.zeros((p, image_3d.shape[1], image_3d.shape[2]), dtype=image_3d.dtype)
        z1 = np.zeros((image_3d.shape[0] + p * 2, p, image_3d.shape[2]), dtype=image_3d.dtype)
        image_3d = np.concatenate((z1, np.concatenate((z0, image_3d, z0), axis=0), z1), axis=1)

        # Generate mesh using marching cubes
        verts, faces, normals, values = measure.marching_cubes(image_3d, 0, step_size=3, allow_degenerate=False)

        # Flip the z-axis to orient the mesh correctly
        for v in verts:
            v[2] *= -1
        surf_mesh = trimesh.Trimesh(verts, faces, validate=True)

        meshfix = mf.MeshFix(pv.wrap(surf_mesh))
        meshfix.repair(verbose=True)
        fixed_mesh = meshfix.mesh

        # Scale the mesh to match the image resolution
        scaled_mesh = fixed_mesh.scale([1.0, 1.0, resolution.z / resolution.y], inplace=False)

        # Smooth the mesh and subdivide it
        smoothed_mesh = scaled_mesh.smooth(n_iter=400)
        final_mesh = smoothed_mesh.subdivide(1, subfilter="loop")
        return final_mesh




