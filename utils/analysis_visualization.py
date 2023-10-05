import os
import numpy as np
import trimesh
from skimage import measure
import pickle
import pyvista as pv
import pymeshfix as mf
from afilament.objects.Node import plot_branching_nodes

# Define constants
class StructureOptions:
    CAP = "cap"
    BOTTOM = "bottom"
    TOTAL = "total"

class VisualizationModes:
    NUCLEUS = "nucleus"
    ACTIN = "actin"

def get_nucleus_mesh(cell):
    """
    Generate a mesh of the nucleus using marching cubes algorithm.

    Args:
        cell: A `Cell` object.

    Returns:
        A `trimesh.Trimesh` object representing the nucleus mesh.
    """
    image_3d = cell.nucleus.nuc_3D_mask

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

    return surf_mesh


def optimize_mesh(mesh, resolution):
    """
    Perform mesh optimization by fixing errors, scaling, smoothing, and subdivision.

    Args:
        mesh: A `trimesh.Trimesh` object.
        resolution: A `Resolution` object representing the resolution of the image.

    Returns:
        A `pyvista.PolyData` object representing the optimized mesh.
    """
    # Fix errors in the mesh
    meshfix = mf.MeshFix(pv.wrap(mesh))
    meshfix.repair(verbose=True)
    fixed_mesh = meshfix.mesh

    # Scale the mesh to match the image resolution
    scaled_mesh = fixed_mesh.scale([1.0, 1.0, resolution.z / resolution.y], inplace=False)

    # Smooth the mesh and subdivide it
    smoothed_mesh = scaled_mesh.smooth(n_iter=400)
    final_mesh = smoothed_mesh.subdivide(1, subfilter="loop")

    return final_mesh


def visualise_actin(afilament_folder_path, image_index, cell_index, min_fiber_thr_microns, node_actin_len_th, show_branching_nodes, structure):
    """
    Visualize the actin fibers of a cell.

    Args:
        afilament_folder_path: A string representing the path to the afilament folder.
        image_index: An integer representing the index of the image to visualize.
        cell_index: An integer representing the index of the cell to visualize.
        min_fiber_thr_pixels: A float representing the minimum fiber length to include in the visualization.
        node_actin_len_th: A float representing the minimum length of actin fibers connecting branching nodes to include in the visualization.
        show_branching_nodes: A boolean indicating whether to highlight branching nodes in the visualization.
        structure: A string representing the actin structure to visualize. Must be one
    """

    # Check input variables
    if structure not in StructureOptions.__dict__.values():
        raise ValueError(f"Invalid structure option: {structure}. Must be one of: {StructureOptions.__dict__.values()}")

    # Load image data
    image_path = os.path.join(afilament_folder_path, f"image_data_{image_index}.pickle")
    with open(image_path, "rb") as f:
        cells_img = pickle.load(f)

    # Get cell object
    if cell_index >= len(cells_img.cells):
        raise ValueError(f"Invalid cell index: {cell_index}. Must be less than: {len(cells_img.cells)}")
    cell = cells_img.cells[cell_index]

    # This code block recalculates the length of the actin for the old format.
    # As per our decision for the LIV paper, we are aligning it with the method used in the KASH paper.
    # The formula used for calculating the old length is: old_len = (fiber.xs[-1] - fiber.xs[0]) * resolution.x
    cell.update_actin_stat_old_format(cells_img.resolution)


    # Get fiber statistics so we can put it on a graph. Fiber object will be updated
    if structure == StructureOptions.CAP:
        cell.actin_cap.create_fibers_aggregated_stat(min_fiber_thr_microns, cells_img.resolution)
    elif structure == StructureOptions.BOTTOM:
        cell.actin_bottom.create_fibers_aggregated_stat(min_fiber_thr_microns, cells_img.resolution)
    else:  # structure == StructureOptions.TOTAL
        cell.actin_total.create_fibers_aggregated_stat(min_fiber_thr_microns, cells_img.resolution)

    if show_branching_nodes:
        cell.find_branching(min_fiber_thr_microns, node_actin_len_th)

        # Generate visualization
        if structure == StructureOptions.CAP:
            plot_branching_nodes(cell.actin_cap, cell.cap_nodes, min_fiber_thr_microns, cells_img.resolution,
                                 structure, image_index, cell_index)
        elif structure == StructureOptions.BOTTOM:
            plot_branching_nodes(cell.actin_bottom, cell.bottom_nodes, min_fiber_thr_microns, cells_img.resolution,
                                 structure, image_index, cell_index)
        else:  # structure == StructureOptions.TOTAL
            plot_branching_nodes(cell.actin_total, cell.total_nodes, min_fiber_thr_microns, cells_img.resolution,
                                 structure, image_index, cell_index)

    else:

        # Generate visualization
        if structure == StructureOptions.CAP:
            cell.actin_cap.plot(min_fiber_thr_microns)
        elif structure == StructureOptions.BOTTOM:
            cell.actin_bottom.plot(min_fiber_thr_microns)
        else:  # structure == StructureOptions.TOTAL
            cell.actin_total.plot(min_fiber_thr_microns)


def visualise_nucleus(afilament_folder_path, image_index, cell_index):
    """
    Visualize the nucleus of a cell.

    Args:
        afilament_folder_path: A string representing the path to the afilament folder.
        image_index: An integer representing the index of the image to visualize.
        cell_index: An integer representing the index of the cell to visualize.
    """
    # Load image data
    image_path = os.path.join(afilament_folder_path, f"image_data_{image_index}.pickle")
    with open(image_path, "rb") as f:
        cells_img = pickle.load(f)

    # Get cell object
    if cell_index >= len(cells_img.cells):
        raise ValueError(f"Invalid cell index: {cell_index}. Must be less than: {len(cells_img.cells)}")
    cell = cells_img.cells[cell_index]

    # Generate and optimize mesh
    tmesh = get_nucleus_mesh(cell)
    final_mesh = optimize_mesh(tmesh, cells_img.resolution)

    # Create and show plot
    p = pv.Plotter()
    p.add_mesh(final_mesh, color=True, show_edges=True)

    p.add_text(f"Volume: {cell.nucleus.nuc_volume:.0f} \u03BCm^3\n\n"
               f"Length: {cell.nucleus.nuc_length:.2f} \u03BCm\n\n"
               f"Width: {cell.nucleus.nuc_width:.2f} \u03BCm\n\n"
               f"Height: {cell.nucleus.nuc_high_alternative:.2f} \u03BCm\n\n",
               position='right_edge',
               color='white',
               font_size=9)
    p.add_title(f"Nucleus of image # {image_index}, cell # {cell_index}", font_size=11)
    final_mesh.save(f'mesh/img_{image_index}__cell_{cell_index}.stl')

    p.show()


if __name__ == '__main__':
    afilament_folder_path = r"D:/BioLab/Current_experiments/afilament/2023.02.14_DAPI_Alexa488_LIV_Experiment/Second_trial/Control_w20/Initial run/img_objects"

    # afilament_folder_path = r"D:\BioLab\Current_experiments\afilament\2023.02.14_DAPI_Alexa488_LIV_Experiment\Second_trial\LIV_w20\Initial run\img_objects"
    image_index = 9
    cell_index = 3
    min_fiber_thr_microns = 7
    node_actin_len_th = 2
    show_branching_nodes = True
    structure = StructureOptions.CAP
    vis_mode = VisualizationModes.ACTIN

    if vis_mode not in VisualizationModes.__dict__.values():
        raise ValueError(
            f"Invalid visualization mode: {vis_mode}. Must be one of: {VisualizationModes.__dict__.values()}")

    if vis_mode == VisualizationModes.ACTIN:
        visualise_actin(afilament_folder_path, image_index, cell_index, min_fiber_thr_microns, node_actin_len_th,
                        show_branching_nodes, structure)
    elif vis_mode == VisualizationModes.NUCLEUS:
        visualise_nucleus(afilament_folder_path, image_index, cell_index)