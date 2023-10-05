import pickle
from afilament.objects.Parameters import TestStructure
from afilament.objects.Fibers import Fibers

if __name__ == '__main__':
    # actin_obj_bottom= pickle.load(open(r"D:\BioLab\Current_experiments\afilament\2023.02.14_DAPI_Alexa488_LIV_Experiment\5_models_analysis_outcome_for_2imgs\1_training_set_img_CP_epoch200_W20\actin_objects\img_num_1__cell_num_0_bottom_3d_actin.obj", "rb"))
    # test_structure = pickle.load(open('test_structure.pickle', "rb"))
    actin_obj_cap= pickle.load(open(r"D:\BioLab\Current_experiments\afilament\2023.02.14_DAPI_Alexa488_LIV_Experiment\2023.03.14_analysis_data_W200__min-fiber-40\Control\analysis_data_Control_cells_0-1-3-15\actin_objects\img_num_0__cell_num_0_bottom_3d_actin.obj", "rb"))

    a = 1
    actin_obj_cap.plot_old()
    # actin_obj_bottom.plot()


    a = 1
    #
    # merged_fibers = Fibers("whole")
    # merged_fibers.merge_fibers(test_structure.fibers, test_structure.nodes, test_structure.pairs, test_structure.resolution)
    # single_fibers_num = len(test_structure.fibers.fibers_list)
    # merged_fibers_num = len(merged_fibers.fibers_list)
    # single_fibers_num_in_merged = 0
    # for merged_fiber in merged_fibers.fibers_list:
    #     single_fibers_num_in_merged += len(merged_fiber.fibers)
    # print(f"single_fibers_num is {single_fibers_num}")
    # print(f"merged_fibers_num is {merged_fibers_num}")
    # print(f"single_fibers_num_in_merged is {single_fibers_num_in_merged}")
    # merged_fibers.plot()
