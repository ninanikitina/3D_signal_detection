import pickle
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

if __name__ == '__main__':
    obj_cap= pickle.load(open(r"D:\BioLab\scr_2.0\afilament\analysis_data\actin_objects\img_num_0__cell_num_3_whole_3d_actin.obj", "rb"))

    df = pd.DataFrame(columns=['Layers'])
    for i, fiber in enumerate(obj_cap.fibers_list):
        fiber_df = pd.DataFrame({'Layers': fiber.last_layer,
                f'Intensity_{i}': fiber.intensities})
        df = pd.merge(df, fiber_df, how="outer", on=["Layers"])

    sns.lineplot(x='Layers', y='value', hue='variable',
                 data=pd.melt(df, ['Layers']), legend=False)
    plt.show()



