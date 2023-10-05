import tkinter as tk
from tkinter.messagebox import showinfo
from tkinter import filedialog
import customtkinter as ctk
import javabridge
import bioformats
import time
import pickle
import logging
import json
import os
from types import SimpleNamespace
from pathlib import Path

from afilament.objects.Parameters import ImgResolution, CellsImg
from afilament.objects.CellAnalyser import CellAnalyser


def run_through_gui_recalculation(input_obj_path, output_stat_folder, fiber_min_layers_theshold):

    start = time.time()
    config_file_path = os.path.join(input_obj_path, "analysis_configurations.json")
    with open(config_file_path, "r") as f:
        config = json.load(f, object_hook=lambda d: SimpleNamespace(**d))
    config.output_analysis_path = output_stat_folder
    config.fiber_min_thr_pixels = fiber_min_layers_theshold
    config.imgs_objects = input_obj_path
    analyser = CellAnalyser(config)

    # Extract statistical data
    # Extract all_cells from images data
    aggregated_stat_list = []
    for file in os.listdir(config.imgs_objects):
        # add check if directory is empty and ask user to specify where get data
        img_path = Path(os.path.join(config.imgs_objects, file))

        # check if it is image file since in this folder we have config file
        if img_path.suffix == ".pickle":
            cells_img = pickle.load(open(img_path, "rb"))
            analyser.img_resolution = cells_img.resolution
            # Save individual cell data to CSV file
            analyser.save_cells_data(cells_img.cells)
            aggregated_stat_list = analyser.add_aggregated_cells_stat(aggregated_stat_list, cells_img.cells)

    # Save aggregated cell statistics to CSV file
    analyser.save_aggregated_cells_stat_list(aggregated_stat_list)
    # # Save current configuration settings to JSON file
    analyser.save_config(is_recalculated=False, img_folder=input_obj_path)

    end = time.time()
    print("Total time is: ")
    print(end - start)


def run_through_gui_new_analysis(bioformat_imgs_path, output_stat_folder,
                                 nuc_channel, actin_channel,
                                 fiber_min_layers_theshold, output_img_obj_folder):
    # Load JSON configuration file. This file can be produced by GUI in the future implementation
    with open("config.json", "r") as f:
        config = json.load(f, object_hook=lambda d: SimpleNamespace(**d))

    # Add user input into configuration file:
    config.actin_channel = actin_channel
    config.nucleus_channel_name = nuc_channel
    config.fiber_min_thr_pixels = fiber_min_layers_theshold
    config.confocal_img = bioformat_imgs_path
    config.output_analysis_path = output_stat_folder
    config.imgs_objects = output_img_obj_folder
    # Specify image numbers to be analyzed
    img_nums = range(1, 5)


    # Start Java virtual machine for Bioformats library
    javabridge.start_vm(class_path=bioformats.JARS)

    # Initialize CellAnalyser object with configuration settings
    analyser = CellAnalyser(config)

    start = time.time()

    # Set up logging to record errors
    logging.basicConfig(filename='myapp.log', level=logging.DEBUG,
                        format='%(asctime)s %(levelname)s %(name)s %(message)s')
    logger = logging.getLogger(__name__)

    for img_num in img_nums:
        try:
            cells = analyser.analyze_img(img_num)

            # Save analyzed image to a pickle file
            image_data_path = os.path.join(config.imgs_objects, "image_data_" + str(img_num) + ".pickle")
            cells_img = CellsImg(analyser.img_resolution, cells)
            with open(image_data_path, "wb") as file_to_save:
                pickle.dump(cells_img, file_to_save)

        except Exception as e:
            # Log error message if analysis fails for an image
            logger.error(f"\n----------- \n Img #{img_num} from file {config.confocal_img} was not analysed. "
                         f"\n Error: {e} \n----------- \n")
            print("An exception occurred")

        # Extract statistical data
        # Extract all_cells from images data
    aggregated_stat_list = []
    for file in os.listdir(config.imgs_objects):
        # add check if directory is empty and ask user to specify where get data
        img_path = Path(os.path.join(config.imgs_objects, file))

        # check if it is image file since in this folder we have config file
        if img_path.suffix == ".pickle":
            cells_img = pickle.load(open(img_path, "rb"))
            analyser.img_resolution = cells_img.resolution
            # Save individual cell data to CSV file
            analyser.save_cells_data(cells_img.cells)
            aggregated_stat_list = analyser.add_aggregated_cells_stat(aggregated_stat_list, cells_img.cells)

    # Save aggregated cell statistics to CSV file
    analyser.save_aggregated_cells_stat_list(aggregated_stat_list)
    # # Save current configuration settings to JSON file
    analyser.save_config(is_recalculated=False, img_folder=config.imgs_objects)

    end = time.time()
    print("Total time is: ")
    print(end - start)

    # Kill Java virtual machine
    javabridge.kill_vm()



def show_data_page():
    type = analysis_page.selected_type.get()
    print(analysis_page.selected_type.get())

    if type == "new":
        print("New analysis")
        new_analysis_frame.tkraise()
    elif type == "recalculate":
        recalculate_frame.tkraise()

    window_width = 370
    window_height = 430

    # get the screen dimension
    screen_width = root.winfo_screenwidth()
    screen_height = root.winfo_screenheight()

    # find the center point
    center_x = int(screen_width / 2 - window_width / 2)
    center_y = int(screen_height / 2 - window_height / 2)

    root.geometry(f'{window_width}x{window_height}+{center_x}+{center_y}')




def run_new_analysis():

    print(f"Input folder: {data_page.input_obj_folder.get()}")
    bioformat_imgs_path = data_page.input_obj_folder.get()

    print(f"Output stat save to: {data_page.output_data_folder.get()}")
    output_stat_folder = data_page.output_data_folder.get()

    print(f"Output img objects save to: {data_page.output_obj_folder.get()}")
    output_objects_folder = data_page.output_obj_folder.get()

    print(f"Mask channel: {data_page.nuc_channel.get()}")
    nuc_channel = data_page.nuc_channel.get()

    print(f"Mask channel: {data_page.actin_channel.get()}")
    actin_channel = data_page.actin_channel.get()

    print(f"Minimal fiber length: {data_page.min_fiber_len.get()}")
    fiber_min_layers_theshold = data_page.min_fiber_len.get()

    if bioformat_imgs_path == "" or output_stat_folder == "" or fiber_min_layers_theshold == "":
        showinfo(
            title='Information',
            message="Please fill all blanks to make the program run"
        )

    else:
        showinfo(
            title='Information',
            message=f"The analysis process is about to begin. Once completed, you can find the results "
                    f"in the designated output folder: {output_stat_folder}. Please note that the analysis "
                    f"can take a significant amount of time, depending on the size of the dataset and "
                    f"the capabilities of your machine."

            )
        try:
            root.destroy()
            run_through_gui_new_analysis(bioformat_imgs_path, output_stat_folder,
                                         int(nuc_channel), int(actin_channel),
                                         int(fiber_min_layers_theshold), output_objects_folder)
        except Exception as e:
            a = 1


def run_recalculation():

    print(f"Input folder: {rec_page.input_obj_folder.get()}")
    input_obj_path = rec_page.input_obj_folder.get()

    print(f"Output stat save to: {rec_page.output_data_folder.get()}")
    output_stat_folder = rec_page.output_data_folder.get()

    print(f"Minimal fiber length: {rec_page.min_fiber_len.get()}")
    fiber_min_layers_theshold = rec_page.min_fiber_len.get()

    if output_stat_folder == "" or fiber_min_layers_theshold == "":
        showinfo(
            title='Information',
            message="Please fill all blanks to make the program run"
        )

    else:
        showinfo(
            title='Information',
            message=f"The analysis process is about to begin. Once completed, you can find the results "
                    f"in the designated output folder: {output_stat_folder}. Please note that the analysis "
                    f"can take a significant amount of time, depending on the size of the dataset and "
                    f"the capabilities of your machine."

            )
        try:
            root.destroy()
            run_through_gui_recalculation(input_obj_path, output_stat_folder,
                                         int(fiber_min_layers_theshold))
        except Exception as e:
            a = 1

class AnalysisTypePage:
    def __init__(self, analysis_root):
        self.selected_type = tk.StringVar()
        analysis_types = (('Run new analysis', 'new'),
                        ('Recalculate statistics', 'recalculate'))

        # label
        label = ctk.CTkLabel(master=analysis_root, text="Choose Type of Analysis: ")
        label.pack(fill='x', padx=10, pady=15)

        # Analysis type radio buttons
        for type in analysis_types:
            r = ctk.CTkRadioButton(
                analysis_root,
                text=type[0],
                value=type[1],
                variable=self.selected_type,

            )
            r.pack(fill='x', padx=50, pady=5)


        # Continue button
        self.continue_button = ctk.CTkButton(
            analysis_root,
            text="Continue")

        self.continue_button.pack(ipadx=5, ipady=5, expand=True, padx=0, pady=30)


class NewAnalysisSettingsPage:
    def __init__(self, data_root):
        self.input_folder = tk.StringVar()
        self.output_data_folder = tk.StringVar()
        self.output_obj_folder = tk.StringVar()
        self.nuc_channel = tk.StringVar()
        self.actin_channel = tk.StringVar()
        self.min_fiber_len = tk.DoubleVar()

        # configure the grid for data root label
        data_root.columnconfigure(0, weight=1)
        data_root.columnconfigure(1, weight=2)

        # input folder path
        ctk.CTkLabel(master=data_root, text='Images path:', anchor='w').grid(column=0, row=1, sticky=tk.W, padx=15, pady=15)
        input_parent = ctk.CTkLabel(master=data_root)
        input_parent.grid(column=1, row=1, sticky=tk.EW)
        input_parent.columnconfigure(0, weight=10)
        input_parent.columnconfigure(1, weight=1)

        input_textbox = ctk.CTkEntry(master=input_parent, textvariable=self.input_folder)
        input_textbox.grid(row=0, column=0, sticky=tk.EW)
        button_br_1 = ctk.CTkButton(master=input_parent, text="...", command=self.input_button, width=30)
        button_br_1.grid(row=0, column=1, padx=5)

        # output stat folder path
        ctk.CTkLabel(master=data_root, text='Save statistics at:', anchor='w').grid(column=0, row=2, sticky=tk.W, padx=15, pady=15)
        input_parent = ctk.CTkLabel(master=data_root)
        input_parent.grid(column=1, row=2, sticky=tk.EW)
        input_parent.columnconfigure(0, weight=10)
        input_parent.columnconfigure(1, weight=1)

        input_textbox = ctk.CTkEntry(master=input_parent, textvariable=self.output_data_folder)
        input_textbox.grid(row=0, column=0, sticky=tk.EW)
        button_br_2 = ctk.CTkButton(master=input_parent, text="...", command=self.output_data_button, width=30)
        button_br_2.grid(row=0, column=1, padx=5)

        # output objects folder path
        ctk.CTkLabel(master=data_root, text='Save analysed objects at*:', anchor='w').grid(column=0, row=3, sticky=tk.W, padx=15, pady=15)
        input_parent = ctk.CTkLabel(master=data_root)
        input_parent.grid(column=1, row=3, sticky=tk.EW)
        input_parent.columnconfigure(0, weight=10)
        input_parent.columnconfigure(1, weight=1)

        input_textbox = ctk.CTkEntry(master=input_parent, textvariable=self.output_obj_folder)
        input_textbox.grid(row=0, column=0, sticky=tk.EW)
        button_br_3 = ctk.CTkButton(master=input_parent, text="...", command=self.output_img_button, width=30)
        button_br_3.grid(row=0, column=1, padx=5)


        # Nucleus channel number
        ctk.CTkLabel(master=data_root, text='Nucleus channel', anchor='w').grid(column=0, row=4, sticky=tk.W, padx=15,
                                                                                pady=15)
        mask_combobox_nuc = ctk.CTkComboBox(master=data_root, values=["0", "1", "2", "3"], variable=self.nuc_channel)
        mask_combobox_nuc['values'] = ["0", "1", "2", "3"]

        # prevent typing a value
        mask_combobox_nuc['state'] = 'readonly'
        mask_combobox_nuc.grid(column=1, row=4, sticky=tk.W, padx=0, pady=15)
        mask_combobox_nuc.get()

        # Actin channel number
        ctk.CTkLabel(master=data_root, text='Actin channel', anchor='w').grid(column=0, row=5, sticky=tk.W, padx=15,
                                                                                pady=15)
        mask_combobox_actin = ctk.CTkComboBox(master=data_root, values=["0", "1", "2", "3"], variable=self.actin_channel)
        mask_combobox_actin['values'] = ["0", "1", "2", "3"]

        # prevent typing a value
        mask_combobox_actin['state'] = 'readonly'
        mask_combobox_actin.grid(column=1, row=5, sticky=tk.W, padx=0, pady=15)
        mask_combobox_actin.get()

        #minimal fiber length
        ctk.CTkLabel(master=data_root, text='Minimal fiber length (pixels)', anchor='w').grid(column=0, row=7, sticky=tk.W, padx=15, pady=15)
        spin_box = ctk.CTkEntry(master=data_root, width=50, textvariable=self.min_fiber_len)
        spin_box.grid(column=1, row=7, sticky=tk.W, padx=0, pady=15)

        self.analize_button = ctk.CTkButton(
            master=data_root,
            text="Analyze")
        self.analize_button.grid(column=1, row=9, sticky=tk.E, padx=15, pady=30)


    def input_button(self):
        self.input_folder.set(filedialog.askdirectory())

    def output_data_button(self):
        self.output_data_folder.set(filedialog.askdirectory())

    def output_img_button(self):
        self.output_obj_folder.set(filedialog.askdirectory())


class RecalculateAnalysisSettingsPage:
    def __init__(self, data_root):
        self.input_obj_folder = tk.StringVar()
        self.output_data_folder = tk.StringVar()
        self.min_fiber_len = tk.DoubleVar()

        # configure the grid for data root label
        data_root.columnconfigure(0, weight=1)
        data_root.columnconfigure(1, weight=2)

        # pre-analysed objects path
        ctk.CTkLabel(master=data_root, text='Analysed objects path:', anchor='w').grid(column=0, row=1, sticky=tk.W, padx=15,
                                                                             pady=15)
        input_parent = ctk.CTkLabel(master=data_root)
        input_parent.grid(column=1, row=1, sticky=tk.EW)
        input_parent.columnconfigure(0, weight=10)
        input_parent.columnconfigure(1, weight=1)

        input_textbox = ctk.CTkEntry(master=input_parent, textvariable=self.input_obj_folder)
        input_textbox.grid(row=0, column=0, sticky=tk.EW)
        button_br_1 = ctk.CTkButton(master=input_parent, text="...", command=self.input_button, width=30)
        button_br_1.grid(row=0, column=1, padx=5)

        # output stat folder path
        ctk.CTkLabel(master=data_root, text='Save statistics at:', anchor='w').grid(column=0, row=2, sticky=tk.W,
                                                                                    padx=15, pady=15)
        input_parent = ctk.CTkLabel(master=data_root)
        input_parent.grid(column=1, row=2, sticky=tk.EW)
        input_parent.columnconfigure(0, weight=10)
        input_parent.columnconfigure(1, weight=1)

        input_textbox = ctk.CTkEntry(master=input_parent, textvariable=self.output_data_folder)
        input_textbox.grid(row=0, column=0, sticky=tk.EW)
        button_br_2 = ctk.CTkButton(master=input_parent, text="...", command=self.output_data_button, width=30)
        button_br_2.grid(row=0, column=1, padx=5)

        # minimal nucleus area
        ctk.CTkLabel(master=data_root, text='Minimal fiber length (pixels)', anchor='w').grid(column=0, row=7,
                                                                                              sticky=tk.W, padx=15,
                                                                                              pady=15)
        spin_box = ctk.CTkEntry(master=data_root, width=50, textvariable=self.min_fiber_len)
        spin_box.grid(column=1, row=7, sticky=tk.W, padx=0, pady=15)

        self.analize_button = ctk.CTkButton(
            master=data_root,
            text="Analyze")
        self.analize_button.grid(column=1, row=9, sticky=tk.E, padx=15, pady=30)

    def input_button(self):
        self.input_obj_folder.set(filedialog.askdirectory())

    def output_data_button(self):
        self.output_data_folder.set(filedialog.askdirectory())


# design customization
ctk.set_appearance_mode("System")  # Modes: system (default), light, dark
ctk.set_default_color_theme("green")  # Themes: blue (default), dark-blue, green

root = ctk.CTk()
root.title('Afilament')
root.iconbitmap(r'..\docs\imgs\favicon.ico')
window_width = 350
window_height = 350

# get the screen dimension
screen_width = root.winfo_screenwidth()
screen_height = root.winfo_screenheight()

# find the center point
center_x = int(screen_width/2 - window_width / 2)
center_y = int(screen_height/2 - window_height / 2)

root.geometry(f'{window_width}x{window_height}+{center_x}+{center_y}')

# First page: Analysis type
analysis_frame = ctk.CTkFrame(master=root)
new_analysis_frame = ctk.CTkFrame(master=root)
recalculate_frame = ctk.CTkFrame(master=root)

for frame in (analysis_frame, new_analysis_frame, recalculate_frame):
    frame.grid(row=0, column=0, sticky='nsew')

analysis_frame.tkraise()

analysis_page = AnalysisTypePage(analysis_frame)
analysis_page.continue_button.configure(command=show_data_page)

data_page = NewAnalysisSettingsPage(new_analysis_frame)
data_page.analize_button.configure(command=run_new_analysis)

rec_page = RecalculateAnalysisSettingsPage(recalculate_frame)
rec_page.analize_button.configure(command=run_recalculation)

root.mainloop()


