####################################################################################################
#                                             main.py                                              #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Created: 07/02/25                                                                                #
#                                                                                                  #
# Purpose: Main script for the BasisREMY tool.                                                     #
#                                                                                                  #
####################################################################################################


#*************#
#   imports   #
#*************#
import matplotlib.pyplot as plt
import numpy as np
import os
import pathlib
import plotly.graph_objs as go
import plotly.io as pio
import threading
import tkinter as tk

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

from multiprocessing import Pool

from oct2py import Oct2Py

from tkinter import ttk
from tkinter import filedialog
from tkinterdnd2 import TkinterDnD, DND_FILES
from tkinterweb import HtmlFrame

from PIL import Image, ImageTk

# own
from remy.MRSinMRS import DataReaders, Table, setup_log, write_log


#**************************************************************************************************#
#                                           Application                                            #
#**************************************************************************************************#
#                                                                                                  #
# The GUI application for the BasisREMY tool. Each tab is a different step in the process,         #
# starting with the data selection and REMY extraction, continuing with the parameter              #
# configuration, and ending with the basis set simulation.                                         #
#                                                                                                  #
#**************************************************************************************************#
class Application(TkinterDnD.Tk):

    def __init__(self):
        super().__init__()

        # initialize data backend
        self.BasisREMY = BasisREMY()

        # define a fixed color palette
        self.main_color = "#607389"   # primary buttons and highlights
        self.bg_1_color = "#f0f0f0"   # window and frame background
        self.bg_2_color = "#e0e0e0"   # secondary panels
        self.bg_3_color = "#f0f0f0"   # other levels
        self.text_color = "#000000"   # always black text

        # set overall widget foreground/background
        self.option_add("*Foreground", self.text_color)
        self.option_add("*Background", self.bg_1_color)
        self.option_add("*Entry.InsertBackground", self.text_color)   # for entry cursor
        self.option_add("*Checkbutton.SelectColor", self.bg_3_color)   # for checkbutton

        # configure ttk theme to a fixed one
        style = ttk.Style(self)
        style.theme_use('clam')  # a neutral, simple theme

        # override ttk widget styling
        style.configure('.', background=self.bg_1_color, foreground=self.text_color)
        style.configure('TFrame', background=self.bg_1_color)
        style.configure('TLabel', background=self.bg_1_color, foreground=self.text_color)
        style.configure('TButton', background=self.main_color, foreground=self.text_color)
        style.configure('TCheckbutton', background=self.bg_1_color, foreground=self.text_color)
        style.configure('TEntry', fieldbackground=self.bg_2_color, foreground=self.text_color)
        style.configure('TCombobox', fieldbackground=self.bg_2_color, foreground=self.text_color)
        style.configure('TNotebook', background=self.bg_1_color)
        style.configure('TNotebook.Tab', font=("Arial", 12, "normal"), background=self.bg_2_color,
                        foreground=self.text_color)

        # setup window
        self.title("BasisREMY")
        self.geometry("1000x700")
        self.configure(bg=self.bg_1_color)

        # create UI elements
        self.create_widgets()

    def create_widgets(self):
        # load and display the header image
        header_image = Image.open("imgs/basisremy_header.png")
        header_image = header_image.resize((1000, 50))
        header_photo = ImageTk.PhotoImage(header_image)
        header_label = tk.Label(self, image=header_photo, bg="#f0f0f0")
        header_label.image = header_photo  # keep a reference to prevent garbage collection
        header_label.pack(fill=tk.X)

        # header text
        header_label = tk.Label(self, text="BasisREMY", font=("Arial bold", 20), bg="#f0f0f0")
        header_label.pack(pady=10)
        header_label = tk.Label(self, text="A tool for generating study-specific basis "
                                           "sets directly from raw MRS data.",
                                font=("Arial", 16), bg=self.bg_1_color)
        header_label.pack(pady=1)

        # notebook (tabbed interface)
        self.notebook = ttk.Notebook(self)
        self.notebook.pack(padx=10, pady=10, expand=True, fill=tk.BOTH)

        # tab frames
        self.tab1 = ttk.Frame(self.notebook)
        self.tab2 = ttk.Frame(self.notebook)
        self.tab3 = ttk.Frame(self.notebook)
        self.notebook.add(self.tab1, text="Data Selection", state="normal")
        self.notebook.add(self.tab2, text="Parameter Configuration", state="disabled")
        self.notebook.add(self.tab3, text="Basis Simulation", state="disabled")

        self.tab1_widgets()
        self.tab2_widgets()
        self.tab3_widgets()

        # apply custom styles to tabs
        style = ttk.Style()
        style.configure("TNotebook.Tab",
                        font=("Arial", 12))  # set the font size and style for the tabs

        self.update()
        self.geometry(f"{self.winfo_reqwidth()}x{self.winfo_reqheight()}")

    def tab1_widgets(self):
        # drag and drop area
        drag_area = tk.Label(self.tab1, text="Select File!\n\nDrag and Drop Files Here...",
                             relief="solid", width=50, height=12, bg=self.bg_2_color,
                             font=("Arial bold", 12))
        drag_area.pack(padx=1, pady=1, expand=True)
        drag_area.drop_target_register(DND_FILES)
        drag_area.dnd_bind('<<Drop>>', self.on_drop)
        drag_area.bind("<Button-1>", self.select_file)  # bind left-click to select_file

        # file label
        self.file_label = tk.Label(self.tab1, text="No file selected.", bg=self.bg_1_color,
                                   font=("Arial italic", 12))
        self.file_label.pack(pady=5)

        # process button
        self.process_button = tk.Button(self.tab1, text="Process File", command=self.process_file,
                                        state=tk.DISABLED, bg=self.main_color, fg=self.text_color,
                                        font=("Arial", 12, "bold"))
        self.process_button.pack(pady=5)

        # skip button
        self.skip_button = tk.Button(self.tab1, text="Skip", command=self.skip_file,
                                     state=tk.NORMAL, bg=self.main_color, fg=self.text_color,
                                     font=("Arial", 12, "bold"))
        self.skip_button.pack(pady=5)

    def tab2_widgets(self):
        # clear existing widgets in tab2 (and 3 for returns)
        self.reset_tab(self.tab2)

        # backend toggle
        backend_frame = tk.Frame(self.tab2)
        backend_frame.pack(fill=tk.X, padx=10, pady=(10, 0))

        tk.Label(backend_frame, text="Select Backend:",
                 font=("Arial", 12, "bold")).pack(side=tk.LEFT, padx=(0, 10))

        self.backend_var = tk.StringVar(value=self.BasisREMY.backend.name)
        backend_options = self.BasisREMY.available_backends

        backend_combo = ttk.Combobox(backend_frame, textvariable=self.backend_var,
                                     values=backend_options, font=("Arial", 12), state="readonly")
        backend_combo.pack(side=tk.LEFT)

        def switch_backend(event=None):
            selected_backend = self.backend_var.get()
            if selected_backend != self.BasisREMY.backend.name:
                self.BasisREMY.set_backend(selected_backend)  # you must define this method
                self.tab2_widgets()  # redraw tab2 widgets for the new backend

        backend_combo.bind("<<ComboboxSelected>>", switch_backend)

        # create a container frame to hold both parameter and metabolite sections
        container = tk.Frame(self.tab2)
        container.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # create the parameters frame (left side)
        params_frame = tk.Frame(container)
        params_frame.grid(row=0, column=0, sticky="nsew", padx=(0, 10))

        # create the metabolites frame (right side)
        metabs_frame = tk.Frame(container)
        metabs_frame.grid(row=0, column=1, sticky="nsew")

        # configure grid weights for responsive resizing
        container.columnconfigure(0, weight=1)
        container.columnconfigure(1, weight=1)
        container.rowconfigure(0, weight=1)

        # function to update mandatory_params when an Entry or Combobox changes
        def update_param(key, var):
            self.BasisREMY.backend.mandatory_params[key] = var.get()
            self.validate_inputs()

        # populate the parameters frame
        row = 0
        for key, value in self.BasisREMY.backend.mandatory_params.items():
            if key in self.BasisREMY.backend.file_selection:
                # label for the path (file)
                label = tk.Label(params_frame, text=f"{key}:", font=("Arial", 12, "bold"))
                label.grid(row=row, column=0, padx=0, pady=5, sticky="e")

                # stringVar for the selected path
                self.file_var = tk.StringVar(value="missing input")

                entry = tk.Entry(params_frame, textvariable=self.file_var, font=("Arial", 12))
                entry.grid(row=row, column=1, padx=0, pady=5, sticky="ew")

                def select_path():
                    # select a file (you could also switch to askdirectory if needed per key)
                    file_path = filedialog.askopenfilename()
                    if file_path:
                        self.file_var.set(file_path)
                        self.BasisREMY.backend.mandatory_params[key] = file_path
                    else:
                        self.file_var.set("missing input")
                        self.BasisREMY.backend.mandatory_params[key] = None
                    self.validate_inputs()

                button = tk.Button(params_frame, text="Browse", command=select_path)
                button.grid(row=row, column=2, padx=0, pady=5)

                # automatically update backend dict on any entry edit
                self.file_var.trace_add('write', lambda *args, key=key,
                                                        var=self.file_var: update_param(key, var))

            elif key == 'Output Path':
                # label for output directory
                label = tk.Label(params_frame, text=f"{key}:", font=("Arial", 12, "bold"))
                label.grid(row=row, column=0, padx=0, pady=5, sticky="e")

                # StringVar for the directory path
                self.out_path_var = tk.StringVar(value="missing input")

                entry = tk.Entry(params_frame, textvariable=self.out_path_var, font=("Arial", 12))
                entry.grid(row=row, column=1, padx=0, pady=5, sticky="ew")

                # button to open directory dialog
                def select_directory():
                    directory = filedialog.askdirectory()
                    if directory:
                        self.out_path_var.set(directory)
                        self.BasisREMY.backend.mandatory_params['Output Path'] = directory
                    else:
                        self.out_path_var.set("missing input")
                        self.BasisREMY.backend.mandatory_params['Output Path'] = None
                    self.validate_inputs()

                button = tk.Button(params_frame, text="Browse", command=select_directory)
                button.grid(row=row, column=2, padx=0, pady=5)

                # trace changes to the StringVar
                self.out_path_var.trace_add('write', lambda *args, key=key,
                                                            var=self.out_path_var: update_param(key, var))

            elif key == 'Metabolites':
                # populate the metabolites frame
                metabs_label = tk.Label(metabs_frame, text="Select Metabolites:", font=("Arial", 12, "bold"))
                metabs_label.grid(row=0, column=0, columnspan=2, pady=5)

                self.metab_vars = {}
                row = 1
                col = 0
                for metab in self.BasisREMY.backend.metabs:
                    # create checkbox for each metabolite
                    var = tk.BooleanVar(value=metab in self.BasisREMY.backend.mandatory_params.get('Metabolites', []))
                    checkbutton = tk.Checkbutton(metabs_frame, text=metab, variable=var)
                    checkbutton.grid(row=row, column=col, sticky="w", padx=5, pady=2)
                    self.metab_vars[metab] = var

                    # arrange checkbuttons in a grid with 5 columns
                    col += 1
                    if col == 5:
                        col = 0
                        row += 1

                # update mandatory_params with selected metabolites
                def update_metabs(*args):
                    selected_metabs = [metab for metab, var in self.metab_vars.items() if var.get()]
                    self.BasisREMY.backend.mandatory_params['Metabolites'] = selected_metabs
                    self.validate_inputs()

                # trace variable changes to update the list
                for var in self.metab_vars.values():
                    var.trace_add('write', update_metabs)

            elif key in self.BasisREMY.backend.dropdown:
                # label for dropdown parameters
                label = tk.Label(params_frame, text=f"{key}:", font=("Arial", 12, "bold"))
                label.grid(row=row, column=0, padx=0, pady=5, sticky="e")

                # StringVar for the Combobox
                var = tk.StringVar(value=str(value) if value is not None else "missing input")

                # Combobox for dropdown parameters
                combobox = ttk.Combobox(params_frame, textvariable=var, values=self.BasisREMY.backend.dropdown[key],
                                        font=("Arial", 12))
                combobox.grid(row=row, column=1, padx=0, pady=5, sticky="ew")

                # trace changes to the StringVar
                var.trace_add('write', lambda *args, key=key, var=var: update_param(key, var))

            else:
                # label for other parameters
                label = tk.Label(params_frame, text=f"{key}:", font=("Arial", 12, "bold"))
                label.grid(row=row, column=0, padx=0, pady=5, sticky="e")

                # StringVar for the Entry
                var = tk.StringVar(value=str(value) if value is not None else "missing input")

                # entry for other parameters
                entry = tk.Entry(params_frame, textvariable=var, font=("Arial", 12))
                entry.grid(row=row, column=1, padx=0, pady=5, sticky="ew")

                # trace changes to the StringVar
                var.trace_add('write', lambda *args, key=key, var=var: update_param(key, var))

            row += 1

        # the Simulate Basis Set button
        self.simulate_button = tk.Button(
            self.tab2,
            text="Simulate Basis Set",
            command=self.simulate_basis,
            bg=self.main_color,
            fg=self.text_color,
            font=("Arial", 12, "bold"),
            state=tk.DISABLED,
        )
        self.simulate_button.pack(pady=5)

        # the back button
        self.back_button = tk.Button(
            self.tab2,
            text="Back",
            command=lambda: self.notebook.select(0),
            bg=self.main_color,
            fg=self.text_color,
            font=("Arial", 12, "bold"),
        )
        self.back_button.pack(pady=5)

    def tab3_widgets(self):
        # clear existing widgets in tab3.
        self.reset_tab(self.tab3)

        # create the simulation progress area
        self.simulation_status_label = tk.Label(self.tab3, text="Simulating basis set...", font=("Arial", 16, "bold"))
        self.simulation_status_label.pack(pady=10)

        self.progress_frame = tk.Frame(self.tab3)
        self.progress_frame.pack(fill=tk.X, padx=10, pady=10)

        self.progress = ttk.Progressbar(self.progress_frame, orient="horizontal", length=400, mode="determinate")
        self.progress.pack(pady=10)

        self.progress_label = tk.Label(self.progress_frame, text="0%", font=("Arial", 12, "bold"))
        self.progress_label.pack(pady=10)

        # container for the plot and checkboxes (populated after simulation)
        self.plot_container = tk.Frame(self.tab3)
        self.plot_container.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)

        # the back button remains visible.
        self.back_button = tk.Button(
            self.tab3,
            text="Back",
            command=lambda: self.notebook.select(1),
            bg=self.main_color,
            fg=self.text_color,
            font=("Arial", 12, "bold")
        )
        self.back_button.pack(pady=5)

    def reset_tab(self, tab):
        # clear existing widgets in the specified tab
        for widget in tab.winfo_children(): widget.destroy()

    def on_drop(self, event):
        # files are dragged and dropped into the drag area
        file_path = event.data
        self.file_label.config(text=f"Selected File: {file_path}")
        self.process_button.config(state=tk.NORMAL)   # enable the process button

    def select_file(self, event=None):
        # files are selected using the file dialog
        file_path = filedialog.askopenfilename(
            title="Select MRS Data File",
            filetypes=[
                ("MRS Data Files", (
                    "*.7", "*.ima", "*.rda", "*.dat",
                    "*.spar", "*.method", "*.nii"
                ))
            ]
        )
        if file_path:
            self.file_label.config(text=f"Selected File: {file_path}")
            self.process_button.config(state=tk.NORMAL)   # enable the process button
        else:
            self.file_label.config(text="No file selected.")
            self.process_button.config(state=tk.DISABLED)

    def process_file(self):
        # process button is clicked -> extract REMY and skip to config tab
        file_path = self.file_label.cget("text").replace("Selected File: ", "")
        if file_path:
            print(f"Processing file: {file_path}")

            # run REMY on the selected file
            MRSinMRS = self.BasisREMY.runREMY(file_path)
            params, opt = self.BasisREMY.parseREMY(MRSinMRS)

            # update the mandatory parameters
            self.BasisREMY.backend.mandatory_params.update(params)
            self.BasisREMY.backend.optional_params.update(opt)

            # update the parameter configuration tab
            self.tab2_widgets()

            # move to the next tab
            self.notebook.tab(1, state="normal")
            self.notebook.select(1)

        else:
            print("No file selected.")

    def skip_file(self):
        # skip to the next tab without processing the file
        self.notebook.tab(1, state="normal")
        self.notebook.select(1)

    def validate_inputs(self):
        # check if all mandatory parameters are filled
        all_params_filled = all(
            self.BasisREMY.backend.mandatory_params[key] not in (None, "", "missing input")
            for key in self.BasisREMY.backend.mandatory_params
            if key != 'Metabolites'
        )

        # check if at least one metabolite is selected
        at_least_one_metab_selected = any(var.get() for var in self.metab_vars.values())

        # enable the simulate button if both conditions are met
        if all_params_filled and at_least_one_metab_selected:
            self.simulate_button.config(state=tk.NORMAL)
        else:
            self.simulate_button.config(state=tk.DISABLED)

    def simulate_basis(self):
        # move to the next tab and simulate the basis set
        self.notebook.tab(2, state="normal")
        self.notebook.select(2)

        print("Simulating basis set with the following parameters:")
        for key, value in self.BasisREMY.backend.mandatory_params.items():
            print(f"{key}: {value}")

        # initialize the progress bar
        self.progress["value"] = 0
        self.progress["maximum"] = len(self.BasisREMY.backend.mandatory_params['Metabolites'])

        # run the simulation in a separate thread to keep the GUI responsive
        threading.Thread(target=self.run_simulation_with_progress, args=(self.on_simulation_complete,)).start()

    def run_simulation_with_progress(self, callback):
        def progress_callback(step, total_steps):
            # update the progress bar
            self.progress["value"] = step
            self.progress_label.config(text=f"{step * 100 // total_steps}%")

            # update the GUI
            self.update_idletasks()

        # run the simulation with the progress callback
        basis = self.BasisREMY.backend.run_simulation_with_progress(self.BasisREMY.backend.mandatory_params, progress_callback)
        self.after(0, callback, basis)

    def on_simulation_complete(self, basis):
        print("Simulation complete.")
        self.basis_set = basis

        # remove simulation progress widgets.
        self.simulation_status_label.pack_forget()
        self.progress_frame.pack_forget()

        # create a container frame for both the plot and the checkboxes.
        plot_frame = tk.Frame(self.plot_container)
        plot_frame.pack(fill=tk.BOTH, expand=True)

        # left frame for the canvas.
        canvas_frame = tk.Frame(plot_frame)
        canvas_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        # right frame for the merged legend (checkboxes with color patches).
        checkbox_frame = tk.Frame(plot_frame)
        checkbox_frame.pack(side=tk.RIGHT, fill=tk.Y, padx=5, pady=1)

        # create a smaller matplotlib figure.
        self.figure = Figure(figsize=(6, 2), dpi=100)
        self.figure.patch.set_alpha(0.0)  # set figure background to be transparent
        self.ax = self.figure.add_subplot(111)
        self.ax.set_facecolor(self.bg_1_color)  # set axes background to be the same as the GUI
        # self.ax.set_title("Basis Simulation Results")
        self.ax.set_xlabel("Chemical Shift [ppm]")

        # no y-axis ticks or labels.
        self.ax.set_yticks([])
        self.ax.set_yticklabels([])

        # embed the Matplotlib canvas in the left frame.
        self.canvas = FigureCanvasTkAgg(self.figure, master=canvas_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.canvas.get_tk_widget().configure(bg=self.bg_1_color)

        # create the toolbar but overlay it on top of the canvas.
        self.toolbar = NavigationToolbar2Tk(self.canvas, canvas_frame, pack_toolbar=False)
        self.toolbar.update()
        self.toolbar.place(in_=self.canvas.get_tk_widget(), relx=0, rely=0)

        # create checkboxes merged with a color patch (serving as the legend).
        self.checkbox_vars = {}
        self.metab_colors = {}
        # get the default color cycle from Matplotlib.
        default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        for i, metab in enumerate(self.basis_set.keys()):
            # assign a color to this metabolite.
            color = default_colors[i % len(default_colors)]
            self.metab_colors[metab] = color

            var = tk.BooleanVar(value=True)
            self.checkbox_vars[metab] = var

            # create a frame row for the color patch and the checkbutton.
            row_frame = tk.Frame(checkbox_frame)
            row_frame.pack(anchor='w')

            # create a small label as a color patch.
            color_patch = tk.Label(row_frame, bg=color, width=2, height=1, font=("Arial", 3))
            color_patch.pack(side=tk.LEFT, padx=5)

            # create the checkbutton for the metabolite.
            checkbutton = tk.Checkbutton(
                row_frame,
                text=metab,
                variable=var,
                command=self.update_plot
            )
            checkbutton.pack(side=tk.LEFT)

        # render the initial plot.
        self.update_plot()

    def update_plot(self):
        # clear the axis and reapply basic settings.
        self.ax.clear()
        # self.ax.set_title("Basis Simulation Results")
        self.ax.set_facecolor(self.bg_1_color)
        self.ax.set_xlabel("Chemical Shift [ppm]")

        # no y-axis ticks or labels.
        self.ax.set_yticks([])
        self.ax.set_yticklabels([])

        # compute ppm axis using cf
        cf = self.BasisREMY.backend.mandatory_params['Center Freq']
        bw = self.BasisREMY.backend.mandatory_params['Bandwidth']
        points = self.BasisREMY.backend.mandatory_params['Samples']
        ppm_axis = 1e6 * np.linspace(-bw/2 / cf, bw/2 / cf, points)
        ppm_axis = np.flip(ppm_axis) + 4.65   # TODO: make this more general

        # plot each metabolite's data if its checkbox is selected.
        for metab, var in self.checkbox_vars.items():
            if var.get() and metab in self.basis_set:
                data = self.basis_set[metab]
                ydata = np.real(np.fft.fftshift(np.fft.fft(data)))
                self.ax.plot(ppm_axis, ydata, color=self.metab_colors[metab])

        # limit the x-axis to 0 - 10 ppm (TODO: this is very proton specific)
        self.ax.set_xlim(0, 10)
        self.ax.invert_xaxis()

        self.canvas.draw()



#**************************************************************************************************#
#                                            BasisREMY                                             #
#**************************************************************************************************#
#                                                                                                  #
# The BasisREMY class is the main class for the BasisREMY tool. It provides the functionality to   #
# extract REMY parameters from MRS data and simulate a basis.                                      #
#                                                                                                  #
#**************************************************************************************************#
class BasisREMY:
    def __init__(self, backend='sLaserSim'):
        self.DRead = DataReaders()
        self.Table = Table()

        self.backends = {
            'LCModel': LCModelBackend(),
            'sLaserSim': sLaserBackend(),
        }
        self.backend = self.backends[backend]
        self.available_backends = list(self.backends.keys())

    def set_backend(self, backend):
        # set the backend to the selected one
        if backend in self.backends:
            old_backend = self.backend
            self.backend = self.backends[backend]
            self.backend.update_from_backend(old_backend)  # update parameters
            print(f"Backend set to: {self.backend.name}")
        else:
            raise ValueError(f"Unknown backend: {backend}. Available backends: {self.available_backends}")

    def run(self, import_fpath, export_fpath=None, method=None, userParams={}, optionalParams={}):
        # run REMY on the selected file
        MRSinMRS = self.runREMY(import_fpath, method)
        params, opt = self.parseREMY(MRSinMRS)
        params['Output Path'] = export_fpath if export_fpath is not None else './'

        # update the mandatory parameters
        self.backend.mandatory_params.update(params)
        self.backend.mandatory_params.update(userParams)

        # update the optional parameters
        self.backend.optional_params.update(opt)
        self.backend.optional_params.update(optionalParams)

        # run fidA simulation
        basis = self.backend.run_simulation(self.mandatory_params)

        # plot the basis set
        import matplotlib.pyplot as plt
        plt.figure()
        for key, value in basis.items():
            plt.plot(np.fft.fft(value), label=key)
        plt.legend()
        plt.show()

    def runREMY(self, import_fpath, method=None):
        # run REMY datareader on the selected file
        if method is None: suf = pathlib.Path(import_fpath).suffix.lower()
        else: suf = method
        log = None
        if suf == '.dat':   # Siemens Twix file
            write_log(log, 'Data Read: Siemens Twix uses pyMapVBVD ')  # log - pyMapVBVD
            MRSinMRS, log = self.DRead.siemens_twix(import_fpath, log)
            vendor_selection = 'Siemens'
        elif suf == '.ima':  # Siemens Dicom file
            write_log(log, 'Data Read: Siemens Dicom uses pydicom ')  # log - pyDicom
            MRSinMRS, log = self.DRead.siemens_ima(import_fpath, log)
            vendor_selection = 'Siemens'
        elif suf == '.rda':  # Siemens RDA file
            write_log(log, 'Data Read: Siemens RDA directly read with RMY ')  # log - pyDicom
            MRSinMRS, log = self.DRead.siemens_rda(import_fpath,    log)
            vendor_selection = 'Siemens'
        elif suf == '.spar':  # Philips SPAR file
            write_log(log, 'Data Read: Philips SPAR uses spec2nii ')  # log - spec2nii
            MRSinMRS, log = self.DRead.philips_spar(import_fpath, log)
            vendor_selection = 'Philips'
        elif suf == '.7':  # GE Pfile
            write_log(log, 'Data Read: GE Pfile uses spec2nii ')  # log - spec2nii
            MRSinMRS, log = self.DRead.ge_7(import_fpath, log)
            vendor_selection = 'GE'
        elif suf == 'bruker_method':  # Bruker Method file
            write_log(log, 'Data Read: Bruker Method uses spec2nii ')  # log - spec2nii
            MRSinMRS, log = self.DRead.bruker_method(import_fpath, log)
            vendor_selection = 'Bruker'
        elif suf == 'bruker_2dseq':  # Bruker 2dseq file
            write_log(log, 'Data Read: Bruker uses BrukerAPI ' +  # log - BrukerAPI
                      'developed by Tomáš Pšorn\n\t' +
                      'github.com/isi-nmr/brukerapi-python')
            MRSinMRS, log = self.DRead.bruker_2dseq(import_fpath, log)
            vendor_selection = 'Bruker'
        elif suf == '.nii' or suf == '.nii.gz':
            write_log(log, 'Data Read: NIfTI json side car')  # log - NIfTI JSON side car
            MRSinMRS, log = self.DRead.nifti_json(import_fpath, log)
            vendor_selection = 'NIfTI'
        else:
            raise ValueError(f'Unknown file format {suf}! Valid formats are:'
                             f' .dat, .ima, .rda, .spar, .7, bruker_method, bruker_2dseq, .nii, .nii.gz')

        # clean the data
        dtype_selection = suf.replace('bruker_', '').replace('.', '')
        MRSinMRS = self.Table.table_clean(vendor_selection, dtype_selection, MRSinMRS)
        return MRSinMRS

    def parseREMY(self, MRSinMRS):
        # extract as much information as possible from the MRSinMRS dict
        mandatory = {
            'Sequence': None,  # TODO: find a way to get this from REMY
            'Samples': MRSinMRS.get('samples', None),
            'Bandwidth': MRSinMRS.get('sample_frequency', None),
            'Bfield': MRSinMRS.get('FieldStrength', None),
            'Linewidth': 1,   # TODO: find how to get from REMY
            'TE': MRSinMRS.get('EchoTime', None),
            'TE2': 0,   # attention! - only holds for SpinEcho or STEAM
                        # TODO: find sound solution!
            'Center Freq': MRSinMRS.get('synthesizer_frequency', None),
        }
        optional = {
            'Nucleus': MRSinMRS.get('Nucleus', None),
            'TR': MRSinMRS.get('RepetitionTime', None),
        }
        return mandatory, optional

    def run_gui(self):
        app = Application()
        app.mainloop()


def initialize_octave():
    # initialize an Octave session with needed paths
    octave = Oct2Py()
    octave.eval("warning('off', 'all');")
    octave.addpath('./fidA/inputOutput/')
    octave.addpath('./fidA/processingTools/')
    octave.addpath('./fidA/simulationTools/')
    return octave

def sim_lcmrawbasis_mp(n, sw, Bfield, lb, metab, tau1, tau2, addref, makeraw, seq, out_path):
    # run the simulation in a separate Octave session for multiprocessing
    octave = initialize_octave()
    result = octave.feval('sim_lcmrawbasis', n, sw, Bfield, lb, metab,
                          tau1, tau2, addref, makeraw, seq, out_path)
    octave.exit()  # ensure the Octave session is properly closed
    return metab, result[:, 0] + 1j * result[:, 1]



#**************************************************************************************************#
#                                             Backend                                              #
#**************************************************************************************************#
#                                                                                                  #
# Defines the backend structure for the simulations. Inherit from this class to create a new       #
# simulation backend with all mandatory attributes and methods.                                    #
#                                                                                                  #
#**************************************************************************************************#
class Backend:
    def __init__(self):
        self.name = None

        # define possible metabolites
        self.metabs = {}

        # dropdown options
        self.dropdown = {}

        # add file selection fields
        self.file_selection = []

        # define dictionary of mandatory parameters
        self.mandatory_params = {}

        # define dictionary of optional parameters
        self.optional_params = {}

    def update_from_backend(self, backend):
        # update the backend parameters from another backend instance
        # TODO: make sure some parameters are reset (see e.g. Sequence in dropdown)
        self.metabs.update({k: v for k, v in backend.metabs.items() if k in self.metabs})
        self.mandatory_params.update({k: v for k, v in backend.mandatory_params.items()
                                      if k in self.mandatory_params})
        self.optional_params.update({k: v for k, v in backend.optional_params.items()
                                     if k in self.optional_params})
        # note: for now no dropdown options are updated, as they are too specific to the backend
        #       at this point

    def run_simulation(self, params):
        raise NotImplementedError("This method should be overridden by subclasses.")

    def run_simulation_with_progress(self, params, progress_callback):
        raise NotImplementedError("This method should be overridden by subclasses.")



#**************************************************************************************************#
#                                          sLaserBackend                                           #
#**************************************************************************************************#
#                                                                                                  #
# Implements the basis set simulation backend for the sLASER sequence.                             #
#                                                                                                  #
#**************************************************************************************************#
class sLaserBackend(Backend):
    def __init__(self):
        super().__init__()
        self.name = 'sLaserSim'

        # init fidA
        self.octave = Oct2Py()
        self.octave.eval("warning('off', 'all');")
        self.octave.addpath('./fidA/inputOutput/')
        self.octave.addpath('./fidA/processingTools/')
        self.octave.addpath('./fidA/simulationTools/')

        self.octave.addpath('./jbss/dependencies/')
        self.octave.addpath('./jbss/')

        # define possible metabolites
        self.metabs = {
            'Ala': False,
            'Asc': True,
            'Asp': False,
            'Bet': False,
            'Ch': False,
            'Cit': False,
            'Cr': True,
            'EtOH': False,
            'GABA': True,
            'GABA_gov': False,
            'GABA_govind': False,
            'GPC': True,
            'GSH': True,
            'GSH_v2': False,
            'Glc': True,
            'Gln': True,
            'Glu': True,
            'Gly': True,
            'H2O': False,
            'Ins': True,
            'Lac': True,
            'NAA': True,
            'NAAG': True,
            'PCh': True,
            'PCr': True,
            'PE': True,
            'Phenyl': False,
            'Ref0ppm': False,
            'Scyllo': True,
            'Ser': False,
            'Tau': True,
            'Tau_govind': False,
            'Tyros': False,
            'bHB': False,
            'bHG': False,
        }

        # dropdown options
        self.dropdown = {
            'System': ['Philips', 'Siemens'],
            'Sequence': ['sLASER'],
        }

        # add file selection fields
        self.file_selection = ['Path to Pulse']

        # define dictionary of mandatory parameters
        self.mandatory_params = {
            "System": "Philips",
            "Sequence": "sLASER",
            "Basis Name": "test.basis",
            "B1max": 22.,
            "Flip Angle": 180.,
            "RefTp": 4.5008,   # duration of the refocusing pulse
            "Samples": None,
            "Bandwidth": None,
            "Linewidth": 1.,
            "Bfield": None,

            # TODO: see that REMY gets the right values
            "thkX": 2.,  # in cm
            "thkY": 2.,
            "fovX": 3.,  # in cm   (if not found, default to +1 slice thickness)
            "fovY": 3.,

            "nX": 64.,
            "nY": 64.,

            "TE": None,
            "Center Freq": None,
            "Metabolites": [key for key, value in self.metabs.items() if value],

            "Tau 1": 15.,   # fake timing
            "Tau 2": 13.,

            "Path to Pulse": None,
            "Output Path": None,
        }


        # define dictionary of optional parameters
        self.optional_params = {
            'Nucleus': None,
            'TR': None,
        }

    def run_simulation(self, params):
        pass

    def run_simulation_with_progress(self, params, progress_callback):
        # create the output directory if it does not exist
        if not os.path.exists(params['Output Path']): os.makedirs(params['Output Path'])

        # fixed parameters
        params.update({
            "Curfolder": os.getcwd() + '/jbss/',
            "Path to FIA-A": os.getcwd() + "/fidA/",
            "Path to Spin System": os.getcwd() + "/jbss/my_mets/",
            "Display": False,
        })

        def sLASER_makebasisset_function(curfolder,pathtofida,system,
                seq_name,basis_name,B1max,flip_angle,refTp,Npts,sw,lw,Bfield,
                thkX,thkY,fovX,fovY,nX,nY,te,centreFreq,spinSysList,tau1,tau2,
                path_to_pulse,path_to_save,path_to_spin_system,display):
            results = self.octave.feval('sLASER_makebasisset_function', curfolder,pathtofida,system,
                                        seq_name,basis_name,B1max,flip_angle,refTp,Npts,sw,lw,Bfield,
                                        thkX,thkY,fovX,fovY,nX,nY,te,centreFreq,spinSysList,tau1,tau2,
                                        path_to_pulse,path_to_save,path_to_spin_system,display)
            return metab, results[:, 0] + 1j * results[:, 1]

        tasks = [(params['Curfolder'], params['Path to FIA-A'], params['System'],
                  params['Sequence'], params['Basis Name'], params['B1max'],
                  params['Flip Angle'], params['RefTp'], params['Samples'],
                  params['Bandwidth'], params['Linewidth'], params['Bfield'],
                  params['thkX'], params['thkY'], params['fovX'], params['fovY'],
                  params['nX'], params['nY'], params['TE'],
                  params['Center Freq'], [metab], params['Tau 1'], params['Tau 2'],
                  params['Path to Pulse'], params['Output Path'],
                  params['Path to Spin System'], params['Display'])
                 for metab in params['Metabolites']]

        # initialize the progress bar
        total_steps = len(tasks)
        progress_step = 100 / total_steps

        # run simulations sequentially
        basis_set = {}
        for i, task in enumerate(tasks):
            metab, data = sLASER_makebasisset_function(*task)
            basis_set[metab] = data
            progress_callback(i + 1, total_steps)   # update the progress bar
        return basis_set



#**************************************************************************************************#
#                                          LCModelBackend                                          #
#**************************************************************************************************#
#                                                                                                  #
# Implements the basis set simulation backend using the FID-A sim_lcmrawbasis.m function. A very   #
# simplified simulation for SE, PRESS, STEAM, and LASER sequences.                                 #
#                                                                                                  #
#**************************************************************************************************#
class LCModelBackend(Backend):
    def __init__(self):
        super().__init__()
        self.name = 'LCModel'

        # init fidA
        self.octave = Oct2Py()
        self.octave.eval("warning('off', 'all');")
        self.octave.addpath('./fidA/inputOutput/')
        self.octave.addpath('./fidA/processingTools/')
        self.octave.addpath('./fidA/simulationTools/')

        # define possible metabolites
        self.metabs = {
            'Ala': False,
            'Asc': True,
            'Asp': False,
            'Ch': False,
            'Cit': False,
            'Cr': True,
            'EtOH': False,
            'GABA': True,
            'GPC': True,
            'GSH': True,
            'Glc': True,
            'Gln': True,
            'Glu': True,
            'Gly': True,
            'H2O': False,
            'Ins': True,
            'Lac': True,
            'Lip': False,
            'NAA': True,
            'NAAG': True,
            'PCh': True,
            'PCr': True,
            'PE': True,
            'Phenyl': False,
            'Ref0ppm': False,
            'Scyllo': True,
            'Ser': False,
            'Tau': True,
            'Tyros': False,
        }

        # dropdown options
        self.dropdown = {
            'Sequence': ['Spin Echo', 'PRESS', 'STEAM', 'LASER'],
            # 'se' for Spin Echo, 'p' for Press, 'st' for Steam, or 'l' for LASER
            'Add Ref.': ['Yes', 'No'],
            'Make .raw': ['Yes', 'No'],
        }

        # define dictionary of mandatory parameters
        self.mandatory_params = {
            'Sequence': None,
            'Samples': None,
            'Bandwidth': None,
            'Bfield': None,
            'Linewidth': 1,
            'TE': None,
            'TE2': None,
            'Add Ref.': 'No',  # default to 'No'
            'Make .raw': 'Yes',  # default to 'Yes' (need for .m script to run properly)
            'Output Path': None,
            'Metabolites': [key for key, value in self.metabs.items() if value],
            'Center Freq': None,   # currently used for plotting on the ppm scale
        }

        # define dictionary of optional parameters
        self.optional_params = {
            'Nucleus': None,
            'TR': None,
        }

    def parse2fidA(self, params):
        # change the parameters to the format used by fidA
        if params['Sequence'] == 'Spin Echo': params['Sequence'] = 'se'
        elif params['Sequence'] == 'PRESS': params['Sequence'] = 'p'
        elif params['Sequence'] == 'STEAM': params['Sequence'] = 'st'
        elif params['Sequence'] == 'LASER': params['Sequence'] = 'l'

        params['Add Ref.'] = params['Add Ref.'][0].lower()
        params['Make .raw'] = params['Make .raw'][0].lower()
        return params

    def run_simulation(self, params, use_multiprocessing=False):
        # create the output directory if it does not exist
        if not os.path.exists(params['Output Path']): os.makedirs(params['Output Path'])

        def sim_lcmrawbasis(n, sw, Bfield, lb, metab, tau1, tau2, addref, makeraw, seq, out_path):
            results = self.octave.feval('sim_lcmrawbasis', n, sw, Bfield, lb, metab,
                                        tau1, tau2, addref, makeraw, seq, out_path + os.sep)
            return metab, results[:, 0] + 1j * results[:, 1]

        # prepare tasks for each metabolite
        params = self.parse2fidA(params)
        tasks = [(params['Samples'], params['Bandwidth'], params['Bfield'],
                  params['Linewidth'], metab, params['TE'], params['TE2'],
                  params['Add Ref.'], params['Make .raw'], params['Sequence'],
                  params['Output Path']) for metab in params['Metabolites']]

        if use_multiprocessing:
            # use multiprocessing to run simulations in parallel
            with Pool() as pool:
                results = pool.starmap(sim_lcmrawbasis_mp, tasks)
        else:
            # run simulations sequentially
            results = [sim_lcmrawbasis(*task) for task in tasks]

        # collect results into a dictionary
        basis_set = {metab: data for metab, data in results}
        return basis_set

    def run_simulation_with_progress(self, params, progress_callback):
        # create the output directory if it does not exist
        if not os.path.exists(params['Output Path']): os.makedirs(params['Output Path'])

        def sim_lcmrawbasis(n, sw, Bfield, lb, metab, tau1, tau2, addref, makeraw, seq, out_path):
            results = self.octave.feval('sim_lcmrawbasis', n, sw, Bfield, lb, metab,
                                        tau1, tau2, addref, makeraw, seq, out_path + os.sep)
            return metab, results[:, 0] + 1j * results[:, 1]

        # prepare tasks for each metabolite
        params = self.parse2fidA(params)
        tasks = [(params['Samples'], params['Bandwidth'], params['Bfield'],
                  params['Linewidth'], metab, params['TE'], params['TE2'],
                  params['Add Ref.'], params['Make .raw'], params['Sequence'],
                  params['Output Path']) for metab in params['Metabolites']]

        # initialize the progress bar
        total_steps = len(tasks)
        progress_step = 100 / total_steps

        # run simulations sequentially
        basis_set = {}
        for i, task in enumerate(tasks):
            metab, data = sim_lcmrawbasis(*task)
            basis_set[metab] = data
            progress_callback(i + 1, total_steps)   # update the progress bar
        return basis_set



if __name__ == "__main__":
    # BasisREMY().run('./example_data/BigGABA_P1P_S01/S01_PRESS_35_act.SPAR', './output/')
    BasisREMY().run_gui()
    
