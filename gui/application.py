####################################################################################################
#                                           application.py                                         #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Created: 08/10/25                                                                                #
#                                                                                                  #
# Purpose: Defines the GUI application for the BasisREMY tool. Each tab is a different             #
#          step in the process, starting with the data selection and REMY extraction,              #
#          continuing with the parameter configuration, and ending with the basis set simulation.  #
#                                                                                                  #
####################################################################################################


#*************#
#   imports   #
#*************#
import matplotlib.pyplot as plt
import numpy as np
import threading
import tkinter as tk

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

from PIL import Image, ImageTk

from tkinter import filedialog, ttk
from tkinterdnd2 import TkinterDnD, DND_FILES
from tkinterweb import HtmlFrame

# own
from core.basisremy import BasisREMY
from gui.help_widget import LabelWithHelp
from gui.export_dialog import ExportDialog


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
        self.windowLength = 800
        self.windowHeight = 700
        self.title("BasisREMY")
        self.geometry("{}x{}".format(self.windowLength, self.windowHeight))
        self.configure(bg=self.bg_1_color)

        # create UI elements
        self.create_widgets()

    def create_widgets(self):
        # create a canvas for the header
        header_canvas = tk.Canvas(self, width=self.windowLength, height=100, highlightthickness=0)
        header_canvas.pack(fill=tk.X)

        # load the header image
        header_image = Image.open("assets/imgs/basisremy_header.png")
        header_image = header_image.resize((1600, 100))
        header_photo = ImageTk.PhotoImage(header_image)

        # place the image on the canvas
        header_canvas.create_image(0, 0, anchor="nw", image=header_photo)
        header_canvas.image = header_photo  # keep a reference to prevent garbage collection

        # add centered text on top of the image
        header_canvas.create_text(
            self.windowLength // 2 + 50, 50,  # x, y position (middle of header + offset for logo)
            text="A tool for generating study-specific basis sets directly from raw MRS data.",
            fill="white",
            font=("Arial", 16, "bold"),
            width=self.windowLength - 300,  # wrap text if too long
            justify="center"
        )

        # add logo
        logo_image = Image.open("assets/imgs/basisremy_logo.png")
        logo_image = logo_image.resize((100, 100))
        logo_photo = ImageTk.PhotoImage(logo_image)
        logo_label = tk.Label(self, image=logo_photo, bg="#f0f0f0")
        logo_label.image = logo_photo  # keep a reference to prevent garbage collection
        logo_label.place(x=0, y=0)

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

        # Track whether the current basis set is still valid for the params
        # on screen.  Set True after a successful simulation, False whenever
        # the user navigates back to tab2 (edit params) from tab3.
        self._basis_set_valid = False

        # Watch tab switches so we can invalidate when going backwards.
        self.notebook.bind("<<NotebookTabChanged>>", self._on_tab_changed)

        self.tab1_widgets()
        self.tab2_widgets()
        self.tab3_widgets()

        # apply custom styles to tabs
        style = ttk.Style()
        style.configure("TNotebook.Tab",
                        font=("Arial", 12))  # set the font size and style for the tabs

        self.update()
        # self.geometry(f"{self.winfo_reqwidth()}x{self.winfo_reqheight()}")

    # ------------------------------------------------------------------ tab navigation
    def _on_tab_changed(self, event=None):
        """Called whenever the notebook switches tabs.

        Rule: if the user navigates *back* to tab2 (Parameter Configuration)
        from tab3 (Basis Simulation), the previously computed basis set is no
        longer guaranteed to match the current parameters — mark it stale and
        lock tab3 so they must re-simulate before going forward again.
        """
        try:
            current = self.notebook.index(self.notebook.select())
        except Exception:
            return

        if current == 1:
            # Arrived at Parameter Configuration — invalidate any existing results.
            if self._basis_set_valid:
                self._basis_set_valid = False
                # Disable tab3 so the user cannot jump straight back to stale results.
                self.notebook.tab(2, state="disabled")
                # Re-enable the Simulate button (it may have been hidden if we got here
                # programmatically; validate_inputs will set the correct state).
                try:
                    self.validate_inputs()
                except Exception:
                    pass

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
        # Reset per-field state that lives on `self`. Without this, switching
        # backend leaves stale StringVars from the previous backend's pulse
        # pickers in `self.file_vars`, which can be picked up by lambdas
        # whose widgets have already been destroyed.
        self.file_vars = {}

        # backend toggle — two-level cascade: Category → Backend
        backend_frame = tk.Frame(self.tab2)
        backend_frame.pack(fill=tk.X, padx=10, pady=(10, 0))

        current_category = self.BasisREMY.get_current_category()
        # Only show categories that have at least one backend registered.
        category_options = [c for c in self.BasisREMY.CATEGORY_ORDER
                            if self.BasisREMY.categories.get(c)]
        # Append any extra categories that may have been registered out of band
        for c in self.BasisREMY.categories:
            if c not in category_options and self.BasisREMY.categories[c]:
                category_options.append(c)

        tk.Label(backend_frame, text="Category:",
                 font=("Arial", 12, "bold")).pack(side=tk.LEFT, padx=(0, 5))
        self.category_var = tk.StringVar(value=current_category)
        category_combo = ttk.Combobox(backend_frame, textvariable=self.category_var,
                                      values=category_options, font=("Arial", 12),
                                      state="readonly", width=12)
        category_combo.pack(side=tk.LEFT, padx=(0, 15))

        # Backend label + combo — only shown when the chosen category has
        # more than one backend option. Single-backend categories (MRSCloud,
        # FSL-MRS, Custom) don't need a sub-selector.
        backend_label_widget = tk.Label(backend_frame, text="Backend:",
                                        font=("Arial", 12, "bold"))
        backend_label_widget.pack(side=tk.LEFT, padx=(0, 5))

        # Display labels in the Backend combo use `display_name`, while we
        # keep a mapping back to the canonical `name` for set_backend().
        def _backends_for(cat):
            names = self.BasisREMY.categories.get(cat, [])
            label_to_name = {}
            labels = []
            for n in names:
                b = self.BasisREMY.backends[n]
                label = getattr(b, 'display_name', None) or b.name
                label_to_name[label] = n
                labels.append(label)
            return labels, label_to_name

        def _set_backend_combo_visibility(labels):
            """Show or hide the Backend widgets depending on choice count."""
            if len(labels) > 1:
                backend_label_widget.pack(side=tk.LEFT, padx=(0, 5))
                backend_combo.pack(side=tk.LEFT)
            else:
                backend_label_widget.pack_forget()
                backend_combo.pack_forget()

        backend_labels, self._backend_label_map = _backends_for(current_category)
        current_label = next(
            (lbl for lbl, nm in self._backend_label_map.items()
             if nm == self.BasisREMY.backend.name),
            backend_labels[0] if backend_labels else '',
        )
        self.backend_var = tk.StringVar(value=current_label)
        backend_combo = ttk.Combobox(backend_frame, textvariable=self.backend_var,
                                     values=backend_labels, font=("Arial", 12),
                                     state="readonly", width=32)
        backend_combo.pack(side=tk.LEFT)
        # Apply initial visibility
        _set_backend_combo_visibility(backend_labels)

        def _do_switch(target_name):
            """Switch backend, checking Octave availability if needed."""
            if target_name == self.BasisREMY.backend.name:
                return True
            new_backend = self.BasisREMY.backends[target_name]
            if new_backend.requires_octave and new_backend.octave is None:
                if not self.check_octave_availability():
                    return False
            self.BasisREMY.set_backend(target_name)
            return True

        def on_category_change(event=None):
            cat = self.category_var.get()
            labels, label_map = _backends_for(cat)
            if not labels:
                return
            backend_combo['values'] = labels
            self._backend_label_map = label_map
            self.backend_var.set(labels[0])
            _set_backend_combo_visibility(labels)
            target_name = label_map[labels[0]]
            if _do_switch(target_name):
                self.tab2_widgets()
            else:
                self.category_var.set(self.BasisREMY.get_current_category())

        def on_backend_change(event=None):
            label = self.backend_var.get()
            target_name = self._backend_label_map.get(label)
            if target_name is None:
                return
            if _do_switch(target_name):
                self.tab2_widgets()
            else:
                cur = next((lbl for lbl, nm in self._backend_label_map.items()
                            if nm == self.BasisREMY.backend.name), label)
                self.backend_var.set(cur)

        category_combo.bind("<<ComboboxSelected>>", on_category_change)
        backend_combo.bind("<<ComboboxSelected>>", on_backend_change)

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
            val = var.get()
            # Store in mandatory or optional params depending on where the key lives
            if key in self.BasisREMY.backend.mandatory_params:
                self.BasisREMY.backend.mandatory_params[key] = val
            elif key in self.BasisREMY.backend.optional_params:
                self.BasisREMY.backend.optional_params[key] = val
            self.validate_inputs()

            # Some keys change which OTHER fields are visible (e.g. picking
            # MEGA reveals editing fields, picking Philips reveals the
            # vendor pulse picker). Backends declare these in
            # `schema_affecting_keys`; rebuild the panel so the UI matches.
            schema_keys = getattr(self.BasisREMY.backend,
                                  'schema_affecting_keys', set())
            if key in schema_keys:
                self.tab2_widgets()

        # Mode selector - only shown if backend has more than one mode
        backend = self.BasisREMY.backend
        has_modes = len(backend.modes) > 1

        if has_modes:
            mode_frame = tk.Frame(params_frame)
            mode_frame.grid(row=0, column=0, columnspan=3, sticky="ew", pady=(0, 10))

            tk.Label(mode_frame, text="Mode:", font=("Arial", 12, "bold")).pack(side=tk.LEFT, padx=(0, 10))

            self.mode_var = tk.StringVar(value=backend.current_mode)
            mode_combo = ttk.Combobox(mode_frame, textvariable=self.mode_var,
                                      values=backend.modes, font=("Arial", 12), state="readonly")
            mode_combo.pack(side=tk.LEFT)

            def on_mode_change(event=None):
                new_mode = self.mode_var.get()
                backend.set_mode(new_mode)
                self.tab2_widgets()  # Rebuild with new mode params

            mode_combo.bind("<<ComboboxSelected>>", on_mode_change)
            row = 1
        else:
            row = 0

        # Get parameters to display based on current mode
        params_to_show = backend.get_params_for_mode()

        # populate the parameters frame
        for key, value in params_to_show.items():
            if key in self.BasisREMY.backend.file_selection:
                # label + help icon for the path (file)
                label = LabelWithHelp(params_frame, key, bg=self.bg_1_color)
                label.grid(row=row, column=0, padx=0, pady=5, sticky="e")

                # stringVar for the selected path. We keep a per-key dict so
                # backends with multiple pulse-file pickers (e.g. MEGA-PRESS:
                # refoc + edit) don't share a single StringVar. `self.file_var`
                # remains for backwards compat; it always points at the most
                # recently created field.
                if not hasattr(self, 'file_vars'):
                    self.file_vars = {}
                self.file_var = tk.StringVar(value="missing input")
                self.file_vars[key] = self.file_var

                entry = tk.Entry(params_frame, textvariable=self.file_var, font=("Arial", 12))
                entry.grid(row=row, column=1, padx=(10, 0), pady=5, sticky="ew")

                # IMPORTANT: bind `key` and `var` as default args so the
                # closure captures THIS iteration's values. Without this
                # the loop's `key` reference is shared across all picker
                # buttons and Browse always writes to the loop's last
                # parameter (the infamous "Tau 2 := pulse path" bug).
                def select_path(key=key, var=self.file_var):
                    file_path = filedialog.askopenfilename()
                    if file_path:
                        var.set(file_path)
                        self.BasisREMY.backend.mandatory_params[key] = file_path
                    else:
                        var.set("missing input")
                        self.BasisREMY.backend.mandatory_params[key] = None
                    self.validate_inputs()

                button = tk.Button(params_frame, text="Browse", command=select_path)
                button.grid(row=row, column=2, padx=0, pady=5)

                # automatically update backend dict on any entry edit
                self.file_var.trace_add('write', lambda *args, key=key,
                                                        var=self.file_var: update_param(key, var))

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
                # label + help icon for dropdown parameters
                label = LabelWithHelp(params_frame, key, bg=self.bg_1_color)
                label.grid(row=row, column=0, padx=0, pady=5, sticky="e")

                # StringVar for the Combobox. Use "Select option" as the
                # placeholder when nothing is set so the user knows they
                # have to pick one (instead of the cryptic "missing input").
                allowed = list(self.BasisREMY.backend.dropdown[key])
                if value is not None and str(value) in allowed:
                    initial = str(value)
                else:
                    initial = "Select option"
                var = tk.StringVar(value=initial)

                # Combobox for dropdown parameters. Readonly so users can
                # only choose values the backend understands (free-typing
                # caused MRSCloud crashes when 'PRESS' was typed into the
                # Sequence box, which expects UnEdited/MEGA/HERMES/HERCULES).
                combobox = ttk.Combobox(params_frame, textvariable=var,
                                        values=allowed,
                                        font=("Arial", 12), state="readonly")
                combobox.grid(row=row, column=1, padx=(10, 0), pady=5, sticky="ew")

                # trace changes to the StringVar
                var.trace_add('write', lambda *args, key=key, var=var: update_param(key, var))

            else:
                # label + help icon for other parameters
                label = LabelWithHelp(params_frame, key, bg=self.bg_1_color)
                label.grid(row=row, column=0, padx=0, pady=5, sticky="e")

                # StringVar for the Entry
                var = tk.StringVar(value=str(value) if value is not None else "missing input")

                # entry for other parameters
                entry = tk.Entry(params_frame, textvariable=var, font=("Arial", 12))
                entry.grid(row=row, column=1, padx=(10, 0), pady=5, sticky="ew")

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

        # All widgets built and bound — refresh the simulate-button state so
        # it lights up immediately when defaults / REMY-parsed values already
        # satisfy validation (otherwise the user would have to touch a field
        # for nothing).
        try:
            self.validate_inputs()
        except Exception:
            pass

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
            filetypes=[("All Files", "*")]
        #     filetypes=[
        #         ("MRS Data Files", (
        #             "*.7", "*.ima", "*.rda", "*.dat",
        #             "*.spar", "*method", "*2dseq", "*.nii", "*.nii.gz"
        #         ))
        #     ]
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
            params, opt = self.BasisREMY.backend.parseREMY(MRSinMRS)

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

    def check_octave_availability(self):
        """
        Check if Octave is available (Docker or local).
        Show a dialog with instructions if not available.

        Returns:
            bool: True if Octave is available, False otherwise
        """
        from core.octave_manager import OctaveManager

        manager = OctaveManager()
        docker_available = manager.check_docker_availability()
        local_available = manager.check_local_octave_availability()

        if docker_available or local_available:
            return True

        # Neither is available - show helpful dialog
        instructions = manager._get_installation_instructions()

        # Create a custom dialog with the instructions
        dialog = tk.Toplevel(self)
        dialog.title("Octave Runtime Required")
        dialog.geometry("700x500")
        dialog.configure(bg=self.bg_1_color)

        # Make it modal
        dialog.transient(self)
        dialog.grab_set()

        # Add text widget with scrollbar for instructions
        frame = tk.Frame(dialog, bg=self.bg_1_color)
        frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        scrollbar = tk.Scrollbar(frame)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        text = tk.Text(frame, wrap=tk.WORD, yscrollcommand=scrollbar.set,
                      font=("Courier", 10), bg="white", fg="black")
        text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.config(command=text.yview)

        text.insert("1.0", instructions)
        text.config(state=tk.DISABLED)

        # Add close button
        close_btn = tk.Button(dialog, text="OK", command=dialog.destroy,
                            bg=self.main_color, fg=self.text_color,
                            font=("Arial", 12, "bold"))
        close_btn.pack(pady=10)

        return False

    def validate_inputs(self):
        # check if all mandatory parameters are filled
        _UNSET = (None, "", "missing input", "Select option")
        all_params_filled = all(
            self.BasisREMY.backend.mandatory_params[key] not in _UNSET
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
        # Check if Octave is available before proceeding
        if self.BasisREMY.backend.requires_octave and self.BasisREMY.backend.octave is None:
            if not self.check_octave_availability():
                return  # Don't proceed if Octave is not available

        # Always reset tab3 to the clean progress-bar state before starting a
        # new run so old results (plot + checkboxes) are never visible while a
        # new simulation is in flight.
        self._basis_set_valid = False
        self.tab3_widgets()

        # move to the next tab
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
        basis = self.BasisREMY.backend.run_simulation(self.BasisREMY.backend.mandatory_params, progress_callback)
        self.after(0, callback, basis)

    def on_simulation_complete(self, basis):
        print("Simulation complete.")
        self.basis_set = basis
        # Mark the basis set as valid for the current parameters.
        # _on_tab_changed will clear this flag if the user goes back to edit params.
        self._basis_set_valid = True

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

        # Export… button (unified exporter, see core/exporters.py)
        self.export_button = tk.Button(
            self.tab3,
            text="Export Basis…",
            command=self.open_export_dialog,
            bg=self.main_color,
            fg=self.text_color,
            font=("Arial", 12, "bold"),
        )
        self.export_button.pack(pady=(8, 0))

        # render the initial plot.
        self.update_plot()

    def open_export_dialog(self):
        """Open the unified export dialog (LCModel / jMRUI / FSL-MRS / Osprey)."""
        if not getattr(self, "basis_set", None):
            return
        ExportDialog(self, self.basis_set, self.BasisREMY.backend.mandatory_params)

    def update_plot(self):
        # clear the axis and reapply basic settings.
        self.ax.clear()
        # self.ax.set_title("Basis Simulation Results")
        self.ax.set_facecolor(self.bg_1_color)
        self.ax.set_xlabel("Chemical Shift [ppm]")

        # no y-axis ticks or labels.
        self.ax.set_yticks([])
        self.ax.set_yticklabels([])

        # compute ppm axis using cf (centre / Larmor frequency in Hz).
        # Some backends (e.g. MRSCloud) intentionally don't expose Center Freq
        # because it's derived from Field Strength internally — fall back to
        # γ × B0 when the explicit value is missing.
        mp = self.BasisREMY.backend.mandatory_params
        cf_raw = mp.get('Center Freq')
        if cf_raw in (None, '', 'missing input'):
            field_str = str(mp.get('Field Strength') or '3T').replace('T', '').strip()
            try:
                b0 = float(field_str)
            except ValueError:
                b0 = 3.0
            cf = 42.577e6 * b0  # Hz
        else:
            cf = float(cf_raw) * (1e6 if float(cf_raw) < 1000 else 1.0)
        bw = float(mp['Bandwidth'])
        # We use the actual FID length (data.size) per metabolite when
        # building the ppm axis below — backends like MRSCloud may return a
        # different length than `Samples` if their internal grid disagrees.

        # plot each metabolite's data if its checkbox is selected.
        for metab, var in self.checkbox_vars.items():
            if var.get() and metab in self.basis_set:
                data = self.basis_set[metab]

                # Ensure data is a proper numpy array
                if not isinstance(data, np.ndarray):
                    try:
                        data = np.array(data, dtype=complex)
                    except:
                        print(f"Warning: Could not convert {metab} data to numpy array")
                        continue

                # Flatten if needed
                if data.ndim > 1:
                    data = data.flatten()

                # Check if data is empty
                if data.size == 0:
                    print(f"Warning: {metab} data is empty")
                    continue

                try:
                    # MRS convention: FFT the FID, then fftshift so frequency
                    # 0 is centred. The ppm axis goes low→high; the x-axis is
                    # inverted below (xlim 10→0) so high ppm appears on the
                    # left, matching radiological MRS display convention.
                    # NOTE: do NOT flip ppm_axis — flipping creates a spectral
                    # mirror by misassigning FFT bins to the wrong ppm values.
                    ydata = np.real(np.fft.fftshift(np.fft.fft(data)))
                    npts = data.size
                    # ppm axis: offset from carrier → absolute chemical shift.
                    # The offset (+4.65) must match the centreFreq used during
                    # simulation (MRSCloud adapter sets centreFreq = 4.65, the
                    # standard water reference). A mismatch causes a rigid shift.
                    ppm_axis = np.linspace(-bw / 2, bw / 2, npts) / cf * 1e6 + 4.65
                    self.ax.plot(ppm_axis, ydata, color=self.metab_colors[metab])
                except Exception as e:
                    print(f"Warning: Could not plot {metab}: {e}")
                    continue

        # MRS convention: high ppm on the left, low ppm on the right.
        # set_xlim(10, 0) with xmin > xmax is all matplotlib needs — do NOT
        # also call invert_xaxis() because that would invert a second time
        # and put the axis back to 0→10.
        self.ax.set_xlim(10, 0)

        self.canvas.draw()