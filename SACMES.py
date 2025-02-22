""" SACMES_SWV.py
Monitor and analyze data from a potentiostat.
"""
#style conventions in code cleanup
# (beyond standard Python conventions, https://peps.python.org/pep-0008/):
# double-quoted strings
# prefix global variables, of which there are many, with global_
# pylint: disable=C0103
# (so long as we have global variables, we ask pylint not to worry about their names)
from typing import NoReturn, Tuple, Union, List, Dict, Optional, Sequence, Set, Any
                                    ########################
                                    ### Import Libraries ###
                                    ########################                                   
import sys
import os
match os.name:
    case 'posix': # Linux/MacOS
        import fcntl
    case 'nt': # Windows
        import msvcrt   
import time
import datetime
import tkinter as tk
from tkinter import filedialog
from tkinter import ttk
import csv
from math import sqrt
from threading import Thread
from queue import Queue, Empty
from enum import Enum
import warnings
import psutil
import numpy as np
from scipy.signal import savgol_filter #type: ignore
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.figure
import matplotlib.lines
import matplotlib.artist
import matplotlib.axes
from matplotlib.patches import Polygon
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import casadi as ca

matplotlib.use("TkAgg")
plt.style.use("ggplot")
#---Clear mac terminal memory---# #TODO: what problem does this solve? is it portable?
#os.system("clear && printf '\e[3J'")
#---Filter out error warnings---#
warnings.simplefilter("ignore", np.RankWarning)         #numpy polyfit_deg warning
warnings.filterwarnings(action="ignore", module="scipy", message="^internal gelsd") #RuntimeWarning
                                    ########################
                                    ### Global Variables with Initializers (including constants)
                                    ### (the remaining global variables are declared after all
                                    # class declarations and before the main program body)
                                    ########################
global_handle_variable: str = "" #file name prefix for input data files
class ElectrodesMode(Enum):
    SINGLE = 0
    MULTIPLE = 1
global_electrodes_mode: ElectrodesMode = ElectrodesMode.SINGLE
class PlotSummaryMode(Enum):
    PHE = 0
    AUC = 1
    def __str__(self) -> str:
        match self:
            case PlotSummaryMode.PHE:
                return "Peak Height Extraction"
            case PlotSummaryMode.AUC:
                return "Area Under the Curve"
global_plot_summary_mode: PlotSummaryMode = PlotSummaryMode.PHE
class PlotTimeReportingMode(Enum):
    EXPERIMENT_TIME = 0
    FILE_NUMBER = 1
    def __str__(self) -> str:
        match self:
            case PlotTimeReportingMode.EXPERIMENT_TIME:
                return "Experiment Time"
            case PlotTimeReportingMode.FILE_NUMBER:
                return "File Number"
STARTING_FILE_NUMBER: int = 1
DEFAULT_NORMALIZATION_FILE_NUMBER: int = 1
# frequencies initially displayed in Frequency Listbox
global_input_frequencies: List[int] = [5, 10, 15, 30, 40, 60, 100, 150, 200]
ELECTRODES_NUMBERS: List[int] = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
########################################
### Polynomial Regression Parameters ###
########################################
### Savitzky-Golay window (in mV range):
global_savitzky_golay_window: float = 5.0
### Savitzky-Golay polynomial degree:
global_savitzky_golay_degree: int = 3
POLYFIT_DEGREE: int = 15        ### degree of polynomial fit
CUTOFF_FREQUENCY: int = 50      ### frequency that separates "low" and "high"
                                ### frequencies for regression analysis and
                                ### smoothing manipulation
###############################
### Gauss Method Parameters ###
###############################
global_gauss_peak: List[List[float]]
global_gauss_baseline: List[List[float]]
global_gauss_maxheight: List[List[float]]
global_low_gauss_peak: float = -0.3
global_high_gauss_peak: float = -0.3
global_low_gauss_baseline: float = 0.1
global_high_gauss_baseline: float = 0.3
global_low_gauss_maxheight: float = 5
global_high_gauss_maxheight: float = 5

global_gauss_solver : List[None]
global_options_Gauss: dict = {
    'ipopt.print_level': 0,
    'print_time': 0,
    'ipopt.tol': 1e-8,
    'ipopt.max_iter': 5000
}
#############################
### Checkpoint Parameters ###
#############################
global_key: int = 0                 ### SkeletonKey
global_poison_pill: bool = False      ### Stop Animation variable
global_found_file_path: bool = False   ### If the user-inputted file is found
global_frequency_warning_label_exists: bool = False
global_analysis_already_initiated: bool = False
global_high_already_reset: bool = False    ### If data for high frequencies has been reset
global_low_already_reset: bool = False      ### If data for low frequencies has been reset
global_already_reset: bool = False
global_analysis_complete: bool = False    ### If analysis has completed, begin PostAnalysis
##################################
### Data Extraction Parameters ###
##################################
global_file_name_pattern: str = "<H><E>_<F>Hz__<N>.txt"
global_current_column: int = 2
global_base_column_index_for_currents: int = global_current_column - 1
global_voltage_column: int = 1
global_voltage_column_index: int = global_voltage_column - 1
global_columns_per_electrode: int = 3
#-- set the initial limit in bytes to filter out preinitialized files < 3000b
#global_byte_limit: int = 2000
#- set the initial byte index to match the checkbutton
#- index in the toolbar menu MainWindow.byte_menu
global_byte_index: int = 1
######################################################
### Low frequency baseline manipulation Parameters ###
######################################################
global_low_frequency_offset: float = 0         ### Vertical offset of normalized data for
                               ### user specified "Low Frequency"
global_low_frequency_slope: float = 0          ### Slope manipulation of norm data for user
                               ### specified "Low Frequency"
###############
### Styling ###
###############
HUGE_FONT = ("Verdana", 18)
LARGE_FONT = ("Verdana", 12)
MEDIUM_FONT = ("Verdana", 10)
SMALL_FONT = ("Verdana", 8)
global_file_path: str = "" #path to the directory containing the input data files
class XBound(Enum):
    START_PLUS = 1
    NOW_MINUS = 2
global_x_left_bound_radiobutton: XBound = XBound.START_PLUS
global_x_right_bound_radiobutton: XBound = XBound.NOW_MINUS
global_x_left_bound_offset: float = 0
global_x_right_bound_offset: float = 0
class YBound(Enum):
    AUTOMATIC = 1
    MANUAL = 2
global_y_norm_radiobutton: YBound = YBound.AUTOMATIC
global_y_kdm_radiobutton: YBound = YBound.AUTOMATIC
FRAME_POST_ANALYSIS: str = "PostAnalysis"
FRAME_LOW_PARAMETER: str = "LowParameter"
FRAME_HIGH_PARAMETER: str = "HighParameter"
FRAME_INPUT: str = "InputFrame"
EPSILON: float = 0.0000001
class Delimiter(Enum):
    SPACE = 1
    TAB = 2
    COMMA = 3
    def __str__(self) -> str:
        match self:
            case Delimiter.SPACE:
                return " "
            case Delimiter.TAB:
                return "\t"
            case Delimiter.COMMA:
                return ","
global_delimiter: Delimiter = Delimiter.SPACE
class FileEncoding(Enum):
    UTF_8 = 1
    UTF_16 = 2
    def __str__(self) -> str:
        match self:
            case FileEncoding.UTF_8:
                return "UTF-8"
            case FileEncoding.UTF_16:
                return "UTF-16"
global_file_encoding: FileEncoding = FileEncoding.UTF_8
                        ########################
                        ### Global Functions ###
                        ########################
def internal_error(internal_error_message: str) -> NoReturn:
    """To be used for situations that should not arise ever.
    """
    print("internal_error:", internal_error_message)
    sys.exit(1)

def file_is_complete(filename: str) -> bool:
    """Heuristic test if the file is complete in the sense
    that the data acquisition software is no longer writing
    the file. In the original SACMES code, this heuristic was
    based on file length, which is not reliable; indeed it is
    wrong in the normal mode of operation.
    In the new implementation, one we know the file exists
    we check if it is open by any running process. If not,
    we consider that it is complete."""
    #return os.path.exists(filename) and os.path.getsize(filename) > global_byte_limit
    #return True
    #Trying a new method since psutil.process_iter() takes a long time.
    if os.path.exists(filename):
        match os.name:
            case 'posix': # Linux/MacOS
                try:
                    with open(filename, 'r+') as f:
                        fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)  # Try locking the file
                        fcntl.flock(f, fcntl.LOCK_UN)  # Unlock immediately
                    return True  # Locking succeeded, file is not in use
                except IOError:
                    return False # File is locked, wait and retry\
            case 'nt': # Windows
                try:
                    with open(filename, 'r+') as f:
                        msvcrt.locking(f.fileno(), msvcrt.LK_NBLCK, 1)  # Try locking the file
                        msvcrt.locking(f.fileno(), msvcrt.LK_UNLCK, 1)  # Unlock immediately
                    return True  # Locking succeeded, file is not in use
                except OSError:
                    return False # File is locked, wait and retry\
        # for proc in psutil.process_iter():
        #     print("STILL SEARCHING PROCESSES")
        #     try:
        #         for file in proc.open_files():
        #             if file.path == filename:
        #                 #print(f"{filename} open by {proc.pid}")
        #                 return False
        #     except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess):
        #         pass
        # #print(f"{filename} assumed not open; length={os.path.getsize(filename)}")
        # return True
    else:
        #print(f"{filename} does not exist")
        return False

def get_time(file_index: int) -> float:
    if file_index == 0:
        file_index = 1
    file_now = global_file_path + make_file_name(file_index,1,global_low_frequency)
    file_1 = global_file_path + make_file_name(1,1,global_low_frequency)
    exp_time = (os.path.getmtime(file_now) - os.path.getmtime(file_1))/3600 # in hours
    return exp_time

def make_file_name(file_index: int, electrode: int, frequency: int) -> str:
    """Instantiate the file name pattern with the given values."""
    name: str = global_file_name_pattern
    match global_electrodes_mode:
        case ElectrodesMode.SINGLE:
            name = name.replace("<H>", global_handle_variable)
            name = name.replace("<E>", "1")
            name = name.replace("<F>", str(frequency))
            name = name.replace("<N>", str(file_index))
        case ElectrodesMode.MULTIPLE:
            name = name.replace("<H>", global_handle_variable)
            name = name.replace("<E>", str(electrode))
            name = name.replace("<F>", str(frequency))
            name = name.replace("<N>", str(file_index))
    return name

def read_data(input_file_name: str, electrode: int) ->\
    Tuple[List[float], List[float], Dict[float, float]]:
    """Extract numerical data for potentials and currents
    from an instrument-produced data file,
    and return them as a tuple
    (potentials, currents, potential_to_current_map).
    """
    potentials: List[float]
    currents: List[float]
    potential_to_current_map: Dict[float, float]
    try:
        with open(input_file_name, "r", encoding=str(global_file_encoding)) as mydata:
            potentials = []
            currents = []
            potential_to_current_map = {}
            for line in mydata:
                check_split_list = line.split(str(global_delimiter))
                while check_split_list[0] == " " or check_split_list[0] == "\t":
                    del check_split_list[0]
                check_split_first_item: str = check_split_list[0].replace(",", "")
                first_item_is_float: bool
                try:
                    float(check_split_first_item)
                    first_item_is_float = True
                except ValueError:
                    first_item_is_float = False
                if first_item_is_float:
                    current_value: float =\
                        1000000 * float(check_split_list[column_index_for_current(electrode)].\
                                        replace(",", ""))
                    currents.append(current_value)
                    potential_value: float =\
                        float(line.split(str(global_delimiter))[global_voltage_column_index].\
                              strip(","))
                    potentials.append(potential_value)
                    potential_to_current_map[potential_value] = current_value
        return potentials, currents, potential_to_current_map
    except FileNotFoundError as exception:
        internal_error("read_data: file not found " + str(exception))

#######################################
### Retrieve the column index value ###
#######################################
def column_index_for_current(electrode: int) -> int:
    """Depending on the type of instrument output file
    (single file for all electrodes, or multiple files),
    guess the column with the data for the given electrode.
    """
    match global_electrodes_mode:
        case ElectrodesMode.SINGLE:
            return global_base_column_index_for_currents +\
                  (electrode-1)*global_columns_per_electrode
        case ElectrodesMode.MULTIPLE:
            return global_base_column_index_for_currents

def _update_global_lists(file: int):
    """Record the file number and sample rate in global lists."""
    if file not in global_file_list:
        global_file_list.append(file)
        sample: float = round(get_time(len(global_file_list)), 3)
        global_sample_list.append(sample)
        global_real_time_sample_label.config(text=str(sample))
        if file != global_number_of_files_to_process:
            global_file_label.config(text=str(file + 1))

def compute_x_bounds(file_index: int) -> Tuple[float, float]:
    """Return the start and end of the time period to be displayed in ratiometric plots.
    """
    left: float #hours OR file number
    right: float #hours OR file number
    match global_x_axis_mode:
        case PlotTimeReportingMode.EXPERIMENT_TIME:
            now = get_time(file_index)
        case PlotTimeReportingMode.FILE_NUMBER:
            now = file_index
    match global_x_left_bound_radiobutton:
        case XBound.START_PLUS:
            left = global_x_left_bound_offset
        case XBound.NOW_MINUS:
            left = now - global_x_left_bound_offset
    match global_x_right_bound_radiobutton:
        case XBound.START_PLUS:
            right = global_x_right_bound_offset
        case XBound.NOW_MINUS:
            right = now - global_x_right_bound_offset
    if left < 0:
        left = 0
    if right < left:
        right = left + 1
    return left, right

def MultiGausFitNoise_LSE(length, options):
    # Define symbolic variables
    c = ca.SX.sym('c', 3)
    mu = ca.SX.sym('mu', 3)
    lambda_ = ca.SX.sym('var', 3)
    L0 = ca.SX.sym('L0')
    noise = ca.SX.sym('noise', length)
    anoise = ca.SX.sym('anoise', length)  
    yB = ca.SX.sym('yB', length)
    v = ca.SX.sym('v', length)
    theta = ca.vertcat(c, mu, lambda_, L0, noise, anoise)
    xB = L0 + c[0] * ca.exp((-(v - mu[0]) ** 2) * lambda_[0]) + c[1] * ca.exp((-(v - mu[1]) ** 2) * lambda_[1]) + c[2] * ca.exp((-(v - mu[2]) ** 2) * lambda_[2]) + noise
    eB = yB - xB
    # Define the objective function
    J = ca.mtimes(eB.T, eB)
    param = ca.vertcat(yB, v)
    constraints = ca.vertcat(
        c[2] * ca.exp((-(mu[0] - mu[2]) ** 2) * lambda_[2]) - 0.05 * c[0],
        c[1] * ca.exp((-(mu[0] - mu[1]) ** 2) * lambda_[1]) - 0.05 * c[0],
        noise - anoise,
        noise + anoise,
        ca.sum1(anoise))
    nlp = {'x': theta, 'f': J, 'g': constraints, 'p': param}
    # Create the solver  
    solver = ca.nlpsol('NLP_canon', 'ipopt', nlp, options)
    return solver
                        ###############################################
                        ########   Graphical User Interface     #######
                        ###############################################
class MainWindow(tk.Tk):
    """A class that contains all of the frames for the GUI.
    """
    #--- Initialize the GUI ---#
    def __init__(self, master: tk.Tk, *args, **kwargs):
        global global_container,\
            global_show_frames,\
            global_high_low_dictionary
        tk.Tk.__init__(self, *args, **kwargs)
        self.withdraw()
        self.master: tk.Tk = master
        self.master.wm_title("SACMES")
        self.master.rowconfigure(0, weight=1)
        self.master.columnconfigure(0, weight=1)
        #--- Create a frame for the UI ---#
        global_container = ttk.Frame(self.master, relief="flat")
        ## container object has UI frame in column 0:
        global_container.grid(row=0, rowspan=11, padx=10, sticky="nsew")
        ## and PlotContainer (visualization) in column 1:
        global_container.rowconfigure(0, weight=1)
        global_container.columnconfigure(0, weight=1)
        #--- Raise the frame for initial UI ---#
        # Key: frame handle / Value: tk.Frame object:
        global_show_frames = {}
        frame = InputFrame(global_container, self.master)
        global_show_frames[FRAME_INPUT] = frame
        frame.grid(row=0, column=0, sticky="nsew")
        self.show_frame(FRAME_INPUT)
        self.create_toolbar()
        #--- High and Low Frequency Dictionary ---#
        global_high_low_dictionary = {}
        self.list_val_entry: ttk.Entry
        self.voltage_column: ttk.Entry
        self.spacing_val_entry: ttk.Entry
        self.space_delimiter: ttk.Radiobutton
        self.tab_delimiter: ttk.Radiobutton
        self.txt_value: ttk.Radiobutton
        self.csv_value: ttk.Radiobutton
        self.dta_value: ttk.Radiobutton
        self.file_name_pattern_entry: ttk.Entry

    def create_toolbar(self) -> None:
        """Create the toolbar for the GUI."""
        menubar = tk.Menu(self.master)
        self.master.config(menu=menubar)
        editmenu = tk.Menu(menubar, tearoff=0)
        editmenu.add_separator()
        editmenu.add_command(label="Customize File Format",\
                             command=self.extraction_adjustment_frame)
        self.delimiter_value = tk.IntVar()
        self.delimiter_value.set(1)
        self.extension_value = tk.IntVar()
        self.extension_value.set(1)
        self.file_encoding_value = tk.IntVar()
        self.file_encoding_value.set(1)
        # self.byte_menu = tk.Menu(menubar)
        # self.byte_menu.add_command(label="   1000", command=lambda: self.set_bytes(1000, 0))
        # self.byte_menu.add_command(label="✓  2000", command=lambda: self.set_bytes(2000, 1))
        # self.byte_menu.add_command(label="   3000", command=lambda: self.set_bytes(3000, 2))
        # self.byte_menu.add_command(label="   4000", command=lambda: self.set_bytes(4000, 3))
        # self.byte_menu.add_command(label="   5000", command=lambda: self.set_bytes(5000, 4))
        # self.byte_menu.add_command(label="   7500", command=lambda: self.set_bytes(7500, 5))
        # editmenu.add_cascade(label="Byte Limit", menu=self.byte_menu)
        menubar.add_cascade(label="Settings", menu=editmenu)

    def extraction_adjustment_frame(self) -> None:
        """Create a frame for adjusting the file extraction parameters.
        """
        win = tk.Toplevel()
        win.wm_title("Customize File Format")
        #-- new frame --#
        row_value = 0
        container = ttk.Frame(win, relief="groove", padding=5)
        container.grid(row=row_value, column=0, columnspan=2, padx=5, pady=5, ipadx=3)
        container_value = 0
        ttk.Label(container, text="Current is in column #:").grid(row=container_value, column=0)
        container_value += 1
        self.list_val_entry = ttk.Entry(container, width=5)
        self.list_val_entry.insert(tk.END, str(global_current_column))
        self.list_val_entry.grid(row=container_value, column=0, pady=5)
        container_value = 0
        ttk.Label(container, text="Voltage is in column #:").grid(row=container_value, column=1)
        container_value += 1
        self.voltage_column = ttk.Entry(container, width=5)
        self.voltage_column.insert(tk.END, str(global_voltage_column))
        self.voltage_column.grid(row=container_value, column=1, pady=5)
        container_value += 1
        ttk.Label(container, text="Multipotentiostat\ncolumns per electrode:").\
            grid(row=container_value, column=0, columnspan=2)
        container_value += 1
        inner_frame = ttk.Frame(container, padding=5)
        inner_frame.grid(row=container_value, column=0, columnspan=2)
        ttk.Label(inner_frame, text="\t         Columns").grid(row=0, column=0)
        self.spacing_val_entry = ttk.Entry(inner_frame, width=4)
        self.spacing_val_entry.insert(tk.END, str(global_columns_per_electrode))
        self.spacing_val_entry.grid(row=0, column=0, pady=1)
        #-- new frame --#
        row_value += 1
        box = ttk.Frame(win, relief="groove", padding=5)
        box.grid(row=row_value, column=0, columnspan=2, pady=7)
        box_value = 0
        ttk.Label(box, text="Delimiter between\ndata columns:").grid(row=box_value, column=0)
        box_value += 1
        ttk.Radiobutton(box, text="Space", variable=self.delimiter_value, value=1).\
            grid(row=box_value, column=0, pady=5)
        box_value += 1
        ttk.Radiobutton(box, text="Tab", variable=self.delimiter_value, value=2).\
            grid(row=box_value, column=0, pady=3)
        box_value += 1
        ttk.Radiobutton(box, text="Comma", variable=self.delimiter_value, value=3).\
            grid(row=box_value, column=0, pady=3)
        row_value += 1
        box2 = ttk.Frame(win, relief="groove", padding=5)
        box2.grid(row=row_value, column=0, columnspan=2, pady=7)
        box2_value = 0
        ttk.Label(box2, text="File encoding:").grid(row=box2_value, column=0)
        box2_value += 1
        ttk.Radiobutton(box2, text="UTF-8", variable=self.file_encoding_value, value=1).\
            grid(row=box2_value, column=0, pady=5)
        box2_value += 1
        ttk.Radiobutton(box2, text="UTF-16", variable=self.file_encoding_value, value=2).\
            grid(row=box2_value, column=0, pady=5)
        box2_value = 0
        row_value += 1
        box3 = ttk.Frame(win, relief="groove", padding=5)
        box3.grid(row=row_value, column=0, columnspan=2, pady=7)
        box3_value = 0
        ttk.Label(box3, text="File name pattern:").grid(row=box3_value, column=1)
        box3_value += 1
        self.file_name_pattern_entry = ttk.Entry(box3, width=30, justify="center", font=LARGE_FONT)
        #self.file_name_pattern_entry.insert(tk.END, global_file_name_pattern)
        self.file_name_pattern_entry.insert(tk.END, "<H><1>_<F>Hz__<N>.txt")
        self.file_name_pattern_entry.grid(row=box3_value, column=1)
        box3_value += 1
        ttk.Label(box3, text="<H> for handle input").grid(row=box3_value, column=1)
        #box3_value += 1
        #ttk.Label(box3, text="<E> for electrode #").grid(row=box3_value, column=1)
        box3_value += 1
        ttk.Label(box3, text="<F> for frequency").grid(row=box3_value, column=1)
        box3_value += 1
        ttk.Label(box3, text="<N> for file #").grid(row=box3_value, column=1)
        #box3_value += 1
        #ttk.Label(box3, text="For Multichannel <E> is omitted").grid(row=box3_value, column=1)
        row_value += 1
        ttk.Button(win, text="Apply", command=self.apply_customize_file_format_settings).\
            grid(row=row_value, column=0, pady=6)
        ttk.Button(win, text="Done", command=win.destroy).\
            grid(row=row_value, column=1, pady=3)

    def apply_customize_file_format_settings(self) -> None:
        """Apply the user's settings for file format.
        """
        global global_current_column,\
            global_base_column_index_for_currents,\
            global_voltage_column,\
            global_voltage_column_index,\
            global_columns_per_electrode,\
            global_delimiter,\
            global_file_encoding,\
            global_file_name_pattern
        global_current_column = int(self.list_val_entry.get())
        global_base_column_index_for_currents = global_current_column - 1
        global_columns_per_electrode = int(self.spacing_val_entry.get())
        global_voltage_column = int(self.voltage_column.get())
        global_voltage_column_index = global_voltage_column - 1
        global_delimiter = Delimiter(self.delimiter_value.get())
        global_file_encoding = FileEncoding(self.file_encoding_value.get())
        global_file_name_pattern = self.file_name_pattern_entry.get()

    # def set_bytes(self, bytes_threshold: int, index: int) -> None:
    #     """Set the byte limit for filtering out preinitialized files.
    #         Here a preinitialized file is defined as a file that is currently being
    #         created by the instrument.
    #     """
    #     global global_byte_limit, global_byte_index
    #     #-- reset the self.byte_menu widgets --#
    #     self.byte_menu.entryconfigure(index, label=f"✓{bytes_threshold}")
    #     self.byte_menu.entryconfigure(global_byte_index, label=f"   {global_byte_limit}")
    #     #-- now change the current data being used --#
    #     global_byte_limit = bytes_threshold
    #     global_byte_index = index

    def show_frame(self, cont: str) -> None:
        """Raise the frame that is currently being displayed.
        """
        frame = global_show_frames[cont]
        frame.tkraise()

class InputFrame(ttk.Frame):
    """
    First frame that is displayed when the program is initialized, and
    contains all of the widgets for the user to input the settings for
    the experiment.
    """
    def __init__(self, parent: ttk.Frame, controller: tk.Tk):
        self.parent = parent
        self.controller = controller
        ttk.Frame.__init__(self, parent)             # initialize the frame
        row_value = 0
        ##############################################
        ### Pack all of the widgets into the frame ###
        ##############################################
        self.select_file_path: ttk.Button =\
            ttk.Button(self, style="Off.TButton", text="Select File Path",\
                       command=lambda: self.set_data_directory(parent))
        self.select_file_path.grid(row=row_value, column=0, columnspan=2)
        row_value += 2
        self.no_selected_path =\
            ttk.Label(self, text="No File Path Selected", font=MEDIUM_FONT, foreground="red")
        self.path_warning_exists = False
        ttk.Label(self, text="Import File Label", font=LARGE_FONT).\
            grid(row=row_value, column=0, columnspan=2)
        self.import_file_entry = ttk.Entry(self)
        self.import_file_entry.grid(row=row_value+1, column=0, columnspan=2, pady=5)
        self.import_file_entry.insert(tk.END, global_handle_variable)
        #--- File Handle Input ---#
        ttk.Label(self, text="Exported File Handle:", font=LARGE_FONT).\
            grid(row=row_value, column=2, columnspan=2)
        self.filehandle = ttk.Entry(self)
        now = datetime.datetime.now()
        self.filehandle.insert(tk.END,
                               f"DataExport_{now.year:4d}_{now.month:02d}_{now.day:02d}_" +\
                                f"{now.hour:02d}_{now.minute:02d}_{now.second:02d}.txt"
                                )
        self.filehandle.grid(row=row_value+1, column=2, columnspan=2, pady=5)
        row_value += 2
        ttk.Label(self, text="", font=LARGE_FONT).\
            grid(row=row_value, rowspan=2, column=0, columnspan=10)
        row_value += 1
        #---File Limit Input---#
        ttk.Label(self, text="Number of Files:", font=LARGE_FONT).\
            grid(row=row_value, column=0, columnspan=2, pady=4)
        self.numfiles = ttk.Entry(self, width=7)
        self.numfiles.insert(tk.END, "50")
        self.numfiles.grid(row=row_value+1, column=0, columnspan=2, pady=6)
        #--- Analysis interval for event callback in ElectrochemicalAnimation ---#
        ttk.Label(self, text="Analysis Interval (ms):", font=LARGE_FONT).\
            grid(row=row_value, column=2, columnspan=2, pady=4)
        self.analysis_interval = ttk.Entry(self, width=7)
        self.analysis_interval.insert(tk.END, "10")
        self.analysis_interval.grid(row=row_value+1, column=2, columnspan=2, pady=6)
        row_value += 2
        #--- Resize interval to update figures
        ttk.Label(self, text="Resize Interval", font=LARGE_FONT).\
            grid(row=row_value, column=2, columnspan=2)
        self.resize_entry = ttk.Entry(self, width=7)
        self.resize_entry.insert(tk.END, "200")
        self.resize_entry.grid(row=row_value+1, column=2, columnspan=2)
        row_value += 2
        ##################################
        ### Select and Edit Electrodes ###
        ##################################
        self.electrode_listbox_frame = ttk.Frame(self)
        self.electrode_listbox_frame.grid(row=row_value, column=0, columnspan=2,\
                                          padx=10, pady=10, ipady=5, sticky="nsew")
        #--- parameters for handling resize ---#
        self.electrode_listbox_frame.rowconfigure(0, weight=1)
        self.electrode_listbox_frame.rowconfigure(1, weight=1)
        self.electrode_listbox_frame.columnconfigure(0, weight=1)
        self.electrode_listbox_frame.columnconfigure(1, weight=1)
        self.electrode_list_exists = False
        self.electrode_label =\
            ttk.Label(self.electrode_listbox_frame, text="Select Electrodes:", font=LARGE_FONT)
        self.electrode_label.grid(row=0, column=0, columnspan=2, sticky="nswe")
        self.electrode_count =\
            tk.Listbox(self.electrode_listbox_frame, exportselection=0,\
                       width=10, font=LARGE_FONT, height=6,\
                        selectmode=tk.MULTIPLE, bd=3)
        self.electrode_count.bind("<<ListboxSelect>>", self.electrode_cur_select)
        self.electrode_count.grid(row=1, column=0, columnspan=2, sticky="nswe")
        for electrode in ELECTRODES_NUMBERS:
            self.electrode_count.insert(tk.END, electrode)
        self.scrollbar = ttk.Scrollbar(self.electrode_listbox_frame, orient="vertical",\
                                       command=self.electrode_count.yview)
        self.scrollbar.grid(row=1, column=1, sticky="nse")
        self.electrode_count.config(yscrollcommand=self.scrollbar.set)
        #--- Option to have data for all electrodes in a single file ---#
        self.single_electrode_file =\
            ttk.Button(self.electrode_listbox_frame, text="Multichannel", style="On.TButton",\
                       command=lambda: self.electrode_select(ElectrodesMode.SINGLE))
        self.single_electrode_file.grid(row=2, column=0,columnspan=2)
        # Disabling Multiplex selection since CH Instruments don't support 'each electrode in a separate file' format
        #--- Option to have data for each electrode in a separate file ---#
        #self.multiple_electrode_files =\
        #    ttk.Button(self.electrode_listbox_frame, text="Multiplex", style="Off.TButton",\
        #               command=lambda: self.electrode_select(ElectrodesMode.MULTIPLE))
        #self.multiple_electrode_files.grid(row=2, column=1)
        #--- Frame for editing electrodes ---#
        self.electrode_settings_frame = ttk.Frame(self, relief="groove")
        self.electrode_settings_frame.grid(row=10, column=0, columnspan=2,\
                                           padx=10, pady=10, sticky="nsew")
        self.electrode_settings_frame.columnconfigure(0, weight=1)
        self.electrode_settings_frame.rowconfigure(0, weight=1)
        self.electrode_settings_frame.rowconfigure(1, weight=1)
        self.electrode_settings_frame.rowconfigure(2, weight=1)
        #####################################################
        ### Select and Edit Frequencies for Data Analysis ###
        #####################################################
        # create a frame to pack in the frequency box and scrollbar:
        self.frequencies_listbox_frame = ttk.Frame(self)
        self.frequencies_listbox_frame.grid(row=row_value, column=2, columnspan=2,\
                                            padx=10, pady=10, sticky="nsew")
        frequencies = global_input_frequencies
        #-- resize ---#
        self.frequencies_listbox_frame.rowconfigure(0, weight=1)
        self.frequencies_listbox_frame.rowconfigure(1, weight=1)
        self.frequencies_listbox_frame.columnconfigure(0, weight=1)
        self.frequency_label = tk.Label(self.frequencies_listbox_frame,\
                                        text="Select Frequencies", font=LARGE_FONT)
        self.frequency_label.grid(row=0, padx=10)
        self.use_scroll_bar_for_frequencies = True
        #--- Variable to check if the frequency_list contains frequencies ---#
        self.frequency_list_exists = False
        self.frequency_list = tk.Listbox(self.frequencies_listbox_frame, relief="groove",\
                                         exportselection=0, width=5, font=LARGE_FONT, height=6,\
                                            selectmode=tk.MULTIPLE, bd=3)
        self.frequency_list.bind("<<ListboxSelect>>", self.frequency_cur_select)
        self.frequency_list.grid(row=1, padx=10, sticky="nswe")
        for frequency in frequencies:
            self.frequency_list.insert(tk.END, frequency)
        #--- Scroll Bar ---#
        self.scrollbar = ttk.Scrollbar(self.frequencies_listbox_frame, orient="vertical",\
                                       command=self.frequency_list.yview)
        #self.scrollbar.config(width=10, command=self.frequency_list.yview)
        self.scrollbar.grid(row=1, sticky="nse")
        self.frequency_list.config(yscrollcommand=self.scrollbar.set)
        ###########################################################
        ### Frame for adding/deleting frequencies from the list ###
        ###########################################################
        manipulate_frequencies_frame = ttk.Frame(self, width=10, relief="groove")
        ttk.Button(self.frequencies_listbox_frame, text="Edit",\
                  command=manipulate_frequencies_frame.tkraise).\
                    grid(row=2, column=0, columnspan=4)
        manipulate_frequencies_frame.grid(row=row_value, column=2, columnspan=2,\
                                          padx=10, pady=10, sticky="nsew")
        ttk.Label(manipulate_frequencies_frame, text="Enter Frequency(s)", font=MEDIUM_FONT).\
            grid(row=0, column=0, columnspan=4)
        self.frequency_entry: ttk.Entry = ttk.Entry(manipulate_frequencies_frame, width=8)
        self.frequency_entry.grid(row=1, column=0, columnspan=4)
        ttk.Button(manipulate_frequencies_frame, text="Add",\
                  command=self.add_frequency).\
                    grid(row=2, column=0)
        ttk.Button(manipulate_frequencies_frame, text="Delete",\
                  command=self.delete_frequency).\
                    grid(row=2, column=1)
        ttk.Button(manipulate_frequencies_frame, text="Clear",\
                  command=self.clear_all_frequencies).\
                    grid(row=3, column=0, columnspan=2)
        ttk.Button(manipulate_frequencies_frame, text="Return",\
                  command=self.return_from_frequencies_frame).\
                    grid(row=4, column=0, columnspan=2)
        manipulate_frequencies_frame.rowconfigure(0, weight=1)
        manipulate_frequencies_frame.rowconfigure(1, weight=1)
        manipulate_frequencies_frame.rowconfigure(2, weight=1)
        manipulate_frequencies_frame.rowconfigure(3, weight=1)
        manipulate_frequencies_frame.rowconfigure(4, weight=1)
        manipulate_frequencies_frame.columnconfigure(0, weight=1)
        manipulate_frequencies_frame.columnconfigure(1, weight=1)
        row_value += 1
        #--- Select KDM Method---#
        ttk.Label(self, text="Select KDM Method", font=LARGE_FONT).\
            grid(row=row_value, column=0, columnspan=2)
        self.kdm_box = tk.Listbox(self, relief="groove", exportselection=0,\
                                      font=LARGE_FONT, height=len(KDMMethod),\
                                        selectmode=tk.SINGLE, bd=3)
        self.kdm_box.bind("<<ListboxSelect>>", self.select_kdm_analysis_method)
        self.kdm_box.grid(row=row_value+1, column=0, columnspan=2)
        for method in KDMMethod:
            self.kdm_box.insert(tk.END, str(method))        
        #--- Select Analysis Method---#
        ttk.Label(self, text="Select Analysis Method", font=LARGE_FONT).\
            grid(row=row_value, column=2, columnspan=2)
        self.methods_box = tk.Listbox(self, relief="groove", exportselection=0,\
                                      font=LARGE_FONT, height=len(AnalysisMethod),\
                                        selectmode=tk.SINGLE, bd=3)
        self.methods_box.bind("<<ListboxSelect>>", self.select_data_analysis_method)
        self.methods_box.grid(row=row_value+1, column=2, columnspan=2)
        # Disabling Frequency Map Selection since its implementation is not final
        #for method in AnalysisMethod:
        #    self.methods_box.insert(tk.END, str(method)) 
        self.methods_box.insert(tk.END, str(AnalysisMethod.CONTINUOUS_SCAN))
        row_value += 2
        #--- Select Data to be Plotted ---#
        ttk.Label(self, text="Select Data to be Plotted", font=LARGE_FONT).\
            grid(row=row_value, column=0, columnspan=2)
        self.plot_options = tk.Listbox(self, relief="groove", exportselection=0,\
                                       font=LARGE_FONT, height=len(PlotSummaryMode),\
                                        selectmode=tk.SINGLE, bd=3)
        self.plot_options.bind("<<ListboxSelect>>", self.select_plot_summary_mode)
        self.plot_options.grid(row=row_value+1, column=0, columnspan=2)
        for summary_option in PlotSummaryMode:
            self.plot_options.insert(tk.END, str(summary_option))
        #--- Select Fitting Method---#
        ttk.Label(self, text="Select Fitting Method", font=LARGE_FONT).\
            grid(row=row_value, column=2, columnspan=2)
        self.peak_box = tk.Listbox(self, relief="groove", exportselection=0,\
                                      font=LARGE_FONT, height=len(PeakMethod),\
                                        selectmode=tk.SINGLE, bd=3)
        self.peak_box.bind("<<ListboxSelect>>", self.select_peak_analysis_method)
        self.peak_box.grid(row=row_value+1, column=2, columnspan=2)
        for method in PeakMethod:
            self.peak_box.insert(tk.END, str(method))
        row_value += 2
        #--- Select units of the X-axis ---#
        ttk.Label(self, text="Select X-axis units", font=LARGE_FONT).\
            grid(row=row_value, column=0, columnspan=4)
        self.x_axis_options = tk.Listbox(self, relief="groove", exportselection=0,\
                                         font=LARGE_FONT, height=len(PlotTimeReportingMode),\
                                            selectmode=tk.SINGLE, bd=3)
        self.x_axis_options.bind("<<ListboxSelect>>", self.select_xaxis_options)
        self.x_axis_options.grid(row=row_value+1, column=0, columnspan=4)
        for reporting_option in PlotTimeReportingMode:
            self.x_axis_options.insert(tk.END, str(reporting_option))
        row_value += 2
        ############################################################
        ### Adjustment of Visualization Parameters: xstart, xend ###
        ############################################################
        #--- Create a frame that will contain all of the widgets ---#
        adjustment_frame = ttk.Frame(self, relief="groove")
        adjustment_frame.grid(row=row_value, column=0, columnspan=4, pady=15)
        row_value += 1
        adjustment_frame.rowconfigure(0, weight=1)
        adjustment_frame.rowconfigure(1, weight=1)
        adjustment_frame.rowconfigure(2, weight=1)
        adjustment_frame.rowconfigure(3, weight=1)
        adjustment_frame.rowconfigure(4, weight=1)
        adjustment_frame.columnconfigure(0, weight=1)
        adjustment_frame.columnconfigure(1, weight=1)
        adjustment_frame.columnconfigure(2, weight=1)
        adjustment_frame.columnconfigure(3, weight=1)
        #--- Y Limit Adjustment Variables ---#
        self.y_limit_parameter_label =\
            ttk.Label(adjustment_frame, text="Select Y Limit Parameters", font=LARGE_FONT)
        self.y_limit_parameter_label.grid(row=0, column=0, columnspan=4, pady=5, padx=5)
        #--- Raw Data Minimum Parameter Adjustment ---#
        self.raw_data_min_parameter_label =\
            ttk.Label(adjustment_frame, text="Raw Min. Factor", font=MEDIUM_FONT)
        self.raw_data_min_parameter_label.grid(row=1, column=0)
        self.raw_data_min = ttk.Entry(adjustment_frame, width=5)
        # initial minimum is set to 0.5*minimum current (baseline) of file 1:
        self.raw_data_min.insert(tk.END, "0.5")
        self.raw_data_min.grid(row=2, column=0, padx=5, pady=2, ipadx=2)
        #--- Raw Data Maximum Parameter Adjustment ---#
        self.raw_data_max_parameter_label =\
            ttk.Label(adjustment_frame, text="Raw Max. Factor", font=MEDIUM_FONT)
        self.raw_data_max_parameter_label.grid(row=3, column=0)
        self.raw_data_max = ttk.Entry(adjustment_frame, width=5)
        # initial adjustment is set to 0.5x the max current (Peak Height) of file 1
        self.raw_data_max.insert(tk.END, "0.5")
        self.raw_data_max.grid(row=4, column=0, padx=5, pady=2, ipadx=2)
        #--- Raw Data Minimum Parameter Adjustment ---#
        self.data_min_parameter_label =\
            ttk.Label(adjustment_frame, text="Data Min. Factor", font=MEDIUM_FONT)
        self.data_min_parameter_label.grid(row=1, column=1)
        self.data_min = ttk.Entry(adjustment_frame, width=5)
        # initial minimum is set to 0.5*minimum current (baseline) of file 1
        self.data_min.insert(tk.END, "0.5")
        self.data_min.grid(row=2, column=1, padx=5, pady=2, ipadx=2)
        #--- Raw Data Maximum Parameter Adjustment ---#
        self.data_max_parameter_label =\
            ttk.Label(adjustment_frame, text="Data Max. Factor", font=MEDIUM_FONT)
        self.data_max_parameter_label.grid(row=3, column=1)
        self.data_max = ttk.Entry(adjustment_frame, width=5)
        # initial adjustment is set to 1.5x the max current (Peak Height) of file 1
        self.data_max.insert(tk.END, "1.5")
        self.data_max.grid(row=4, column=1, padx=5, pady=2, ipadx=2)
        #--- Normalized Data Minimum Parameter Adjustment ---#
        self.norm_data_min_parameter_label =\
            ttk.Label(adjustment_frame, text="Norm. Min.", font=MEDIUM_FONT)
        self.norm_data_min_parameter_label.grid(row=1, column=2)
        self.norm_data_min = ttk.Entry(adjustment_frame, width=5)
        self.norm_data_min.insert(tk.END, "0.5")
        self.norm_data_min.grid(row=2, column=2, padx=5, pady=2, ipadx=2)
        #--- Normalized Data Maximum Parameter Adjustment ---#
        self.norm_data_max_parameter_label =\
            ttk.Label(adjustment_frame, text="Norm. Max.", font=MEDIUM_FONT)
        self.norm_data_max_parameter_label.grid(row=3, column=2)
        self.norm_data_max = ttk.Entry(adjustment_frame, width=5)
        self.norm_data_max.insert(tk.END, "2")
        self.norm_data_max.grid(row=4, column=2, padx=5, pady=2, ipadx=2)
        #--- Raw Data Minimum Parameter Adjustment ---#
        self.kdm_min_label = ttk.Label(adjustment_frame, text="KDM Min.", font=MEDIUM_FONT)
        self.kdm_min_label.grid(row=1, column=3)
        self.kdm_min = ttk.Entry(adjustment_frame, width=5)
        # initial minimum is set to 0.5*minimum current (baseline) of file 1:
        self.kdm_min.insert(tk.END, "0.5")
        self.kdm_min.grid(row=2, column=3, padx=5, pady=2, ipadx=2)
        #--- Raw Data Maximum Parameter Adjustment ---#
        self.kdm_max_label = ttk.Label(adjustment_frame, text="KDM Max. ", font=MEDIUM_FONT)
        self.kdm_max_label.grid(row=3, column=3)
        self.kdm_max = ttk.Entry(adjustment_frame, width=5)
        # initial adjustment is set to 2x the max current (Peak Height) of file 1:
        self.kdm_max.insert(tk.END, "2")
        self.kdm_max.grid(row=4, column=3, padx=5, pady=2, ipadx=2)
        #--- Ask the User whether to export the data to a .txt file ---#
        self.export_data_flag = tk.BooleanVar()
        self.export_data_flag.set(False)
        ttk.Checkbutton(self, variable=self.export_data_flag, onvalue=True, offvalue=False,\
                       text="Export Data").\
                        grid(row=row_value, column=0, columnspan=2)
        #--- Ask the User if this is an injection experiment ---#
        self.injection_flag = tk.BooleanVar()
        self.injection_flag.set(False)
        ttk.Checkbutton(self, variable=self.injection_flag, onvalue=True, offvalue=False,\
                       text="Injection Experiment?").\
                        grid(row=row_value, column=2, columnspan=2)
        row_value += 1
        #--- Quit Button ---#
        ttk.Button(self, width=9, text="Quit",\
                   command=lambda: os._exit(0)).\
                    grid(row=row_value, column=0, columnspan=2, pady=10, padx=10)
        #--- Button to Initialize Data Analysis --#
        ttk.Button(self, width=9, text="Initialize",\
                   command=self.initialize_data_analysis).\
                    grid(row=row_value, column=2, columnspan=2, pady=10, padx=10)
        row_value += 1
        for row in range(row_value):
            row += 1
            self.rowconfigure(row, weight=1)
        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)
        self.columnconfigure(2, weight=1)
        self.columnconfigure(3, weight=1)
        ### Raise the initial frame for Electrode and Frequency Selection ###
        self.frequencies_listbox_frame.tkraise()
        self.electrode_listbox_frame.tkraise()

    #################################################
    ### Functions to track Selections and Entries ###
    #################################################
    def add_frequency(self) -> None:
        """Add the frequency or frequencies entered by the user to the list of frequencies."""
        frequencies_entered: Optional[str] = self.frequency_entry.get()
        self.frequency_entry.delete(0, tk.END)
        if frequencies_entered is not None:
            frequency_list: List[str] = frequencies_entered.split(" ")
            for frequency in frequency_list:
                if int(frequency) not in global_input_frequencies:
                    global_input_frequencies.append(int(frequency))
            global_input_frequencies.sort()
            self.frequency_list.delete(0, 1)
            self.frequency_list.delete(0, tk.END)
            for frequency_int in global_input_frequencies:
                self.frequency_list.insert(tk.END, frequency_int)

    def delete_frequency(self) -> None:
        """Delete the frequency or frequencies entered by the user from the list of frequencies."""
        frequencies_entered: Optional[str] = self.frequency_entry.get()
        self.frequency_entry.delete(0, tk.END)
        if frequencies_entered is not None:
            frequency_list: List[str] = frequencies_entered.split(" ")
            for frequency in frequency_list:
                frequency_int: int = int(frequency)
                if frequency_int in global_input_frequencies:
                    global_input_frequencies.remove(frequency_int)
                self.frequency_list.delete(0, tk.END)
                for frequency_int in global_input_frequencies:
                    self.frequency_list.insert(tk.END, frequency_int)

    def clear_all_frequencies(self) -> None:
        """Clear all of the frequencies from the list of frequencies."""
        global global_input_frequencies
        self.frequency_list.delete(0, tk.END)
        global_input_frequencies = []

    def return_from_frequencies_frame(self) -> None:
        """Abandon the frame for selecting frequencies and go back to the input frame."""
        self.frequencies_listbox_frame.tkraise()
        self.frequency_entry.delete(0, tk.END)

    def electrode_select(self, selector: ElectrodesMode) -> None:
        """Select whether the data for each electrode will be in a separate file or all electrode
        data in a single file.
        """
        global global_electrodes_mode
        match selector:
            case ElectrodesMode.SINGLE:
                global_electrodes_mode = ElectrodesMode.SINGLE
                self.single_electrode_file.config(style="On.TButton")
                self.multiple_electrode_files.config(style="Off.TButton")
            case ElectrodesMode.MULTIPLE:
                global_electrodes_mode = ElectrodesMode.MULTIPLE
                self.single_electrode_file.config(style="Off.TButton")
                self.multiple_electrode_files.config(style="On.TButton")

    def set_data_directory(self, parent) -> None:
        """Set the directory containing the files for data analysis.
        Input files should be placed there by the instrument control software.
        Output files will be placed in the same directory.
        """
        global global_file_path,\
            global_export_path,\
            global_found_file_path,\
            global_data_directory
        try:
            ### prompt the user to select a  ###
            ### directory for  data analysis ###
            dialog_outcome = filedialog.askdirectory(parent=parent)
            if not (dialog_outcome == () or dialog_outcome == ""):
                global_file_path = dialog_outcome + "/"
                ### Path for directory in which the    ###
                ### exported .txt file will be placed  ###
                local_export_path: List[str] = global_file_path.split("/")
                #-- change the text of the find file button to the folder the user chose --#
                global_data_directory = f"{local_export_path[-3]}/{local_export_path[-2]}"
                self.select_file_path.config(style="On.TButton")
                self.select_file_path.config(text=global_data_directory)
                del local_export_path[-1]
                del local_export_path[-1]
                global_export_path = "/".join(local_export_path) + "/"
                ## Indicates that the user has selected a File Path ###
                global_found_file_path = True
                if self.path_warning_exists:
                    self.no_selected_path.config(text="")
                    self.no_selected_path.grid_forget()
        except FileNotFoundError:
            print("find_file: Could Not Find File Path")

    #--- Analysis Method ---#
    def select_data_analysis_method(self, event: tk.Event) -> None:
        """Select the method for data analysis (continuous scan or frequency map)."""
        global global_analysis_method
        global_analysis_method = AnalysisMethod(self.methods_box.curselection()[0])

    def select_plot_summary_mode(self, event: tk.Event) -> None:
        """Select the plot summary mode (peak height or area under the curve)."""
        global global_plot_summary_mode
        #global_plot_summary_mode = self.plot_options.get(self.plot_options.curselection())
        global_plot_summary_mode = PlotSummaryMode(self.plot_options.curselection()[0])

    def select_xaxis_options(self, event: tk.Event) -> None:
        """Select the units for the x-axis (time or file number)."""
        global global_x_axis_mode
        #global_x_axis_mode = str((self.x_axis_options.get(self.x_axis_options.curselection())))
        global_x_axis_mode = PlotTimeReportingMode(self.x_axis_options.curselection()[0])

    #--- KDM Method ---#    
    def select_kdm_analysis_method(self, event: tk.Event) -> None:
        """Select the method for kdm analysis (old or new)."""
        global global_kdm_method
        global_kdm_method = KDMMethod(self.kdm_box.curselection()[0])

    #--- Peak Method ---#    
    def select_peak_analysis_method(self, event: tk.Event) -> None:
        """Select the method for peak analysis (poly fit or gauss)."""
        global global_peak_method
        global_peak_method = PeakMethod(self.peak_box.curselection()[0])

    #--- Electrode Selection ---#
    def electrode_cur_select(self, event: tk.Event) -> None:
        """Select the electrodes for data analysis."""
        global global_electrode_count,\
            global_electrode_list,\
            global_electrode_dict
        global_electrode_list = [int(self.electrode_count.get(idx))\
                                 for idx in self.electrode_count.curselection()]
        global_electrode_count = len(global_electrode_list)
        index = 0
        global_electrode_dict = {}
        for electrode in global_electrode_list:
            global_electrode_dict[electrode] = index
            index += 1
        if global_electrode_count == 0:
            self.electrode_list_exists = False
            print("electrode_cur_select: No Electrodes Selected")
            self.electrode_label.config(foreground="red")
        else:
            self.electrode_list_exists = True
            self.electrode_label.config(foreground="black")

    #--- Frequency Selection ---#
    def frequency_cur_select(self, event: tk.Event) -> None:
        """Select the frequencies for data analysis."""
        global global_frequency_list,\
            global_frequency_dict,\
            global_low_frequency,\
            global_high_frequency
        global_frequency_list = [self.frequency_list.get(idx)\
                                 for idx in self.frequency_list.curselection()]
        if len(global_frequency_list) == 0:
            self.frequency_list_exists = False
            self.frequency_label.config(foreground="red")
        else:
            self.frequency_list_exists = True
            self.frequency_label.config(foreground="black")
            # Initial Low Frequency for KDM/Ratiometric analysis:
            global_low_frequency = min(global_frequency_list)
            # Initial High Frequency for KDM/Ratiometric analysis:
            global_high_frequency = max(global_frequency_list)
            global_high_low_dictionary[HighLow.HIGH] = global_high_frequency
            global_high_low_dictionary[HighLow.LOW] = global_low_frequency
            #--- Frequency Dictionary ---#
            global_frequency_dict = {}
            count = 0
            for frequency in global_frequency_list:
                global_frequency_dict[frequency] = count
                count += 1

    def initialize_data_analysis(self) -> None:
        """Check to see if the user has filled out all required fields.
        If so, initialize the program.
        """
        global global_file_handle,\
            global_export_file_path
        #########################################################
        ### Initialize Canvases and begin tracking animation  ###
        #########################################################
        global_file_handle = self.filehandle.get() # handle for exported .txt file
        global_export_file_path = global_export_path + global_file_handle
        if global_export_path == "":
            self.no_selected_path.config(text="Please select a file path")
            self.no_selected_path.grid(row=0, column=2, columnspan=2)
            self.path_warning_exists = True
        else:
            if self.path_warning_exists:
                self.no_selected_path.grid_forget()
                self.path_warning_exists = False
        if self.frequency_list_exists:
            self.frequency_label.config(foreground="black")
        else:
            self.frequency_label.config(foreground="red")
        if self.electrode_list_exists:
            self.electrode_label.config(foreground="black")
        else:
            print("initialize_data_analysis: No Electrodes Selected")
            self.electrode_label.config(foreground="red")
        if not self.path_warning_exists:
            if self.frequency_list_exists:
                self.start_program()
            else:
                print("Could Not Start Program")

    def start_program(self) -> None:
        """Initialize the program and begin data acquisition, analysis, and animation."""
        global global_handle_variable,\
            global_search_interval,\
            global_resize_interval,\
            global_injection_point,\
            global_injection_selected,\
            global_min_kdm,\
            global_max_kdm,\
            global_min_norm,\
            global_max_norm,\
            global_min_raw,\
            global_max_raw,\
            global_min_data,\
            global_max_data,\
            global_high_frequency,\
            global_low_frequency,\
            global_initialized_normalization,\
            global_ratiometric_check,\
            global_normalization_vault,\
            global_text_file_export_activated,\
            global_number_of_files_to_process,\
            global_sample_rate,\
            global_normalization_point,\
            global_queue
        #---Get the User Input and make it globally accessible---#
        global_sample_rate = 20
        match global_analysis_method:
            case AnalysisMethod.CONTINUOUS_SCAN:
                global_number_of_files_to_process = int(self.numfiles.get())
            case AnalysisMethod.FREQUENCY_MAP:
                global_number_of_files_to_process = 1
        global_queue = Queue()
        global_injection_point = None
        # tracks if the data has been normalized to the starting normalization point:
        global_initialized_normalization = False
        # tracks changes to high and low frequencies:
        global_ratiometric_check = False
        global_normalization_point = DEFAULT_NORMALIZATION_FILE_NUMBER
        # tracks if text file export has been activated:
        global_text_file_export_activated = self.export_data_flag.get()
        #tracks if injection was selected:
        global_injection_selected = self.injection_flag.get()
        # interval at which xaxis of plots resizes
        global_resize_interval = int(self.resize_entry.get())
        # string handle used for the input file
        global_handle_variable = self.import_file_entry.get()
        #--- Y Limit Adjustment Parameters ---#
        global_min_norm = float(self.norm_data_min.get())
        global_max_norm = float(self.norm_data_max.get())
        global_min_raw = float(self.raw_data_min.get())
        global_max_raw = float(self.raw_data_max.get())
        global_min_data = float(self.data_min.get())
        global_max_data = float(self.data_max.get())
        global_min_kdm = float(self.kdm_min.get())
        global_max_kdm = float(self.kdm_max.get())
        #############################################################
        ### Interval at which the program searches for files (ms) ###
        #############################################################
        global_search_interval = int(self.analysis_interval.get())
        ## set the resizeability of the container ##
        ## frame to handle PlotContainer resize   ##
        global_container.columnconfigure(1, weight=1)
        #--- High and Low Frequency Selection for Drift Correction (KDM) ---#
        global_high_frequency = max(global_frequency_list)
        global_low_frequency = min(global_frequency_list)
        global_high_low_dictionary[HighLow.HIGH] = global_high_frequency
        global_high_low_dictionary[HighLow.LOW] = global_low_frequency
        #--- Create a timevault for normalization variables if the chosen normalization point
        #  has not yet been analyzed ---#
        # timevault for Normalization Points:
        global_normalization_vault = []
        # append the starting normalization point:
        global_normalization_vault.append(global_normalization_point)
        ################################################################
        ### If all checkpoints have been met, initialize the program ###
        ################################################################
        if global_found_file_path:
            CheckPoint(self.parent, self.controller)
#-------------------------------------------------------------------------------------------------#
class CheckPoint():
    """Check to see if the user's settings are accurate and the data files are present and valid.
    """
    def __init__(self, parent: ttk.Frame, controller: tk.Tk):
        #-- Check to see if the user's settings are accurate
        #-- Search for the presence of the files. If they exist,
        #-- initialize the functions and frames for Real Time Analysis
        self.win = tk.Toplevel(parent)
        self.win.wm_title("CheckPoint")    
        ttk.Label(self.win, text="Searching for files...", font=HUGE_FONT).\
            grid(row=0, column=0, columnspan=2, pady=10, padx=10, sticky="news")
        self.parent = parent
        self.win.transient()
        self.win.attributes("-topmost", "true")
        self.controller = controller
        row_value = 1
        self.frame_dict = {}
        self.already_verified: Dict[int, Dict[int, bool]] = {}
        for electrode in global_electrode_list:
            electrode_label = ttk.Label(self.win, text=f"{global_handle_variable}{electrode}", font=LARGE_FONT)
            electrode_label.grid(row=row_value, column=0, pady=5, padx=5)
            frame = ttk.Frame(self.win, relief="groove")
            frame.grid(row=row_value, column=1, pady=5, padx=5)
            self.frame_dict[electrode] = frame
            self.already_verified[electrode] = {}
            row_value += 1
            column_value = 0
            match global_analysis_method:
                case AnalysisMethod.CONTINUOUS_SCAN:
                    for frequency in global_frequency_list:
                        label = ttk.Label(frame, text=f"{frequency}Hz", foreground="red")
                        label.grid(row=0, column=column_value, padx=5, pady=5)
                        self.already_verified[electrode][frequency] = False
                        column_value += 1
                case AnalysisMethod.FREQUENCY_MAP:
                    electrode_label = ttk.Label(frame, text=f"E{electrode}", font=HUGE_FONT)
                    electrode_label.grid(row=row_value, column=column_value, pady=5, padx=5)
                    self.already_verified[electrode][global_frequency_list[0]] = False
                    if column_value == 1:
                        column_value = 0
                        row_value += 1
                    else:
                        column_value = 1
            if electrode != 1:
                frame.destroy()
                electrode_label.destroy()
        ttk.Button(self.win, text="Stop", command=self.stop).\
            grid(row=row_value, column=0, columnspan=2, pady=5)
        self.stop_search = False
        self.num = 0
        self.count = 0
        self.analysis_count = 0
        self.analysis_limit = global_electrode_count * len(global_frequency_list)
        self.electrode_limit = global_electrode_count - 1
        self.frequency_limit = len(global_frequency_list) - 1
        root.after(50, self.verify)
        #declarations for fields that will be initialized in other methods
        self.electrode: int

    def verify(self) -> None:
        """Checks the existence of the first data file (file index 1)
        for all electrodes and frequencies.
        If the files aren't all there, wait and retry.
        A file must be longer than a threshold - this seems to be done
        as a heuristic to avoid reading files currently
        being written by an instrument.
        The file contents will then be checked in the method verify_multi
        (but only for ELECTRODES_MODE_SINGLE).
        If all is well, i.e., there is at least the first set of data files,
        calls the method proceed to launch a window
        for monitoring the instrument.
        """
        self.electrode = global_electrode_list[self.num]
        if not self.stop_search:
            match global_analysis_method:
                case AnalysisMethod.CONTINUOUS_SCAN:
                    for frequency in global_frequency_list:
                        filename = make_file_name(1, self.electrode, frequency)
                        myfile = global_file_path + filename
                        print("FILENAME")
                        if file_is_complete(myfile):
                            print("INSIDE FILE IS COMPLETE IF")
                            match global_electrodes_mode:
                                case ElectrodesMode.SINGLE:
                                    check_ = self.verify_multi(myfile)
                                case ElectrodesMode.MULTIPLE:
                                    check_ = True
                            if check_:
                                if not self.already_verified[self.electrode][frequency]:
                                    self.already_verified[self.electrode][frequency] = True
                                    if not self.stop_search:
                                        self.analysis_count += 1
                            if self.analysis_count == self.analysis_limit:
                                if not self.stop_search:
                                    self.stop_search = True
                                    self.win.destroy()
                                    root.after(10, self.proceed)
                            print("SKIPPED FILE IS COMPLETE IF")        
                    if self.num < self.electrode_limit:
                        self.num += 1
                    else:
                        self.num = 0
                    if self.analysis_count < self.analysis_limit:
                        print("INSIDE STILL SEARCH")
                        if not self.stop_search:
                            root.after(100, self.verify)
                case AnalysisMethod.FREQUENCY_MAP:
                    frequency = global_frequency_list[0]
                    filename = make_file_name(1, self.electrode, frequency)
                    myfile = global_file_path + filename
                    if file_is_complete(myfile):
                        match global_electrodes_mode:
                            case ElectrodesMode.SINGLE:
                                check_ = self.verify_multi(myfile)
                            case ElectrodesMode.MULTIPLE:
                                check_ = True
                        if check_:
                            if not self.already_verified[self.electrode][frequency]:
                                self.already_verified[self.electrode][frequency] = True
                                if not self.stop_search:
                                    self.analysis_count += 1
                        if self.analysis_count == global_electrode_count:
                            if not self.stop_search:
                                self.stop_search = True
                                self.win.destroy()
                                root.after(10, self.proceed)
                    if self.num < self.electrode_limit:
                        self.num += 1
                    else:
                        self.num = 0
                    if self.analysis_count < self.analysis_limit:
                        if not self.stop_search:
                            root.after(200, self.verify)

    def verify_multi(self, myfile: str) -> bool:
        """This method ought to parse the file and ensure it is a
        valid instrument data file (from an
        instrument that records data from multiple electrodes in a
        single file, with several columns per
        electrode), such that it can subsequently by read with read_data.
        As it is, this method just looks for a line that has enough
        columns and starts with a number.
        """
        total_columns: int = 0
        # changing the column index
        #---Set the electrode index value---#
        check_ = False
        with open(myfile, "r", encoding=str(global_file_encoding)) as mydata:
            for line in mydata:
                # a comment here was "delete any tabs that may come before the
                # first value" but that is not what the code does
                check_split_list = line.split(str(global_delimiter))
                while check_split_list[0] == "" or check_split_list[0] == " ":
                    del check_split_list[0]
                check_split = check_split_list[0]
                check_split = check_split.replace(",", "")
                try:
                    _ = float(check_split)
                    check_split_is_float = True
                except ValueError:
                    check_split_is_float = False
                if check_split_is_float:
                    total_columns = len(check_split_list)
                    check_ = True
                    # note that here we declare all to be well even though
                    # we have only checked that the first item is a float whereas
                    # the rest of the line could be rubbish
                    ##observation: total_columns is the number of values in the
                    # first data line discovered in the file
                    break
        if check_:
            list_val = global_current_column + (self.electrode-1)*global_columns_per_electrode
            return list_val <= total_columns
        print("verify_multi: could not find a line that begins with a number")
        return False

    def proceed(self) -> None:
        """Actually initialize the program and begin data acquisition, analysis, and animation."""
        global global_wait_time,\
            global_track,\
            global_data_normalization,\
            global_post_analysis
        frame: ttk.Frame
        self.win.destroy()
        ##############################
        ### Synchronization Classes ###
        ##############################
        global_wait_time = WaitTime()
        global_track = Track()
        ######################################################
        ### Matplotlib Canvas, Figure, and Artist Creation ###
        ######################################################
        match global_analysis_method:
            case AnalysisMethod.CONTINUOUS_SCAN:
                InitializeContinuousCanvas()
                #################################
                ### Data Normalization Module ###
                #################################
                global_data_normalization = DataNormalization()
                ############################
                ### Post Analysis Module ###
                ############################
                global_post_analysis = PostAnalysis(self.parent, self.controller)
                global_show_frames[FRAME_POST_ANALYSIS] = global_post_analysis
                global_post_analysis.grid(row=0, column=0, sticky="nsew")
                ################################################
                ### Initialize the RealTimeManipulationFrame ###
                ################################################
                frame = ContinuousScanManipulationFrame(global_container)
                global_show_frames[str(global_analysis_method)] = frame
                frame.grid(row=0, column=0, sticky="nsew")
            case AnalysisMethod.FREQUENCY_MAP:
                InitializeFrequencyMapCanvas()
                ################################################
                ### Initialize the RealTimeManipulationFrame ###
                ################################################
                frame = FrequencyMapManipulationFrame(global_container)
                global_show_frames[str(global_analysis_method)] = frame
                frame.grid(row=0, column=0, sticky="nsew")
        #---When initialized, raise the Start Page and the plot for electrode 1---#
        # raises the frame for real-time data manipulation:
        self.show_frame(str(global_analysis_method))
        # raises the figure for electrode 1:
        self.show_plot(global_plot_values[0])

    def stop(self) -> None:
        """Callback for the Stop button in the user interface; stop searching for file 1,
        and allow the user to enter different settings, because presumably the current settings
        are wrong.
        """
        self.stop_search = True
        self.win.destroy()

    #--- Function to switch between visualization frames ---#
    def show_plot(self, frame: ttk.Frame) -> None:
        """Raise the given frame."""
        frame.tkraise()

    def show_frame(self, cont: str) -> None:
        """Raise the named frame."""
        frame = global_show_frames[cont]
        frame.tkraise()
#-------------------------------------------------------------------------------------------------#
############################################################
### Frame displayed during experiment with widgets and   ###
### functions for Real-time Data Manipulation            ###
############################################################
class ContinuousScanManipulationFrame(ttk.Frame):
    """Frame for real-time data manipulation during a continuous scan experiment.
    A single instance of this class is created whenever the user clicks Initialize.
    """
    def __init__(self, parent: ttk.Frame):
        global global_low_frequency_entry,\
            global_high_xstart_entry,\
            global_low_xstart_entry,\
            global_high_xend_entry,\
            global_low_xend_entry,\
            global_xstart_entry,\
            global_xend_entry,\
            global_high_frequency_entry,\
            global_norm_warning,\
            global_file_label,\
            global_real_time_sample_label,\
            global_set_point_norm,\
            global_normalization_var
        ttk.Frame.__init__(self, parent)
        #--- Display the file number ---#
        ttk.Label(self, text="File Number", font=MEDIUM_FONT,).\
            grid(row=0, column=0, padx=5, pady=5)
        global_file_label = ttk.Label(self, text="1", font=MEDIUM_FONT, style="Fun.TButton")
        global_file_label.grid(row=0, column=1, padx=5, pady=5)
        #--- Display the experiment duration as a function of the user-provided Sample Rate ---#
        ttk.Label(self, text="Experiment Time (h)", font=MEDIUM_FONT).\
            grid(row=0, column=2, padx=5, pady=5)
        global_real_time_sample_label =\
            ttk.Label(self, text="0", font=MEDIUM_FONT, style="Fun.TButton")
        global_real_time_sample_label.grid(row=0, column=3, padx=5, pady=5)
        #--- Real-time Normalization Variable ---#
        set_point_norm_label = ttk.Label(self, text="Set Normalization Point", font=MEDIUM_FONT)
        global_normalization_var = tk.StringVar()
        norm_string = str(DEFAULT_NORMALIZATION_FILE_NUMBER)
        global_normalization_var.set(norm_string)
        self.set_point_norm = ttk.Entry(self, textvariable=global_normalization_var, width=4)
        global_set_point_norm = self.set_point_norm
        #--- Button to apply any changes to the normalization variable ---#
        normalize_button = ttk.Button(self, text="Apply",\
                                     command=self.real_time_normalization, width=10)
        self.norm_warning = ttk.Label(self, text="", foreground="red", font=MEDIUM_FONT)
        global_norm_warning = self.norm_warning
        if global_injection_selected:
            set_point_norm_label.grid(row=2, column=0, pady=2, sticky="nsew")
            self.set_point_norm.grid(row=3, column=0, pady=2, padx=2)
            normalize_button.grid(row=4, column=0, pady=2, padx=2)
            self.norm_warning.grid(row=5, column=0, pady=0)
        else:
            set_point_norm_label.grid(row=2, column=0, columnspan=4, pady=2, sticky="nsew")
            self.set_point_norm.grid(row=3, column=0, columnspan=4, pady=2, padx=2)
            normalize_button.grid(row=4, column=0, columnspan=4, pady=2, padx=2)
            self.norm_warning.grid(row=5, column=0, columnspan=4, pady=0)
        #--- Real-time Injection tracking ---#
        set_injection_label = ttk.Label(self, text="Set Injection Range", font=MEDIUM_FONT)
        injection_button = ttk.Button(self, text="Apply",\
                                     command=self.real_time_injection, width=10)
        self.set_injection_point = ttk.Entry(self, width=8)
        ## If this is an injection experiment, grid the widgets ##
        if global_injection_selected:
            self.set_injection_point.grid(row=3, column=1, pady=2, padx=5)
            injection_button.grid(row=4, column=1, pady=2, padx=2)
            set_injection_label.grid(row=2, column=1, pady=2, sticky="nsew")
        row_value = 6
        if len(global_frequency_list) > 1:
            self.frequency_frame = ttk.Frame(self, relief="groove")
            self.frequency_frame.grid(row=row_value, column=0, columnspan=4,\
                                      pady=2, padx=3, ipady=2)
            #--- Drift Correction Title ---#
            ttk.Label(self.frequency_frame, text="Drift Correction", font=MEDIUM_FONT).\
                grid(row=0, column=0, columnspan=3, pady=1, padx=5)
            #--- High Frequency Selection for KDM and Ratiometric Analysis ---#
            ttk.Label(self.frequency_frame, text="High Frequency", font=MEDIUM_FONT).\
                grid(row=1, column=1, pady=3, padx=5)
            global_high_frequency_entry = ttk.Entry(self.frequency_frame, width=7)
            global_high_frequency_entry.insert(tk.END, str(global_high_frequency))
            global_high_frequency_entry.grid(row=2, column=1, padx=5)
            #--- Low Frequency Selection for KDM and Ratiometric Analysis ---#
            ttk.Label(self.frequency_frame, text="Low Frequency", font=MEDIUM_FONT).\
                grid(row=1, column=0, pady=3, padx=5)
            global_low_frequency_entry = ttk.Entry(self.frequency_frame, width=7)
            global_low_frequency_entry.insert(tk.END, str(global_low_frequency))
            global_low_frequency_entry.grid(row=2, column=0, padx=5)
            ttk.Label(self.frequency_frame, text="Low Frequency Offset", font=MEDIUM_FONT).\
                grid(row=3, column=0, pady=2, padx=2)
            self.low_frequency_offset = ttk.Entry(self.frequency_frame, width=7)
            self.low_frequency_offset.insert(tk.END, str(global_low_frequency_offset))
            self.low_frequency_offset.grid(row=4, column=0, padx=2, pady=2)
            ttk.Label(self.frequency_frame,\
                     text="Low Frequency Slope Manipulation", font=MEDIUM_FONT).\
                        grid(row=3, column=1, pady=2, padx=2)
            self.low_frequency_slope = ttk.Entry(self.frequency_frame, width=7)
            self.low_frequency_slope.insert(tk.END, str(global_low_frequency_slope))
            self.low_frequency_slope.grid(row=4, column=1, padx=2, pady=2)
            ttk.Button(self.frequency_frame, text="Apply",\
                      command=self.real_time_kdm).\
                        grid(row=5, column=0, columnspan=4, pady=3, padx=5)
            row_value += 1
        #################################################
        ### Nested Frame for Real-Time adjustment     ###
        ### of voltammogram and polynomial regression ###
        ### or Gauss method                           ###
        #################################################
        regression_frame = ttk.Frame(self, relief="groove", padding=5)
        regression_frame.grid(row=row_value, column=0, columnspan=4,\
                              pady=5, padx=5, ipadx=3, sticky="ns")
        regression_frame.rowconfigure(0, weight=1)
        regression_frame.rowconfigure(1, weight=1)
        regression_frame.rowconfigure(2, weight=1)
        regression_frame.columnconfigure(0, weight=1)
        regression_frame.columnconfigure(1, weight=1)
        row_value += 1

        match global_peak_method:
            case PeakMethod.POLY:
                #--- Title ---#
                ttk.Label(regression_frame,\
                        text="Poly Fit Savitzky-Golay smoothing window (mV)", font=MEDIUM_FONT).\
                            grid(row=0, column=0, columnspan=4, pady=5, padx=5)
                ###################################################################
                ### Real Time Manipulation of Savitzky-Golay Smoothing Function ###
                ###################################################################
                self.smoothing_entry = ttk.Entry(regression_frame, width=10)
                self.smoothing_entry.grid(row=2, column=0, columnspan=4, pady=3)
                self.smoothing_entry.insert(tk.END, str(global_savitzky_golay_window))
                
                parameter_frame = [([ttk.Frame]*len(global_frequency_list)) for i in range(global_electrode_count)]
                self.xstart_entry = [([ttk.Entry]*len(global_frequency_list)) for i in range(global_electrode_count)]
                self.xend_entry = [([ttk.Entry]*len(global_frequency_list)) for i in range(global_electrode_count)]

                for elec in range(global_electrode_count):
                    for freq in range(len(global_frequency_list)):
                        parameter_frame[elec][freq] = ttk.Frame(regression_frame, padding=5)
                        parameter_frame[elec][freq].grid(row=3, column=0, columnspan=6, sticky="nsew")
                        parameter_frame[elec][freq].rowconfigure(0, weight=1)
                        parameter_frame[elec][freq].rowconfigure(1, weight=1)
                        parameter_frame[elec][freq].rowconfigure(2, weight=1)
                        parameter_frame[elec][freq].columnconfigure(0, weight=1)
                        parameter_frame[elec][freq].columnconfigure(1, weight=1)
                        global_show_frames[f"{elec}{freq}"] = parameter_frame[elec][freq]
                        #--- points discarded at the beginning of the voltammogram, xstart ---#
                        ttk.Label(parameter_frame[elec][freq], text="xstart (V)", font=MEDIUM_FONT).\
                            grid(row=0, column=0)
                        self.xstart_entry[elec][freq] = ttk.Entry(parameter_frame[elec][freq], width=7)
                        self.xstart_entry[elec][freq].insert(tk.END, str(global_xstart[elec][freq]))
                        self.xstart_entry[elec][freq].grid(row=1, column=0)
                        global_xstart_entry[elec][freq] = self.xstart_entry[elec][freq]
                        #--- points discarded at the beginning of the voltammogram, xend ---#
                        ttk.Label(parameter_frame[elec][freq], text="xend (V)", font=MEDIUM_FONT).\
                            grid(row=0, column=max(1,len(global_frequency_list)-1))
                        self.xend_entry[elec][freq] = ttk.Entry(parameter_frame[elec][freq], width=7)
                        self.xend_entry[elec][freq].insert(tk.END, str(global_xend[elec][freq]))
                        self.xend_entry[elec][freq].grid(row=1, column=max(1,len(global_frequency_list)-1))
                        global_xend_entry[elec][freq] = self.xend_entry[elec][freq]

                        col = 0
                        for i in range(len(global_frequency_list)):
                            if i == freq:
                                styl = "On.TButton"
                            else:
                                styl = "Off.TButton"
                            frame_str = str(elec)+str(i)
                            ttk.Button(parameter_frame[elec][freq], style=styl,\
                                    text=f"{global_frequency_list[i]}Hz-E{global_electrode_list[elec]}",\
                                        command=lambda frame_str=frame_str: self.show_frame(frame_str)).\
                                            grid(row=4, column=col)
                            col += 1
                self.show_frame("00")
                #--- Button to apply adjustments ---#
                self.adjust_parameter_button = ttk.Button(regression_frame, text="Apply",\
                                                        command=self.adjust_parameters)
                self.adjust_parameter_button.grid(row=5, column=0, columnspan=4, pady=10, padx=10)
            case PeakMethod.GAUSS:
                ttk.Label(regression_frame, text="Gauss Method", font=MEDIUM_FONT).\
                        grid(row=0, column=0, columnspan=4, pady=5, padx=5)
                
                parameter_frame = [([ttk.Frame]*len(global_frequency_list)) for i in range(global_electrode_count)]
                self.gauss_peak_entry = [([ttk.Entry]*len(global_frequency_list)) for i in range(global_electrode_count)]
                self.gauss_baseline_entry = [([ttk.Entry]*len(global_frequency_list)) for i in range(global_electrode_count)]
                self.gauss_maxheight_entry = [([ttk.Entry]*len(global_frequency_list)) for i in range(global_electrode_count)]

                for elec in range(global_electrode_count):
                    for freq in range(len(global_frequency_list)):
                        parameter_frame[elec][freq] = ttk.Frame(regression_frame, padding=5)
                        parameter_frame[elec][freq].grid(row=3, column=0, columnspan=6, sticky="nsew")
                        parameter_frame[elec][freq].rowconfigure(0, weight=1)
                        parameter_frame[elec][freq].rowconfigure(1, weight=1)
                        parameter_frame[elec][freq].rowconfigure(2, weight=1)
                        parameter_frame[elec][freq].columnconfigure(0, weight=1)
                        parameter_frame[elec][freq].columnconfigure(1, weight=1)
                        global_show_frames[f"{elec}{freq}"] = parameter_frame[elec][freq]
                        #--- peak location of the voltammogram ---#
                        ttk.Label(parameter_frame[elec][freq], text="Exp. Peak Location (V)", font=MEDIUM_FONT).\
                            grid(row=0, column=0)
                        self.gauss_peak_entry[elec][freq] = ttk.Entry(parameter_frame[elec][freq], width=7)
                        self.gauss_peak_entry[elec][freq].insert(tk.END, str(global_gauss_peak[elec][freq]))
                        self.gauss_peak_entry[elec][freq].grid(row=1, column=0)
                        #--- baseline of the voltammogram ---#
                        ttk.Label(parameter_frame[elec][freq], text="Exp. Baseline (uA)", font=MEDIUM_FONT).\
                            grid(row=0, column=1)
                        self.gauss_baseline_entry[elec][freq] = ttk.Entry(parameter_frame[elec][freq], width=7)
                        self.gauss_baseline_entry[elec][freq].insert(tk.END, str(global_gauss_baseline[elec][freq]))
                        self.gauss_baseline_entry[elec][freq].grid(row=1, column=1)
                        #--- max peak height of the voltammogram ---#
                        ttk.Label(parameter_frame[elec][freq], text="Exp. Max Peak Height (uA)", font=MEDIUM_FONT).\
                            grid(row=0, column=2)
                        self.gauss_maxheight_entry[elec][freq] = ttk.Entry(parameter_frame[elec][freq], width=7)
                        self.gauss_maxheight_entry[elec][freq].insert(tk.END, str(global_gauss_maxheight[elec][freq]))
                        self.gauss_maxheight_entry[elec][freq].grid(row=1, column=2)

                        col = 0
                        for i in range(len(global_frequency_list)):
                            if i == freq:
                                styl = "On.TButton"
                            else:
                                styl = "Off.TButton"
                            frame_str = str(elec)+str(i)
                            ttk.Button(parameter_frame[elec][freq], style=styl,\
                                    text=f"{global_frequency_list[i]}Hz-E{global_electrode_list[elec]}",\
                                        command=lambda frame_str=frame_str: self.show_frame(frame_str)).\
                                            grid(row=4, column=col)
                            col += 1
                self.show_frame("00")
                #--- Button to apply adjustments ---#
                self.adjust_parameter_button = ttk.Button(regression_frame, text="Apply",\
                                                        command=self.adjust_parameters)
                self.adjust_parameter_button.grid(row=5, column=0, columnspan=4, pady=10, padx=10)
                
        self.x_axis_control_left = tk.IntVar()
        self.x_axis_control_left.set(XBound.START_PLUS.value)
        self.x_axis_control_right = tk.IntVar()
        self.x_axis_control_right.set(XBound.NOW_MINUS.value)
        x_axis_control_frame = ttk.Frame(self, relief="groove", padding=5)
        x_axis_control_frame.grid(row=row_value, column=0, columnspan=2,\
                                  pady=3, padx=3, ipadx=3, sticky="ns")
        x_axis_control_frame.rowconfigure(0, weight=1)
        x_axis_control_frame.rowconfigure(1, weight=1)
        x_axis_control_frame.rowconfigure(2, weight=1)
        x_axis_control_frame.rowconfigure(3, weight=1)
        x_axis_control_frame.rowconfigure(4, weight=1)
        x_axis_control_frame.columnconfigure(0, weight=1)
        x_axis_control_frame.columnconfigure(1, weight=1)
        x_axis_control_frame.columnconfigure(2, weight=1)
        x_axis_control_frame.columnconfigure(3, weight=1)
        ttk.Label(x_axis_control_frame, text="x-axis display range", font=MEDIUM_FONT).\
            grid(row=0, column=0, columnspan=4, pady=3, padx=5)
        ttk.Label(x_axis_control_frame, text="from:", font=MEDIUM_FONT).\
            grid(row=1, column=0, columnspan=2, pady=3, padx=5)
        ttk.Label(x_axis_control_frame, text="to:", font=MEDIUM_FONT).\
            grid(row=1, column=2, columnspan=2, pady=3, padx=5)
        ttk.Radiobutton(x_axis_control_frame, text="start +",\
                       variable=self.x_axis_control_left, value=XBound.START_PLUS.value).\
                        grid(row=2, column=0)
        self.entry_for_x_axis_control_left_start_plus = ttk.Entry(x_axis_control_frame, width=3)
        self.entry_for_x_axis_control_left_start_plus.grid(row=2, column=1)
        ttk.Radiobutton(x_axis_control_frame, text="now \u2014",\
                       variable=self.x_axis_control_left, value=XBound.NOW_MINUS.value).\
                        grid(row=3, column=0)
        self.entry_for_x_axis_control_left_now_minus = ttk.Entry(x_axis_control_frame, width=3)
        self.entry_for_x_axis_control_left_now_minus.grid(row=3, column=1)
        ttk.Radiobutton(x_axis_control_frame, text="start +",\
                       variable=self.x_axis_control_right, value=XBound.START_PLUS.value).\
                        grid(row=2, column=2)
        self.entry_for_x_axis_control_right_start_plus = ttk.Entry(x_axis_control_frame, width=3)
        self.entry_for_x_axis_control_right_start_plus.grid(row=2, column=3)
        ttk.Radiobutton(x_axis_control_frame, text="now \u2014",\
                       variable=self.x_axis_control_right, value=XBound.NOW_MINUS.value).\
                        grid(row=3, column=2)
        self.entry_for_x_axis_control_right_now_minus = ttk.Entry(x_axis_control_frame, width=3)
        self.entry_for_x_axis_control_right_now_minus.grid(row=3, column=3)
        ttk.Button(x_axis_control_frame, text="Apply", width=6,\
                  command=self.apply_x_axis_settings).\
                    grid(row=4, column=1)
        self.y_axis_control_norm = tk.IntVar()
        self.y_axis_control_norm.set(global_y_norm_radiobutton.value)
        self.y_axis_control_kdm = tk.IntVar()
        self.y_axis_control_kdm.set(global_y_kdm_radiobutton.value)
        y_axis_control_frame = ttk.Frame(self, relief="groove", padding=5)
        y_axis_control_frame.grid(row=row_value, column=2, columnspan=2,\
                                  pady=3, padx=3, ipadx=3, sticky="ns")
        y_axis_control_frame.rowconfigure(0, weight=1)
        y_axis_control_frame.rowconfigure(1, weight=1)
        y_axis_control_frame.rowconfigure(2, weight=1)
        y_axis_control_frame.rowconfigure(3, weight=1)
        y_axis_control_frame.rowconfigure(4, weight=1)
        y_axis_control_frame.rowconfigure(5, weight=1)
        y_axis_control_frame.rowconfigure(6, weight=1)
        y_axis_control_frame.columnconfigure(0, weight=1)
        y_axis_control_frame.columnconfigure(1, weight=2)
        y_axis_control_frame.columnconfigure(2, weight=1)
        y_axis_control_frame.columnconfigure(3, weight=2)
        ttk.Label(y_axis_control_frame, text="y-axis display range", font=MEDIUM_FONT).\
            grid(row=0, column=1, columnspan=2)
        ttk.Label(y_axis_control_frame, text="Norm. Ratio", font=MEDIUM_FONT).\
            grid(row=1, column=0, columnspan=2, pady=3, padx=1)
        ttk.Label(y_axis_control_frame, text="KDM", font=MEDIUM_FONT).\
            grid(row=1, column=2, columnspan=2, pady=3, padx=1)
        ttk.Radiobutton(y_axis_control_frame, text="auto",\
                       variable=self.y_axis_control_norm, value=YBound.AUTOMATIC.value).\
                        grid(row=2, column=0)
        ttk.Radiobutton(y_axis_control_frame, text="manual",\
                       variable=self.y_axis_control_norm, value=YBound.MANUAL.value).\
                        grid(row=3, column=0)
        ttk.Label(y_axis_control_frame, text="min", font=MEDIUM_FONT).\
            grid(row=4, column=0)
        self.entry_for_y_axis_control_norm_manual_min = ttk.Entry(y_axis_control_frame, width=3)
        self.entry_for_y_axis_control_norm_manual_min.insert(tk.END, str(global_min_norm))
        self.entry_for_y_axis_control_norm_manual_min.grid(row=4, column=1)
        ttk.Label(y_axis_control_frame, text="max", font=MEDIUM_FONT).\
            grid(row=5, column=0)
        self.entry_for_y_axis_control_norm_manual_max = ttk.Entry(y_axis_control_frame, width=3)
        self.entry_for_y_axis_control_norm_manual_max.insert(tk.END, str(global_max_norm))
        self.entry_for_y_axis_control_norm_manual_max.grid(row=5, column=1)
        ttk.Radiobutton(y_axis_control_frame, text="auto",\
                       variable=self.y_axis_control_kdm, value=YBound.AUTOMATIC.value).\
                        grid(row=2, column=2)
        ttk.Radiobutton(y_axis_control_frame, text="manual",\
                       variable=self.y_axis_control_kdm, value=YBound.MANUAL.value).\
                        grid(row=3, column=2)
        ttk.Label(y_axis_control_frame, text="min", font=MEDIUM_FONT).\
            grid(row=4, column=2)
        self.entry_for_y_axis_control_kdm_manual_min = ttk.Entry(y_axis_control_frame, width=3)
        self.entry_for_y_axis_control_kdm_manual_min.insert(tk.END, str(global_min_kdm))
        self.entry_for_y_axis_control_kdm_manual_min.grid(row=4, column=3)
        ttk.Label(y_axis_control_frame, text="max", font=MEDIUM_FONT).\
            grid(row=5, column=2)
        self.entry_for_y_axis_control_kdm_manual_max = ttk.Entry(y_axis_control_frame, width=3)
        self.entry_for_y_axis_control_kdm_manual_max.insert(tk.END, str(global_max_kdm))
        self.entry_for_y_axis_control_kdm_manual_max.grid(row=5, column=3)
        ttk.Button(y_axis_control_frame, text="Apply", width=6,\
                  command=self.apply_y_axis_settings).\
                    grid(row=6, column=1)
        row_value += 1
        #---Buttons to switch between electrode frames---#
        frame_value = 0
        column_value = 0
        for _ in global_plot_values:
            ttk.Button(self, text=global_frame_list[frame_value],\
                       command=lambda frame_value=frame_value:\
                        self.show_plot(global_plot_values[frame_value],frame_value)).\
                            grid(row=row_value, column=column_value, pady=2, padx=5)
            ## allows .grid() to alternate between
            ## packing into column 1 and column 2
            if column_value == 1:
                column_value = 0
                row_value += 1
            ## if gridding into the 1st column,
            ## grid the next into the 2nd column
            else:
                column_value += 1
            frame_value += 1
        row_value += 1
        #--- Start ---#
        ttk.Button(self, text="Start", style="Fun.TButton", command=self.skeleton_key).\
            grid(row=row_value, column=0, pady=5, padx=5)
        #--- Reset ---#
        ttk.Button(self, text="Reset", style="Fun.TButton", command=self.reset).\
            grid(row=row_value, column=1, pady=5, padx=5)
        #row_value += 1
        #--- Quit ---#
        ttk.Button(self, text="Quit", command=lambda: os._exit(0)).\
            grid(row=row_value, column=3, columnspan=4, pady=5)
        for row in range(row_value):
            row += 1
            self.rowconfigure(row, weight=1)
        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)
        self.columnconfigure(2, weight=1)
        self.columnconfigure(3, weight=1)
                                        ###################################################
                                        ###################################################
                                        ###### Real Time Data Manipulation Functions ######
                                        ###################################################
                                        ###################################################
    #####################################
    ### Manipulation of the Injection ###
    ### Point for visualization       ###
    #####################################
    def real_time_injection(self) -> None:
        """Adjust the injection point in real time."""
        global global_injection_point
        global_injection_point = int(self.set_injection_point.get())
        print(f"New Injection Point: {global_injection_point}")

    def adjust_parameters(self) -> None:
        """Adjust the parameters for the polynomial regression."""
        global global_low_xstart,\
            global_high_xstart,\
            global_low_xend,\
            global_high_xend,\
            global_xstart,\
            global_xend,\
            global_savitzky_golay_window,\
            global_low_gauss_peak,\
            global_high_gauss_peak,\
            global_gauss_peak,\
            global_low_gauss_baseline,\
            global_high_gauss_baseline,\
            global_gauss_baseline,\
            global_low_gauss_maxheight,\
            global_high_gauss_maxheight,\
            global_gauss_maxheight
        match global_peak_method:
            case PeakMethod.POLY:
                ###############################################
                ### Polynomial Regression Range Parameters ###
                ###############################################
                for elec in range(global_electrode_count):
                    for freq in range(len(global_frequency_list)):
                        global_xstart[elec][freq] = float(self.xstart_entry[elec][freq].get())
                        global_xend[elec][freq] = float(self.xend_entry[elec][freq].get())
                #######################################
                ### Savitzky-Golay Smoothing Window ###
                #######################################
                global_savitzky_golay_window = float(self.smoothing_entry.get())
            case PeakMethod.GAUSS:
                for elec in range(global_electrode_count):
                    for freq in range(len(global_frequency_list)):
                        global_gauss_peak[elec][freq] = float(self.gauss_peak_entry[elec][freq].get())
                        global_gauss_baseline[elec][freq] = float(self.gauss_baseline_entry[elec][freq].get())
                        global_gauss_maxheight[elec][freq] = float(self.gauss_maxheight_entry[elec][freq].get())

    #########################################################
    ### Real-time adjustment of High and Low frequencies  ###
    ### used for KDM and ratiometric analysis             ###
    #########################################################
    def real_time_kdm(self) -> None:
        """Adjust the high and low frequencies used for KDM and ratiometric analysis."""
        global global_low_frequency_offset,\
            global_low_frequency_slope,\
            global_frequency_warning_label_exists,\
            global_wrong_frequency_label,\
            global_ratiometric_check
        temp_high_frequency = int(global_high_frequency_entry.get())
        temp_low_frequency = int(global_low_frequency_entry.get())
        global_low_frequency_offset = float(self.low_frequency_offset.get())
        global_low_frequency_slope = float(self.low_frequency_slope.get())
        #--- Reset the variable for the Warning Label (WrongFrequencyLabel) ---#
        if temp_low_frequency in global_frequency_list\
            and temp_high_frequency not in global_frequency_list:
            if global_frequency_warning_label_exists:
                global_wrong_frequency_label.grid_forget()
            global_wrong_frequency_label =\
                ttk.Label(self.frequency_frame, text="High Frequency Does Not Exist",\
                          foreground="red")
            global_wrong_frequency_label.grid(row=6, column=0, columnspan=4)
            if not global_frequency_warning_label_exists:
                global_frequency_warning_label_exists = True
        #--- if only the low_frequency does not exist ---#
        elif temp_low_frequency not in global_frequency_list\
            and temp_high_frequency in global_frequency_list:
            if global_frequency_warning_label_exists:
                global_wrong_frequency_label.grid_forget()
            global_wrong_frequency_label =\
                ttk.Label(self.frequency_frame, text="Low Frequency Does Not Exist",\
                          foreground="red")
            global_wrong_frequency_label.grid(row=6, column=0, columnspan=4)
            if not global_frequency_warning_label_exists:
                global_frequency_warning_label_exists = True
        #--- if both the high_frequency and low_frequency do not exist ---#
        elif temp_low_frequency not in global_frequency_list\
            and temp_high_frequency not in global_frequency_list:
            if global_frequency_warning_label_exists:
                global_wrong_frequency_label.grid_forget()
            global_wrong_frequency_label =\
                  ttk.Label(self.frequency_frame,\
                           text="High and Low Frequencies Do Not Exist", foreground="red")
            global_wrong_frequency_label.grid(row=6, column=0, columnspan=4)
            if not global_frequency_warning_label_exists:
                global_frequency_warning_label_exists = True
        #--- else, if they both exist, remove the warning label ---#
        else:
            global_high_low_dictionary[HighLow.HIGH] = temp_high_frequency
            global_high_low_dictionary[HighLow.LOW] = temp_low_frequency
            global_data_normalization.reset_ratiometric_data()
            #--- if a warning label exists, forget it ---#
            if global_frequency_warning_label_exists:
                global_wrong_frequency_label.grid_forget()
            #--- Tells RawVoltammogramVisualization to revisualize data
            # for new High and Low frequencies ---#
            global_ratiometric_check = True
            if global_analysis_complete:
                global_post_analysis._adjust_data()

    #--- Function for Real-time Normalization ---#
    def real_time_normalization(self) -> None:
        """Adjust the normalization point in real time."""
        global global_normalization_point
        global_normalization_point = int(self.set_point_norm.get())
        file = int(global_file_label.cget("text"))
        if file >= global_normalization_point:
            global_wait_time.normalization_wait_time()
        else:
            global_norm_warning.config(foreground="red")
            global_norm_warning.config(text=f"File {global_normalization_point}" +\
                                       f" has not been analyzed yet")
        if global_analysis_complete:
            global_post_analysis._adjust_data()

    def reset(self) -> None:
        """Reset the experiment."""
        global global_key,\
            global_poison_pill,\
            global_analysis_complete,\
            global_analysis_already_initiated,\
            global_low_already_reset,\
            global_high_already_reset,\
            global_already_reset,\
            global_x_left_bound_radiobutton,\
            global_x_left_bound_offset,\
            global_x_right_bound_radiobutton,\
            global_x_right_bound_offset,\
            global_gauss_solver
        global_key = 0
        global_poison_pill = True
        global_analysis_already_initiated = False # reset the start variable
        global_already_reset = True
        # Raise the initial user input frame
        self.show_frame(FRAME_INPUT)
        self.close_frame(str(global_analysis_method))
        global_post_analysis._reset()
        ## Take resize weight away from the Visualization Canvas
        global_container.columnconfigure(1, weight=0)
        global_analysis_complete = False
        global_x_left_bound_radiobutton = XBound.START_PLUS
        global_x_left_bound_offset = 0
        global_x_right_bound_radiobutton = XBound.NOW_MINUS
        global_x_right_bound_offset = 0

    def show_frame(self, cont: str) -> None:
        """Raise the frame to the front of the canvas."""
        frame = global_show_frames[cont]
        frame.tkraise()

    def skeleton_key(self) -> None:
        """Start the data analysis and visualization."""
        global global_key, global_poison_pill, global_analysis_already_initiated
        if global_analysis_already_initiated:
            print("Program has already been initiated")
        else:
            ######################################################################
            ### Initialize Animation (Visualization) for each electrode figure ###
            ######################################################################
            fig_count = 0                   # index value for the frame
            for figure in global_figures:
                fig, _ = figure
                electrode = global_electrode_list[fig_count]
                global_animations.append(\
                    ElectrochemicalAnimation(fig,\
                                            electrode,\
                                            resize_interval=global_resize_interval,\
                                            fargs=None))
                fig_count += 1
            global_analysis_already_initiated = True
            #--- reset poison pill variables --#
            global_poison_pill = False
            # tells Generate() to start data analysis
            if global_key == 0:
                global_key += 100

    def show_plot(self, frame,frame_value) -> None:
        """Raise the frame for the specific electrode to the front of the canvas."""
        frame.tkraise()
        self.show_frame(f"{frame_value}0")

    def close_frame(self, cont: str) -> None:
        """Destroy the frame."""
        frame = global_show_frames[cont]
        frame.grid_forget()
        # close all matplotlib figures
        plt.close("all")
        # destroy the frames holding the figures
        for frame in global_plot_values:
            frame.destroy()
        # destroy the container holding those frames
        global_plot_container.destroy()

    def apply_x_axis_settings(self) -> None:
        """Apply the x-axis settings."""
        global global_x_left_bound_radiobutton,\
            global_x_right_bound_radiobutton,\
            global_x_left_bound_offset,\
            global_x_right_bound_offset
        match XBound(self.x_axis_control_left.get()):
            case XBound.START_PLUS:
                try:
                    offset = float(self.entry_for_x_axis_control_left_start_plus.get())
                except ValueError:
                    offset = 0
                global_x_left_bound_radiobutton = XBound.START_PLUS
                global_x_left_bound_offset = offset
                self.entry_for_x_axis_control_left_now_minus.delete(0, tk.END)
            case XBound.NOW_MINUS:
                try:
                    offset = float(self.entry_for_x_axis_control_left_now_minus.get())
                except ValueError:
                    offset = 0
                global_x_left_bound_radiobutton = XBound.NOW_MINUS
                global_x_left_bound_offset = offset
                self.entry_for_x_axis_control_left_start_plus.delete(0, tk.END)
        match XBound(self.x_axis_control_right.get()):
            case XBound.START_PLUS:
                try:
                    offset = float(self.entry_for_x_axis_control_right_start_plus.get())
                except ValueError:
                    offset = 0
                global_x_right_bound_radiobutton = XBound.START_PLUS
                global_x_right_bound_offset = offset
                self.entry_for_x_axis_control_right_now_minus.delete(0, tk.END)
            case XBound.NOW_MINUS:
                try:
                    offset = float(self.entry_for_x_axis_control_right_now_minus.get())
                except ValueError:
                    offset = 0
                global_x_right_bound_radiobutton = XBound.NOW_MINUS
                global_x_right_bound_offset = offset
                self.entry_for_x_axis_control_right_start_plus.delete(0, tk.END)

    def apply_y_axis_settings(self) -> None:
        """Apply the y-axis settings."""
        global global_y_norm_radiobutton,\
            global_y_kdm_radiobutton,\
            global_max_norm,\
            global_min_norm,\
            global_max_kdm,\
            global_min_kdm
        match YBound(self.y_axis_control_norm.get()):
            case YBound.AUTOMATIC:
                high_bound = 1.6 #placeholder only
                low_bound = 0.4 #placeholder only
                global_y_norm_radiobutton = YBound.AUTOMATIC
                global_max_norm = high_bound
                global_min_norm = low_bound
                self.entry_for_y_axis_control_norm_manual_max.delete(0, tk.END)
                self.entry_for_y_axis_control_norm_manual_min.delete(0, tk.END)
            case YBound.MANUAL:
                try:
                    high_bound = float(self.entry_for_y_axis_control_norm_manual_max.get())
                except ValueError:
                    high_bound = 2
                try:
                    low_bound = float(self.entry_for_y_axis_control_norm_manual_min.get())
                except ValueError:
                    low_bound = 0
                if high_bound <= low_bound:
                    high_bound = 2
                    low_bound = 0
                global_y_norm_radiobutton = YBound.MANUAL
                global_max_norm = high_bound
                global_min_norm = low_bound
        match YBound(self.y_axis_control_kdm.get()):
            case YBound.AUTOMATIC:
                high_bound = 1.6 #placeholder only
                low_bound = 0.4 #placeholder only
                global_y_kdm_radiobutton = YBound.AUTOMATIC
                global_max_kdm = high_bound
                global_min_kdm = low_bound
                self.entry_for_y_axis_control_kdm_manual_max.delete(0, tk.END)
                self.entry_for_y_axis_control_kdm_manual_min.delete(0, tk.END)
            case YBound.MANUAL:
                try:
                    high_bound = float(self.entry_for_y_axis_control_kdm_manual_max.get())
                except ValueError:
                    high_bound = 2
                try:
                    low_bound = float(self.entry_for_y_axis_control_kdm_manual_min.get())
                except ValueError:
                    low_bound = 0
                if high_bound <= low_bound:
                    high_bound = 2
                    low_bound = 0
                global_y_kdm_radiobutton = YBound.MANUAL
                global_max_kdm = high_bound
                global_min_kdm = low_bound

class FrequencyMapManipulationFrame(ttk.Frame):
    """Frame for Real-time Data Manipulation."""
    def __init__(self, parent: ttk.Frame):
        global global_high_xstart_entry,\
            global_low_xstart_entry,\
            global_high_xend_entry,\
            global_low_xend_entry
        ttk.Frame.__init__(self, parent)         # Initialize the frame
        
        #################################################
        ### Nested Frame for Real-Time adjustment     ###
        ### of voltammogram, polynomial regression,   ###
        ### and savitzky-golay smoothing              ###
        #################################################
        regression_frame = ttk.Frame(self, relief="groove")
        regression_frame.grid(row=7, column=0, columnspan=4, pady=5, padx=5, ipadx=3, sticky="ns")
        regression_frame.rowconfigure(0, weight=1)
        regression_frame.rowconfigure(1, weight=1)
        regression_frame.rowconfigure(2, weight=1)
        regression_frame.columnconfigure(0, weight=1)
        regression_frame.columnconfigure(1, weight=1)
        match global_peak_method:
            case PeakMethod.POLY:
                #--- Title ---#
                ttk.Label(regression_frame, text="Real Time Analysis Manipulation", font=LARGE_FONT).\
                    grid(row=0, column=0, columnspan=4, pady=5, padx=5)
                ###################################################################
                ### Real Time Manipulation of Savitzky-Golay Smoothing Function ###
                ###################################################################
                ttk.Label(regression_frame, text="Savitzky-Golay Window (mV)", font=LARGE_FONT).\
                    grid(row=1, column=0, columnspan=4, pady=1)
                self.smoothing_entry = ttk.Entry(regression_frame, width=10)
                self.smoothing_entry.grid(row=2, column=0, columnspan=4, pady=3)
                self.smoothing_entry.insert(tk.END, str(global_savitzky_golay_window))
                #--- Check for the presence of high and low frequencies ---#
                self.high = global_frequency_list[-1] > 50
                self.low = global_frequency_list[0] <= 50
                ###################################################
                ### If a frequency <= 50Hz exists, grid a frame ###
                ### for low frequency data manipulation         ###
                ###################################################
                if self.low:
                    low_parameter_frame = ttk.Frame(regression_frame)
                    low_parameter_frame.grid(row=3, column=0, columnspan=4, sticky="nsew")
                    low_parameter_frame.rowconfigure(0, weight=1)
                    low_parameter_frame.rowconfigure(1, weight=1)
                    low_parameter_frame.rowconfigure(2, weight=1)
                    low_parameter_frame.columnconfigure(0, weight=1)
                    low_parameter_frame.columnconfigure(1, weight=1)
                    global_show_frames[FRAME_LOW_PARAMETER] = low_parameter_frame
                    #--- points discarded at the beginning of the voltammogram, xstart ---#
                    ttk.Label(low_parameter_frame, text="xstart (V)", font=MEDIUM_FONT).\
                        grid(row=0, column=0)
                    self.low_xstart_entry = ttk.Entry(low_parameter_frame, width=7)
                    self.low_xstart_entry.insert(tk.END, str(global_low_xstart))
                    self.low_xstart_entry.grid(row=1, column=0)
                    global_low_xstart_entry = self.low_xstart_entry
                    #--- points discarded at the beginning of the voltammogram, xend ---#
                    ttk.Label(low_parameter_frame, text="xend (V)", font=MEDIUM_FONT).\
                        grid(row=0, column=1)
                    self.low_xend_entry = ttk.Entry(low_parameter_frame, width=7)
                    self.low_xend_entry.insert(tk.END, str(global_low_xend))
                    self.low_xend_entry.grid(row=1, column=1)
                    global_low_xend_entry = self.low_xend_entry
                ##################################################
                ### If a frequency > 50Hz exists, grid a frame ###
                ### for high frequency data manipulation       ###
                ##################################################
                if self.high:
                    high_parameter_frame = ttk.Frame(regression_frame)
                    high_parameter_frame.grid(row=3, column=0, columnspan=4, sticky="nsew")
                    high_parameter_frame.rowconfigure(0, weight=1)
                    high_parameter_frame.rowconfigure(1, weight=1)
                    high_parameter_frame.rowconfigure(2, weight=1)
                    high_parameter_frame.columnconfigure(0, weight=1)
                    high_parameter_frame.columnconfigure(1, weight=1)
                    global_show_frames[FRAME_HIGH_PARAMETER] = high_parameter_frame
                    #--- points discarded at the beginning of the voltammogram, xstart ---#
                    ttk.Label(high_parameter_frame, text="xstart (V)", font=MEDIUM_FONT).\
                        grid(row=0, column=0)
                    self.high_xstart_entry = ttk.Entry(high_parameter_frame, width=7)
                    self.high_xstart_entry.insert(tk.END, str(global_high_xstart))
                    self.high_xstart_entry.grid(row=1, column=0)
                    global_high_xstart_entry = self.high_xstart_entry
                    #--- points discarded at the beginning of the voltammogram, xend ---#
                    ttk.Label(high_parameter_frame, text="xend (V)", font=MEDIUM_FONT).\
                        grid(row=0, column=1)
                    self.high_xend_entry = ttk.Entry(high_parameter_frame, width=7)
                    self.high_xend_entry.insert(tk.END, str(global_high_xend))
                    self.high_xend_entry.grid(row=1, column=1)
                    global_high_xend_entry = self.high_xend_entry
                ############################################################
                ### If both high and low frequencies are being analyzed, ###
                ### create buttons to switch between the two             ###
                ############################################################
                if self.high and self.low:
                    self.select_low_parameters =\
                        ttk.Button(regression_frame, style="Off.TButton", text="f <= 50Hz",\
                                    command=lambda: self.show_frame(FRAME_LOW_PARAMETER))
                    self.select_low_parameters.grid(row=4, column=0, pady=5, padx=5)
                    self.select_high_parameters =\
                        ttk.Button(regression_frame, style="On.TButton", text="f > 50Hz",\
                                    command=lambda: self.show_frame(FRAME_HIGH_PARAMETER))
                    self.select_high_parameters.grid(row=4, column=1, pady=5, padx=5)
                #--- Button to apply adjustments ---#
                ttk.Button(regression_frame, text="Apply Adjustments",\
                        command=self.adjust_parameters).\
                            grid(row=5, column=0, columnspan=4, pady=10, padx=10)
            case PeakMethod.GAUSS:
                #--- Title ---#
                ttk.Label(regression_frame, text="Expected Peak Location (V)", font=LARGE_FONT).\
                    grid(row=1, column=0, columnspan=4, pady=1)
                self.gauss_peak_entry = ttk.Entry(regression_frame, width=10)
                self.gauss_peak_entry.grid(row=2, column=0, columnspan=4, pady=3)
                self.gauss_peak_entry.insert(tk.END, str(global_gauss_peak))
                #--- Button to apply adjustments ---#
                ttk.Button(regression_frame, text="Apply Adjustments",\
                        command=self.adjust_parameters).\
                            grid(row=5, column=0, columnspan=4, pady=10, padx=10)
        #---Buttons to switch between electrode frames---#
        frame_value = 0
        row_value = 8
        column_value = 0
        for _ in global_plot_values:
            ttk.Button(self, text=global_frame_list[frame_value],\
                       command=lambda frame_value=frame_value:\
                        self.show_plot(global_plot_values[frame_value])).\
                            grid(row=row_value, column=column_value, pady=2, padx=5)
            ## allows .grid() to alternate between
            ## packing into column 1 and column 2
            if column_value == 1:
                column_value = 0
                row_value += 1
            ## if gridding into the 1st column,
            ## grid the next into the 2nd column
            else:
                column_value += 1
            frame_value += 1
        row_value += 1
        #--- Start ---#
        ttk.Button(self, text="Start", style="Fun.TButton", command=self.skeleton_key).\
            grid(row=row_value, column=0, pady=5, padx=5)
        #--- Reset ---#
        ttk.Button(self, text="Reset", style="Fun.TButton", command=self.reset).\
            grid(row=row_value, column=1, pady=5, padx=5)
        row_value += 1
        #--- Quit ---#
        ttk.Button(self, text="Quit Program", command=lambda: os._exit(0)).\
            grid(row=row_value, column=0, columnspan=4, pady=5)
        for row in range(row_value):
            row += 1
            self.rowconfigure(row, weight=1)
        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)
                                        ###################################################
                                        ###################################################
                                        ###### Real Time Data Manipulation Functions ######
                                        ###################################################
                                        ###################################################
    #################################################
    ### Adjustment of points discarded at the     ###
    ### beginning and end of Regression Analysis  ###
    #################################################
    def adjust_parameters(self) -> None:
        """Adjust the parameters used to visualize the raw voltammogram,
        smoothed currents, and polynomial fit.
        """
        global global_low_xstart,\
            global_high_xstart,\
            global_low_xend,\
            global_high_xend,\
            global_savitzky_golay_window,\
            global_gauss_peak
        match global_peak_method:
            case PeakMethod.POLY:
                ###############################################
                ### Polynomical Regression Range Parameters ###
                ###############################################
                if self.low:
                    #--- parameters for frequencies equal or below 50Hz ---#
                    # xstart/xend adjust the points at the start and end of
                    # the voltammogram/smoothed currents, respectively
                    global_low_xend = float(self.low_xend_entry.get())
                    global_low_xstart = float(self.low_xstart_entry.get())
                if self.high:
                    #--- parameters for frequencies above 50Hz ---#
                    global_high_xstart = float(self.high_xstart_entry.get())
                    global_high_xend = float(self.high_xend_entry.get())
                #######################################
                ### Savitzky-Golay Smoothing Window ###
                #######################################
                global_savitzky_golay_window = float(self.smoothing_entry.get())
                print(f"adjust_parameters: SG_Window (mV) {global_savitzky_golay_window}")
            case PeakMethod.GAUSS:
                global_gauss_peak = float(self.gauss_peak_entry.get())

    def reset(self) -> None:
        """Reset the experiment."""
        global global_key,\
            global_poison_pill,\
            global_analysis_already_initiated,\
            global_x_left_bound_radiobutton,\
            global_x_left_bound_offset,\
            global_x_right_bound_radiobutton,\
            global_x_right_bound_offset
        global_key = 0
        global_poison_pill = True
        global_analysis_already_initiated = False # reset the start variable
        global_x_left_bound_radiobutton = XBound.START_PLUS
        global_x_left_bound_offset = 0
        global_x_right_bound_radiobutton = XBound.NOW_MINUS
        global_x_right_bound_offset = 0
        # Raise the initial user input frame
        self.show_frame(FRAME_INPUT)
        self.close_frame(str(global_analysis_method))

    def show_frame(self, cont: str) -> None:
        """Raise the frame to the front of the canvas."""
        frame = global_show_frames[cont]
        frame.tkraise()
        if cont == FRAME_LOW_PARAMETER:
            self.select_low_parameters.config(style="On.TButton")
            self.select_high_parameters.config(style="Off.TButton")
        elif cont == FRAME_HIGH_PARAMETER:
            self.select_low_parameters.config(style="Off.TButton")
            self.select_high_parameters.config(style="On.TButton")

    def skeleton_key(self) -> None:
        """Start the data analysis and visualization."""
        global global_key, global_poison_pill, global_analysis_already_initiated
        if global_analysis_already_initiated:
            print("Program has already been initiated")
        else:
            ######################################################################
            ### Initialize Animation (Visualization) for each electrode figure ###
            ######################################################################
            fig_count = 0                   # index value for the frame
            for figure in global_figures:
                fig, _ = figure
                electrode = global_electrode_list[fig_count]
                global_animations.append(ElectrochemicalAnimation(fig, electrode,\
                                                                  resize_interval=0, fargs=None))
                fig_count += 1
            global_analysis_already_initiated = True
            #--- reset poison pill variables --#
            global_poison_pill = False
            # tells Generate() to start data analysis:
            if global_key == 0:
                global_key += 100

    def show_plot(self, frame: ttk.Frame) -> None:
        """Raise the frame for the specific electrode to the front of the canvas."""
        frame.tkraise()

    def close_frame(self, cont: str) -> None:
        """Destroy the frame."""
        frame = global_show_frames[cont]
        frame.grid_forget()
        for value in global_plot_values:
            value.destroy()
        global_plot_container.destroy()
#-------------------------------------------------------------------------------------------------#
class ContinuousScanVisualizationFrame(ttk.Frame):
    """Electrode Frame Class for data visualization, displayed next to
    the RealTimeManipulationFrame. Embeds a canvas within the tkinter
    MainWindow containing figures that visualize the data for that electrode.
    """
    def __init__(self, electrode_frame: str, count: int, parent: ttk.Frame):
        ttk.Frame.__init__(self, parent)
        #--- for resize ---#
        self.columnconfigure(0, weight=2)
        self.columnconfigure(1, weight=1)
        self.rowconfigure(2, weight=2)
        ttk.Label(self, text=electrode_frame, font=HUGE_FONT).\
            grid(row=0, column=0, pady=5, sticky="n")
        ttk.Label(self, text="", font=MEDIUM_FONT).\
            grid(row=0, column=1, pady=3, sticky="ne")
        #--- Voltammogram, Raw Peak Height, and Normalized Figure and Artists ---#
        # Call the figure and artists for the electrode
        # and place the artists within the frame
        fig, _ = global_figures[count]
        canvas = FigureCanvasTkAgg(fig, self)
        # initial draw call to create the artists that will be blitted:
        canvas.draw()
        # does not affect size of figure within plot container:
        canvas.get_tk_widget().grid(row=1, columnspan=2, pady=6, ipady=5, sticky="news")
        if len(global_frequency_list) > 1:
            #--- Ratiometric Figure and Artists ---#
            # Call the figure and artists for the electrode
            # and place the artists within the frame
            fig, _ = global_ratiometric_figures[count]
            canvas = FigureCanvasTkAgg(fig, self)
            canvas.draw()
            # does not affect size of figure within plot container:
            canvas.get_tk_widget().grid(row=2, columnspan=2, pady=6, ipady=5, sticky="sew")
            #--- add weight to the second row for resizing ---#
            self.rowconfigure(2, weight=2)

class FrequencyMapVisualizationFrame(ttk.Frame):
    """Displayed next to the RealTimeManipulationFrame.
    """
    def __init__(self, electrode_frame: str, count: int, parent: ttk.Frame):
        ttk.Frame.__init__(self, parent)
        #--- for resize ---#
        self.columnconfigure(0, weight=2)
        self.columnconfigure(1, weight=1)
        self.rowconfigure(2, weight=2)
        ttk.Label(self, text=electrode_frame, font=HUGE_FONT).\
            grid(row=0, column=0, pady=5, sticky="n")
        ttk.Label(self, text="", font=MEDIUM_FONT).\
            grid(row=0, column=1, pady=3, sticky="ne")
        #--- Voltammogram, Raw Peak Height, and Normalized Figure and Artists ---#
        #Call the figure and artists for the electrode:
        fig, _ = global_figures[count]
        #and place the artists within the frame:
        canvas = FigureCanvasTkAgg(fig, self)
        #initial draw call to create the artists that will be blitted:
        canvas.draw()
        #does not affect size of figure within plot container:
        canvas.get_tk_widget().grid(row=1, columnspan=2, pady=6, ipady=5, sticky="news")
                                    #############################################################
                                    #############################################################
                                    ###                   End of GUI Classes                  ###
                                    #############################################################
                                    #############################################################
#-------------------------------------------------------------------------------------------------#
                                    #############################################################
                                    #############################################################
                                    ### Creation of Matplotlib Canvas, Figures, Axes, Artists ###
                                    ### and all other decorators (e.g., axis labels, titles)  ###
                                    #############################################################
                                    #############################################################
class InitializeContinuousCanvas():
    """Initialize the canvas, figure, axes, and artists for each electrode."""
    def __init__(self):
        global global_text_file_export,\
            global_offset_normalized_data_list,\
            global_animations,\
            global_frame_reference,\
            global_plot_container,\
            global_kdm_list,\
            global_ratiometric_plots,\
            global_ratiometric_figures,\
            global_normalized_ratiometric_data_list,\
            global_normalized_data_list,\
            global_data_list,\
            global_plot_list_continuous_scan,\
            global_file_list,\
            global_figures,\
            global_frame_list,\
            global_sample_list,\
            global_plot_frames,\
            global_plot_values,\
            global_gauss_solver,\
            global_peak_list,\
            global_xstart,\
            global_xend,\
            global_xstart_entry,\
            global_xend_entry,\
            global_gauss_peak,\
            global_gauss_baseline,\
            global_gauss_maxheight        
        ##############################################
        ### Generate global lists for data storage ###
        ##############################################
        self.length = len(global_frequency_list)
        #--- Animation list ---#
        global_animations = []
        #--- Figure lists ---#
        global_figures = []
        global_ratiometric_figures = []
        ############################################
        ### Create global lists for data storage ###
        ############################################
        # xstart xend
        global_xstart = [[0.0]]*global_electrode_count
        global_xend = [[0.0]]*global_electrode_count
        global_xstart_entry = [[0.0]]*global_electrode_count
        global_xend_entry = [[0.0]]*global_electrode_count
        # gauss method solver
        global_gauss_solver = []
        # gauss peak, baseline, maxheight
        global_gauss_peak = [[-0.3]]*global_electrode_count
        global_gauss_baseline = [[0.3]]*global_electrode_count
        global_gauss_maxheight = [[5.0]]*global_electrode_count
        # gauss peak location
        global_peak_list = [[[0.0]]]*global_electrode_count
        # Peak Height/AUC data (after smoothing and polynomial regression or gauss method):
        #global_data_list = [0]*global_electrode_count
        global_data_list = [[[0.0]]]*global_electrode_count
        # normalized data:
        global_normalized_data_list = [[[0.0]]]*global_electrode_count
        # to hold data with low frequency offset:
        global_offset_normalized_data_list = [[0]]*global_electrode_count
        # uses ratio of normalized peak heights:
        global_normalized_ratiometric_data_list = []
        # to hold the data for kinetic differential measurements:
        global_kdm_list = []
        for num in range(global_electrode_count):
            # xstart xend
            global_xstart[num] = [0.0]*self.length
            global_xend[num] = [0.0]*self.length
            global_xstart_entry[num] = [0.0]*self.length
            global_xend_entry[num] = [0.0]*self.length
            # gauss peak, baseline, maxheight
            global_gauss_peak[num] = [-0.3]*self.length
            global_gauss_baseline[num] = [0.3]*self.length
            global_gauss_maxheight[num] = [5.0]*self.length
            # gauss peak location
            global_peak_list[num] = [[0.0]]*self.length
            # a data list for each electrode:
            global_data_list[num] = [[0.0]]*self.length
            global_normalized_data_list[num] = [[0.0]]*self.length
            global_offset_normalized_data_list[num] = [0.0]*global_number_of_files_to_process
            # a data list for each frequency for that electrode:
            for count in range(self.length):
                # use [0]*number_of_files_to_process to preallocate list space:
                global_peak_list[num][count] = [0.0]*global_number_of_files_to_process
                global_data_list[num][count] = [0.0]*global_number_of_files_to_process
                global_normalized_data_list[num][count] = [0.0]*global_number_of_files_to_process
        for num in range(global_electrode_count):
            global_normalized_ratiometric_data_list.append([])
            global_kdm_list.append([])
        #--- Lists of Frames and Artists ---#
        global_plot_list_continuous_scan = []
        global_ratiometric_plots = []
        global_frame_list = []
        #--- Misc Lists ---#
        global_file_list = []
        global_sample_list = []        # For plotting Peak Height vs. sample rate
        ######################################################
        ### Create a figure and artists for each electrode ###
        ######################################################
        for num in range(global_electrode_count):
            electrode = global_electrode_list[num]
            figure = self.make_figure(electrode)
            global_figures.append(figure)
            if len(global_frequency_list) > 1:
                ratio_figure = self.make_ratiometric_figure(electrode)
                global_ratiometric_figures.append(ratio_figure)
        #####################################################
        ### Create a frame for each electrode and embed   ###
        ### within it the figure containing its artists   ###
        #####################################################
        global_plot_frames = {}                # Dictionary of frames for each electrode
        global_plot_values = []                # create a list of frames
        #--- Create a container that can be created and destroyed
        # when Start() or Reset() is called, respectively ---#
        global_plot_container = ttk.Frame(global_container, relief="groove")
        global_plot_container.grid(row=0, column=1, sticky="nsew")
        global_plot_container.rowconfigure(0, weight=1)
        global_plot_container.columnconfigure(0, weight=1)
        frame_count = 0
        # Iterate through the frame of each electrode:
        for electrode_frame in global_frame_list:
            #--- create an instance of the frame and append it to the global frame dictionary ---#
            ## PlotContainer is the 'parent' frame
            global_frame_reference = ContinuousScanVisualizationFrame(electrode_frame,\
                                                                        frame_count,\
                                                                        global_plot_container)
            # sticky must be "nsew" so it expands and contracts with resize:
            global_frame_reference.grid(row=0, column=0, sticky="nsew")
            global_plot_frames[electrode_frame] = global_frame_reference
            frame_count += 1
        #--- Create a list containing the Frame objects for each electrode ---#
        for _, frame in global_plot_frames.items():
            global_plot_values.append(frame)
        #################################
        ### Initiate .txt File Export ###
        #################################
        #--- If the user has indicated that text file export should be activated ---#
        if global_text_file_export_activated:
            print("Initializing Text File Export")
            global_text_file_export = TextFileExport().initialize()

    def make_figure(self,\
                electrode: int)\
                -> Tuple[matplotlib.figure.Figure, List[List[matplotlib.axes.Axes]]]:
        """Create the figure and artist objects for the given electrode."""
        global global_empty_plots
        try:
            ########################
            ### Setup the Figure ###
            ########################
            length = self.length
            ## figsize=(width, height)
            fig, npaxs = plt.subplots(nrows=3, ncols=length, squeeze=False, figsize=(9, 4.5))
            axs: List[List[matplotlib.axes.Axes]] = npaxs.tolist()
            ### adjust the spacing between subplots:
            plt.subplots_adjust(bottom=0.1, hspace=0.6, wspace=0.3)
            #######################
            ### Set axis labels ###
            #######################
            axs[0][0].set_ylabel("Current\n(µA)", fontweight="bold")
            match global_plot_summary_mode:
                case PlotSummaryMode.PHE:
                    axs[1][0].set_ylabel("Peak Height\n(µA)", fontweight="bold")
                case PlotSummaryMode.AUC:
                    axs[1][0].set_ylabel("AUC (a.u.)", fontweight="bold")
            axs[2][0].set_ylabel("Normalized", fontweight="bold")
            ##########################################
            ### Set subplot axes for each frequency ###
            ##########################################
            electrode_plot: List[Tuple[matplotlib.lines.Line2D,\
                                                    matplotlib.lines.Line2D,\
                                                    matplotlib.lines.Line2D,\
                                                    matplotlib.lines.Line2D,\
                                                    matplotlib.lines.Line2D,\
                                                    matplotlib.lines.Line2D,\
                                                    Polygon,\
                                                    matplotlib.lines.Line2D]] = []
            for freq in range(length):
                frequency = global_frequency_list[freq]
                axs[0][freq].set_xlabel("Potential (V)")
                #--- if the resize interval is larger than the number of files, ---#
                #--- make the x lim the number of files (& vice versa)          ---#
                if global_resize_interval > global_number_of_files_to_process:
                    xlim_factor = global_number_of_files_to_process
                else:
                    xlim_factor = global_resize_interval
                match global_x_axis_mode:
                    case PlotTimeReportingMode.EXPERIMENT_TIME:
                        axs[1][freq].set_xlim(0, xlim_factor*global_sample_rate/3600\
                                                + global_sample_rate/7200)
                        axs[2][freq].set_xlim(0, xlim_factor*global_sample_rate/3600\
                                                + global_sample_rate/7200)
                    case PlotTimeReportingMode.FILE_NUMBER:
                        axs[1][freq].set_xlim(-0.05, xlim_factor+0.1)
                        axs[2][freq].set_xlim(-0.05, xlim_factor+0.1)
                        axs[2][freq].set_xlabel("File Number")
                #################################################################################
                ###       Analyze the first file and create the Y limits of the subplots      ###
                ###               depending on the data range of the first file               ###
                #################################################################################
                self.initialize_subplots(axs, frequency, electrode, freq)
                #---Set Subplot Title---#
                axs[0][freq].set_title(str(frequency) + " Hz", fontweight="bold")
                #---Initiate the subplots---#
                # this assigns a Line2D artist object to the artist object (Axes)
                smooth, = axs[0][freq].plot([], [], "ko", markersize=2)
                regress, = axs[0][freq].plot([], [], "b-")
                center, = axs[0][freq].plot([], [], "r-")
                peak, = axs[1][freq].plot([], [], "ko", markersize=1)
                peak_injection, = axs[1][freq].plot([], [], "bo", markersize=1)
                normalization, = axs[2][freq].plot([], [], "ko", markersize=1)
                norm_injection, = axs[2][freq].plot([], [], "ro", markersize=1)
                #--- shading for AUC ---#
                verts = [(0, 0), *zip([], []), (0, 0)]
                poly = Polygon(verts, alpha=0.5)
                axs[0][freq].add_patch(poly)
                #####################################################
                ### Create a list of the primitive artists        ###
                ### (Line2D objects) that will be returned        ###
                ### to ElectrochemicalAnimation to be visualized  ###
                #####################################################
                # this is the list that will be returned as _drawn_artists
                #  to the Funcanimation class
                plots: Tuple[matplotlib.lines.Line2D,\
                                                    matplotlib.lines.Line2D,\
                                                    matplotlib.lines.Line2D,\
                                                    matplotlib.lines.Line2D,\
                                                    matplotlib.lines.Line2D,\
                                                    matplotlib.lines.Line2D,\
                                                    Polygon,\
                                                    matplotlib.lines.Line2D] =\
                                                          (smooth,\
                                                           regress,\
                                                            peak,\
                                                            peak_injection,\
                                                            normalization,\
                                                            norm_injection,\
                                                            poly,\
                                                            center)
                #--- And append that list to keep a global reference ---#
                # 'plots' is a list of artists that are passed to animate:
                electrode_plot.append(plots)
                electrode_frame = f"Electrode {electrode}"
                if electrode_frame not in global_frame_list:
                    global_frame_list.append(electrode_frame)
                #--- Create empty plots to return to animate for initializing---#
                global_empty_plots = [smooth, regress, peak, normalization,center]
            # global_plot_list is a list of lists containing 'plots' for each electrode:
            global_plot_list_continuous_scan.append(electrode_plot)
            #-- Return both the figure and the axes to be stored as global variables --#
            return fig, axs
        except Exception as exception:
            internal_error("exception in make_figure:" + str(exception))

    def make_ratiometric_figure(self,\
                    electrode: int) -> Tuple[matplotlib.figure.Figure, List[matplotlib.axes.Axes]]:
        """Create the figure and artist objects for
        the ratiometric data for the given electrode
        (e.g., Kinetic Differential Measurement, Normalized Ratio).
        """
        try:
            figure, npaxs = plt.subplots(nrows=1, ncols=2, squeeze=False, figsize=(8.5, 1.85))
            axs: List[matplotlib.axes.Axes] = npaxs.flatten().tolist()
            ### adjust the spacing between subplots:
            plt.subplots_adjust(bottom=0.3, hspace=0.6, wspace=0.3)
            ################################################
            ## Set the X and Y axes for the Ratiometric  ##
            ## Plots (KDM and Norm Ratio)                 ##
            ################################################
            axs[0].set_ylabel("% Signal", fontweight="bold")
            axs[1].set_ylabel("% Signal", fontweight="bold")
            x_left_bound, x_right_bound = compute_x_bounds(1)
            axs[0].set_xlim(x_left_bound, x_right_bound)
            axs[1].set_xlim(x_left_bound, x_right_bound)
            match global_x_axis_mode:
                case PlotTimeReportingMode.EXPERIMENT_TIME:
                    axs[0].set_xlabel("Time (h)")
                    axs[1].set_xlabel("Time (h)")
                case PlotTimeReportingMode.FILE_NUMBER:
                    axs[0].set_xlabel("File Number")
                    axs[1].set_xlabel("File Number")
            axs[0].set_ylim(100*global_min_norm, 100*global_max_norm)
            axs[1].set_ylim(100*global_min_kdm, 100*global_max_kdm)
            axs[0].set_title("Normalized Ratio")
            axs[1].set_title("KDM")
            #####################################################
            ### Create the primitive artists (Line2D objects) ###
            ### that will contain the data that will be       ###
            #### visualized by ElectrochemicalAnimation       ###
            #####################################################
            # normalized ratio of high and low freq's:
            norm_ratiometric_plot, = axs[0].plot([], [], "ro", markersize=1)
            kdm, = axs[1].plot([], [], "ro", markersize=1)
            # if injection_point =! None, these will
            # visualize the points after the injection
            norm_injection, = axs[0].plot([], [], "bo", markersize=1)
            kdm_injection, = axs[1].plot([], [], "bo", markersize=1)
            ratio_plots = [norm_ratiometric_plot, norm_injection, kdm, kdm_injection]
            global_ratiometric_plots.append(ratio_plots)
            return figure, axs
        except Exception as exception:
            internal_error("error in make_ratiometric_figure:" + str(exception))

    def initialize_subplots(self,\
                            axs: List[List[matplotlib.axes.Axes]],\
                            frequency: int,\
                            electrode: int,\
                            subplot_count: int) -> None:
        """Initialize the y limits of each figure depending on the y values of the first file."""
        try:
            filename = make_file_name(1, electrode, frequency)
            myfile = global_file_path + filename
            if file_is_complete(myfile):
                self.run_initialization(myfile, axs, subplot_count, electrode, frequency)
        except FileNotFoundError:
            print(f"could not find file for electrode {electrode}")
            #--- If search time has not met the search limit keep searching ---#
            root.after(1000, self.initialize_subplots, axs, frequency, electrode, subplot_count)

    def run_initialization(self,\
                            myfile: str,\
                            axs: List[List[matplotlib.axes.Axes]],\
                            subplot_count: int,\
                            electrode: int,\
                            frequency: int) -> None:
        """Initialize subplots."""
        global global_high_xstart,\
            global_high_xend,\
            global_low_xstart,\
            global_low_xend,\
            global_xstart,\
            global_xend,\
            global_gauss_solver
        
        try:
            #########################
            ### Retrieve the data ###
            #########################
            potentials, currents, _ = read_data(myfile, electrode)
            ##########################################
            ### Set the x axes of the voltammogram ###
            ##########################################
            min_potential = min(potentials)
            max_potential = max(potentials)
            #-- Reverse voltammogram to match the 'Texas' convention --#
            axs[0][subplot_count].set_xlim(max_potential, min_potential)
            #######################################
            ### Get the high and low potentials ###
            #######################################
            # if not global_already_reset:
            global_xstart[global_electrode_dict[electrode]][global_frequency_dict[frequency]] = max(potentials)
            global_xend[global_electrode_dict[electrode]][global_frequency_dict[frequency]] = min(potentials)
            xstart = global_xstart[global_electrode_dict[electrode]][global_frequency_dict[frequency]]
            xend = global_xend[global_electrode_dict[electrode]][global_frequency_dict[frequency]]

            cut_value = 0
            for value in potentials:
                if abs(value) < EPSILON:
                    cut_value += 1
            if cut_value > 0:
                potentials = potentials[:-cut_value]
                currents = currents[:-cut_value]
            adjusted_potentials = [value for value in potentials if xend <= value <= xstart]
            #################################
            ### Create Gauss solvers once ###
            #################################
            length = len(adjusted_potentials)
            global_gauss_solver.append(MultiGausFitNoise_LSE(length, global_options_Gauss))
            #########################################
            ### Savitzky-Golay smoothing          ###
            #########################################
            smooth_currents = savgol_filter(currents, 15, global_savitzky_golay_degree)
            potential_to_current_map_smoothed = dict(zip(potentials, smooth_currents))
            #######################################
            ### adjust the smooth currents to   ###
            ### match the adjusted potentials   ###
            #######################################
            adjusted_currents: List[float]
            adjusted_currents = []
            for potential in adjusted_potentials:
                adjusted_currents.append(potential_to_current_map_smoothed[potential])
            ######################
            ### Polynomial fit ###
            ######################
            polynomial_coeffs = np.polyfit(adjusted_potentials, adjusted_currents, POLYFIT_DEGREE)
            eval_regress = np.polyval(polynomial_coeffs, adjusted_potentials).tolist()
            fit_half = round(len(eval_regress)/2)
            min1 = min(eval_regress[:-fit_half])
            min2 = min(eval_regress[fit_half:])
            max1 = max(eval_regress[:-fit_half])
            max2 = max(eval_regress[fit_half:])
            match global_plot_summary_mode:
                case PlotSummaryMode.PHE:
                    peak_height = max(max1, max2) - min(min1, min2)
                    data = peak_height
                case PlotSummaryMode.AUC:
                    auc_index = 1
                    auc: float
                    auc = 0
                    auc_potentials = adjusted_potentials
                    auc_min = min(adjusted_currents)
                    auc_currents = [y - auc_min for y in adjusted_currents]
                    while auc_index <= len(auc_currents) - 1:
                        auc_height = (auc_currents[auc_index] + auc_currents[auc_index - 1])/2
                        auc_width = auc_potentials[auc_index] - auc_potentials[auc_index - 1]
                        auc += (auc_height * auc_width)
                        auc_index += 1
                    data = auc
            #--- calculate the baseline current ---#
            minimum_current = min(min1, min2)
            maximum_current = max(max1, max2)
            #- Voltammogram -#
            axs[0][subplot_count].set_ylim(minimum_current-abs(global_min_raw*minimum_current),\
                                            maximum_current+abs(global_max_raw*maximum_current))
            #- PHE/AUC Data -#
            axs[1][subplot_count].set_ylim(data-abs(global_min_data*data),\
                                            data+abs(global_max_data*data))
            #- Normalized Data -#
            axs[2][subplot_count].set_ylim(global_min_norm, global_max_norm)
        except Exception as exception:
            internal_error("run_initialization" + str(exception))

class InitializeFrequencyMapCanvas():
    """Real Time Data Animation Canvas for Frequency Map Analysis."""
    def __init__(self):
        global global_text_file_export,\
            global_file_list,\
            global_animations,\
            global_frame_reference,\
            global_plot_container,\
            global_data_list,\
            global_plot_list_frequency_map,\
            global_figures,\
            global_frame_list,\
            global_plot_frames,\
            global_plot_values,\
            global_gauss_solver,\
            global_peak_list
        ##############################################
        ### Generate global lists for data storage ###
        ##############################################
        self.length = len(global_frequency_list)
        #--- Animation list ---#
        global_animations = []
        #--- file list ---#
        global_file_list = [0]*global_number_of_files_to_process
        #--- Figure lists ---#
        global_figures = []
        ############################################
        ### Create global lists for data storage ###
        ############################################
         # gauss method solver:
        global_gauss_solver = []
        # gauss peak location
        global_peak_list = [[[0.0]]]*global_electrode_count
        # Peak Height/AUC data (after smoothing and polynomial regression):
        global_data_list = [[[0.0]]]*global_electrode_count
        for num in range(global_electrode_count):
            # a data list for each electrode:
            global_peak_list[num] = [[0.0]]*self.length
            global_data_list[num] = [[0.0]]*self.length
            # a data list for each frequency for that electrode:
            for count in range(self.length):
                global_peak_list[num][count] = [0.0]*global_number_of_files_to_process
                global_data_list[num][count] = [0.0]*global_number_of_files_to_process
        #--- Lists of Frames and Artists ---#
        global_plot_list_frequency_map = []
        global_frame_list = []
        ######################################################
        ### Create a figure and artists for each electrode ###
        ######################################################
        for num in range(global_electrode_count):
            electrode = global_electrode_list[num]
            figure = self.make_figure(electrode)
            global_figures.append(figure)
        #####################################################
        ### Create a frame for each electrode and embed   ###
        ### within it the figure containing its artists   ###
        #####################################################
        global_plot_frames = {}                # Dictionary of frames for each electrode
        global_plot_values = []                # create a list of frames
        #--- Create a container that can be created and destroyed
        # when Start() or Reset() is called, respectively ---#
        global_plot_container = ttk.Frame(global_container, relief="groove")
        global_plot_container.grid(row=0, column=1, sticky="nsew")
        global_plot_container.rowconfigure(0, weight=1)
        global_plot_container.columnconfigure(0, weight=1)
        frame_count = 0
        # Iterate through the frame of each electrode
        for electrode_frame in global_frame_list:
            #--- create an instance of the frame and append it to the global frame dictionary ---#
            # PlotContainer is the 'parent' frame
            global_frame_reference = FrequencyMapVisualizationFrame(electrode_frame,\
                                                                    frame_count,\
                                                                     global_plot_container)
            # sticky must be 'nsew' so it expands and contracts with resize
            global_frame_reference.grid(row=0, column=0, sticky="nsew")
            global_plot_frames[electrode_frame] = global_frame_reference
            frame_count += 1
        #--- Create a list containing the Frame objects for each electrode ---#
        for _, frame in global_plot_frames.items():
            global_plot_values.append(frame)
        #################################
        ### Initiate .txt File Export ###
        #################################
        #--- If the user has indicated that text file export should be activated ---#
        if global_text_file_export_activated:
            global_text_file_export = TextFileExport().initialize()
        #declarations for fields that will be initialized in other methods
        self.list_val: int

    def make_figure(self,\
                electrode: int)\
                -> Tuple[matplotlib.figure.Figure, List[List[matplotlib.axes.Axes]]]:
        """Create the figure and artist objects for the given electrode."""
        global global_empty_plots
        try:
            ##########################################
            ### Setup the Figure for voltammograms ###
            ##########################################
            fig, npaxs = plt.subplots(nrows=2, ncols=1, squeeze=False, figsize=(9, 4.5))
            axs: List[List[matplotlib.axes.Axes]] = npaxs.tolist()
            ### adjust the spacing between subplots
            plt.subplots_adjust(bottom=0.2, hspace=0.6, wspace=0.3)
            #######################
            ### Set axis labels ###
            #######################
            axs[0][0].set_ylabel("Current (µA)", fontweight="bold")
            axs[0][0].set_xlabel("Voltage (V)", fontweight="bold")
            axs[1][0].set_ylabel("Charge (µC)", fontweight="bold")
            axs[1][0].set_xlabel("Frequency (Hz)", fontweight="bold")
            ##########################################
            ### Set suplot axes for each frequency ###
            ##########################################
            electrode_plot = []
            axs[1][0].set_xscale("log")
            #################################################################################
            #################################################################################
            ###       Analyze the first file and create the Y limits of the subplots      ###
            ###               depending on the data range of the first file               ###
            #################################################################################
            self.initialize_subplots(axs, electrode)
            #################################################################################
            #################################################################################
            #---Initiate the subplots---#
            # this assigns a Line2D artist object to the artist object (Axes)
            smooth, = axs[0][0].plot([], [], "ko", markersize=2)
            regress, = axs[0][0].plot([], [], "r-")
            charge, = axs[1][0].plot([], [], "ko", markersize=1)
            #####################################################
            ### Create a list of the primitive artists        ###
            ### (Line2D objects) that will be returned        ###
            ### to ElectrochemicalAnimation to be visualized  ###
            #####################################################
            # this is the list that will be returned as _drawn_artists to the Funcanimation class
            plots: List[matplotlib.artist.Artist] = [smooth, regress, charge]
            #--- And append that list to keep a global reference ---#
            # 'plots' is a list of artists that are passed to animate
            electrode_plot.append(plots)
            electrode_frame = f"Electrode {electrode}"
            if electrode_frame not in global_frame_list:
                global_frame_list.append(electrode_frame)
            #--- Create empty plots to return to animate for initializing---#
            global_empty_plots = [smooth, regress, charge]
            global_plot_list_frequency_map.append((smooth, regress, charge))
            #-- Return both the figure and the axes to be stored as global variables --#
            return fig, axs
        except Exception as exception:
            internal_error("make_figure" + str(exception))

    def initialize_subplots(self, axs: List[List[matplotlib.axes.Axes]], electrode: int) -> None:
        """Initialize the y limits of each figure depending on the y values of the first file."""
        self.list_val = column_index_for_current(electrode)
        try:
            frequency = global_frequency_list[0]
            filename = make_file_name(1, electrode, frequency)
            myfile = global_file_path + filename
            if file_is_complete(myfile):
                self.run_initialization(myfile, axs, electrode,frequency)
        except FileNotFoundError:
            print(f"could not find file for electrode {electrode}")
            #--- If search time has not met the search limit keep searching ---#
            root.after(1000, self.initialize_subplots, axs, electrode)

    def run_initialization(self,\
                        myfile: str,\
                        axes: List[List[matplotlib.axes.Axes]],\
                        electrode: int,\
                        frequency: int) -> None:
        """Initialize subplots."""
        global global_high_xstart, global_high_xend, global_low_xstart, global_low_xend, global_gauss_solver
        try:
            #########################
            ### Retrieve the data ###
            #########################
            potentials, currents, _ = read_data(myfile, electrode)
            ##########################################
            ### Set the x axes of the voltammogram ###
            ##########################################
            min_potential = min(potentials)
            max_potential = max(potentials)
            #-- Reverse voltammogram to match the 'Texas' convention --#
            axes[0][0].set_xlim(max_potential, min_potential)
            #######################################
            ### Get the high and low potentials ###
            #######################################
            #-- set the local variables to the global ---#
            xstart = max(potentials)
            xend = min(potentials)
            global_low_xstart = xstart
            global_high_xstart = xstart
            global_low_xend = xend
            global_high_xend = xend
            cut_value = 0
            for value in potentials:
                if abs(value) < EPSILON:
                    cut_value += 1
            if cut_value > 0:
                potentials = potentials[:-cut_value]
                currents = currents[:-cut_value]
            adjusted_potentials = [value for value in potentials if xend <= value <= xstart]
            #################################
            ### Create Gauss solvers once ###
            #################################
            length = len(adjusted_potentials)
            global_gauss_solver.append(MultiGausFitNoise_LSE(length, global_options_Gauss))
            #########################################
            ### Savitzky-Golay smoothing          ###
            #########################################
            smooth_currents = savgol_filter(currents, 15, global_savitzky_golay_degree)
            potential_to_current_map_smoothed = dict(zip(potentials, smooth_currents))
            #######################################
            ### adjust the smooth currents to   ###
            ### match the adjusted potentials   ###
            #######################################
            adjusted_currents = []
            for potential in adjusted_potentials:
                adjusted_currents.append(potential_to_current_map_smoothed[potential])
            ######################
            ### Polynomial fit ###
            ######################
            polynomial_coeffs = np.polyfit(adjusted_potentials, adjusted_currents, POLYFIT_DEGREE)
            eval_regress = np.polyval(polynomial_coeffs, adjusted_potentials).tolist()
            fit_half = round(len(eval_regress)/2)
            min1 = min(eval_regress[:-fit_half])
            min2 = min(eval_regress[fit_half:])
            max1 = max(eval_regress[:-fit_half])
            max2 = max(eval_regress[fit_half:])
            match global_plot_summary_mode:
                case PlotSummaryMode.PHE:
                    pass
                case PlotSummaryMode.AUC:
                    auc_index = 1
                    auc = 0
                    auc_potentials = [abs(potential) for potential in adjusted_potentials]
                    auc_min = min(adjusted_currents)
                    auc_currents = [y - auc_min for y in adjusted_currents]
                    while auc_index <= len(auc_currents) - 1:
                        auc_height = (auc_currents[auc_index] + auc_currents[auc_index - 1])/2
                        auc_width = auc_potentials[auc_index] - auc_potentials[auc_index - 1]
                        auc += (auc_height * auc_width)
                        auc_index += 1
            #--- calculate the baseline current ---#
            minimum_current = min(min1, min2)
            maximum_current = max(max1, max2)
            peak_current = maximum_current - minimum_current
            charge = peak_current/(global_frequency_list[0])
            ## Reverse voltammogram to match the 'Texas' convention ##
            axes[0][0].set_xlim(max_potential, min_potential)
            # voltammogram:
            axes[0][0].set_ylim(minimum_current-abs(global_min_raw*minimum_current),\
                                maximum_current+abs(global_max_raw*maximum_current))
            ## set the limits of the Lovric plot ##
            axes[1][0].set_ylim(charge-abs(global_min_data*charge),\
                                charge+abs(global_max_data*charge))
            axes[1][0].set_xlim(int(global_frequency_list[0]), int(global_frequency_list[-1]))
        except Exception as exception:
            print("Error in run_initialization", str(exception))
                                #############################################################
                                #############################################################
                                ###              END OF INITIALIZATION FUNCTIONS          ###
                                #############################################################
                                #############################################################
#-------------------------------------------------------------------------------------------------#
class ElectrochemicalAnimation():
    """Animation class for real time data visualization."""
    def __init__(self, fig: matplotlib.figure.Figure, electrode: int,\
                  generator=None, func=None, resize_interval=0, fargs=None):
        self.electrode = electrode                       # Electrode for this class instance
        self.num = global_electrode_dict[self.electrode] # Electrode index value
        self.file = STARTING_FILE_NUMBER                 # Starting File
        self.index = 0                                   # File Index Value
        self.count = 0                                   # Frequency index value
        self.frequency_limit = len(global_frequency_list) - 1 # -1 so it matches the index value
        ### Lists for sample rate (time passed)  ###
        ### and file count for each electrode    ###
        self.sample_list: List[float] = []
        self.file_list: List[int] = []
        self.frequency_axis: List[int] = []
        self.charge_axis: List[float] = []
        ##############################
        ## Set the generator object ##
        ##############################
        if generator is None:
            self.generator = self._raw_generator
        else:
            self.generator = generator
        ################################
        ## Set the animation function ##
        ################################
        if func is None:
            match global_analysis_method:
                case AnalysisMethod.CONTINUOUS_SCAN:
                    self._func = self._continuous_func
                case AnalysisMethod.FREQUENCY_MAP:
                    self._func = self._frequency_map_func
        else:
            self._func = func
        self.resize_interval = resize_interval
        self.resize_limit = self.resize_interval
        if fargs:
            self._args = fargs
        else:
            self._args = ()
        self._fig = fig
        # Disables blitting for backends that don't support it.  This
        # allows users to request it if available, but still have a
        # fallback that works if it is not.
        self._blit: bool = fig.canvas.supports_blit # type: ignore
        # Instead of starting the event source now, we connect to the figure's
        # draw_event, so that we only start once the figure has been drawn.
        self._first_draw_id = fig.canvas.mpl_connect("draw_event", self._start)
        # Connect to the figure's close_event so that we don't continue to
        # fire events and try to draw to a deleted figure.
        self._close_id = self._fig.canvas.mpl_connect("close_event", self._stop)
        self._setup_blit()
        #declarations for fields that will be initialized in other methods
        self._drawn_artists: List[matplotlib.artist.Artist]
        self._resize_id: int

    def _start(self, *args) -> None:
        """Start the animation.
        Add the draw frame command to the GUI; call show to start the event loop.
        """
        # First disconnect our draw event handler
        self._fig.canvas.mpl_disconnect(self._first_draw_id)
        #self._first_draw_id = None  # So we can check on save #but this is never checked
        # Now do any initial draw
        self._init_draw()
        class ThreadedAnimation(Thread):
            """Create a thread to obtain the file from a Queue
            and analyze the data.
            """
            def __init__(self, queue: Queue):
                global global_poison_pill
                Thread.__init__(self)     # initiate the thread
                self.queue = queue
                #-- set the poison pill event for Reset --#
                global_poison_pill = bool(tk.Event())
                self.file = 1
                root.after(10, self.start)                       # initiate the run() method
            def run(self):
                while True:
                    try:
                        task = self.queue.get(block=False)
                    except Empty:
                        break
                    else:
                        if not global_poison_pill:
                            root.after(global_search_interval, task)
                if not global_analysis_complete:
                    if not global_poison_pill:
                        root.after(10, self.run)
        ThreadedAnimation(queue=global_queue)
        self._step()

    def _stop(self, *args) -> None:
        """Stop the animation."""
        # On stop we disconnect all of our events.
        self._fig.canvas.mpl_disconnect(self._resize_id)
        self._fig.canvas.mpl_disconnect(self._close_id)

    def _setup_blit(self) -> None:
        """Set up the blit."""
        # Setting up the blit requires a cache of the background for the axes
        self._blit_cache: Dict[matplotlib.axes.Axes,\
                               Any] = {}
                                #2023-06-05: pylance lacks type information for matplotlib.backends
                                #on Linux, this is matplotlib.backends._backend_agg.BufferRegion
        self._drawn_artists = []
        self._resize_id = self._fig.canvas.mpl_connect("resize_event",
                                                       self._handle_resize)
        self._post_draw(True)

    def _blit_clear(self, artists: List[matplotlib.artist.Artist],\
                    bg_cache: Dict[matplotlib.axes.Axes,\
                                   Any]) -> None:
        """Get a list of the axes that need clearing from the artists that
        have been drawn. Grab the appropriate saved background from the
        cache and restore.
        """
        axs: Set[matplotlib.axes.Axes] = {artist.axes for artist in artists}
        for ax in axs:
            if ax in bg_cache:
                canvas = ax.figure.canvas
                canvas.restore_region(bg_cache[ax]) #type: ignore
                #2023-06-05: pylance lacks type information for matplotlib.backends

    def _init_draw(self) -> None:
        """Initialize the drawing by returning a sequence of blank artists."""
        self._drawn_artists = list(global_empty_plots)
        for artist in self._drawn_artists:
            artist.set_animated(self._blit)

    def _redraw_figures(self) -> None:
        """Resize raw plots."""
        _, axs = global_figures[self.num]
        for count in range(len(global_frequency_list)):
            match global_x_axis_mode:
                case PlotTimeReportingMode.EXPERIMENT_TIME:
                    axs[1][count].set_xlim(0, self.resize_limit*global_sample_rate/3600\
                                                + global_sample_rate/7200)
                    axs[2][count].set_xlim(0, self.resize_limit*global_sample_rate/3600\
                                                + global_sample_rate/7200)
                case PlotTimeReportingMode.FILE_NUMBER:
                    axs[1][count].set_xlim(0, self.resize_limit+0.1)
                    axs[2][count].set_xlim(0, self.resize_limit+0.1)
        self._post_draw(True)

    def _handle_resize(self, *args) -> None:
        """This is called when the window is resized."""
        # On resize, we need to disable the resize event handling so we don't
        # get too many events. Also stop the animation events, so that
        # we're paused. Reset the cache and re-init. Set up an event handler
        # to catch once the draw has actually taken place.
        #################################################
        ### Stop the event source and clear the cache ###
        #################################################
        self._fig.canvas.mpl_disconnect(self._resize_id)
        self._blit_cache.clear()
        self._init_draw()
        self._resize_id = self._fig.canvas.mpl_connect("draw_event",
                                                       self._end_redraw)

    def _end_redraw(self, event) -> None:
        # Now that the redraw has happened, do the post draw flushing and
        # blit handling. Then re-enable all of the original events.
        self._post_draw(True)
        self._fig.canvas.mpl_disconnect(self._resize_id)
        self._resize_id = self._fig.canvas.mpl_connect("resize_event",
                                                       self._handle_resize)

    def _draw_next_frame(self, framedata, fargs=None) -> None:
        # Breaks down the drawing of the next frame into steps of pre- and
        # post- draw, as well as the drawing of the frame itself.
        self._pre_draw(framedata)
        self._draw_frame(framedata, fargs)
        self._post_draw(False)

    def _pre_draw(self, framedata) -> None:
        # Perform any cleaning or whatnot before the drawing of the frame.
        # This default implementation allows blit to clear the frame.
        self._blit_clear(self._drawn_artists, self._blit_cache)

    def _draw_frame(self,
                    framedata: Tuple[List[float],\
                                    List[float],\
                                    List[float],\
                                    List[float],\
                                    List[float]],\
                    fargs) -> None:
        """Retrieve the data from _ratiometric_animation and blit the data onto the canvas."""
        # Ratiometric #
        if fargs:
            if fargs == "ratiometric_analysis":
                self._drawn_artists = self._ratiometric_animation(framedata, *self._args)
                self._drawn_artists = sorted(self._drawn_artists,
                                             key=lambda x: x.get_zorder())
                for artist in self._drawn_artists:
                    artist.set_animated(self._blit)
        else:
            self._drawn_artists = self._func(framedata, *self._args)
            if self._drawn_artists is None:
                raise RuntimeError("The animation function must return a "
                                   "sequence of Artist objects.")
            self._drawn_artists = sorted(self._drawn_artists,
                                         key=lambda x: x.get_zorder())
            for artist in self._drawn_artists:
                artist.set_animated(self._blit)

    def _post_draw(self, redraw: bool) -> None:
        # After the frame is rendered, this handles the actual flushing of
        # the draw, which can be a direct draw_idle() or make use of the
        # blitting.
        if redraw:
            # Data plots #
            self._fig.canvas.draw()
            # ratiometric plots
            match global_analysis_method:
                case AnalysisMethod.CONTINUOUS_SCAN:
                    if len(global_frequency_list) > 1:
                        ratio_fig, _ = global_ratiometric_figures[self.num]
                        ratio_fig.canvas.draw()
                case AnalysisMethod.FREQUENCY_MAP:
                    pass
        elif self._drawn_artists:
            self._blit_draw(self._drawn_artists, self._blit_cache)

    # The rest of the code in this class is to facilitate easy blitting
    def _blit_draw(self,\
                   artists: List[matplotlib.artist.Artist],\
                    bg_cache: Dict[matplotlib.axes.Axes,\
                                   Any]) -> None:
                                #2023-06-05: pylance lacks type information for matplotlib.backends
                                # on Linux, this is a matplotlib.backends._backend_agg.BufferRegion
        # Handles blitted drawing, which renders only the artists given instead
        # of the entire figure.
        updated_axs: List[matplotlib.axes.Axes] = []
        for artist in artists:
            # If we haven't cached the background for this axes object, do
            # so now. This might not always be reliable, but it's an attempt
            # to automate the process.
            if artist.axes not in bg_cache:
                canvas = artist.figure.canvas #type: ignore
                #2023-06-05: pylance lacks type information for
                # matplotlib.backends.backend_tkagg.FigureCanvasTkAgg
                bg_cache[artist.axes] = canvas.copy_from_bbox(artist.axes.bbox) #type: ignore
                #here artist is either a matplotlib.lines.Line2D or a matplotlib.patches.Polygon
            artist.axes.draw_artist(artist)
            updated_axs.append(artist.axes)
        # After rendering all the needed artists, blit each axes individually.
        for ax in set(updated_axs):
            ax.figure.canvas.blit(ax.bbox) #type: ignore
            #2023-06-05: pylance lacks type information for matplotlib.axes.Axes;
            # the actual type of ax.bbox is matplotlib.transforms.TransformedBbox

    ## callback that is called every 'interval' ms ##
    def _step(self) -> None:
        """Perform the next iteration of the animation."""
        
        frequency = int(global_frequency_list[self.count])
        match global_analysis_method:
            case AnalysisMethod.CONTINUOUS_SCAN:
                self.electrode = global_electrode_list[self.num]
                filename = make_file_name(self.file, self.electrode, frequency)
                myfile = global_file_path + filename
            case AnalysisMethod.FREQUENCY_MAP:
                filename = make_file_name(self.file, self.electrode, frequency)
                myfile = global_file_path + filename
        if file_is_complete(myfile):
            if self.file not in self.file_list:
                self.file_list.append(self.file)
                self.sample_list.append(get_time(len(self.file_list)))
            global_queue.put(lambda: self._run_analysis(myfile, frequency))
        else:
            if not global_poison_pill:
                print("Waiting for the next File")
                root.after(100, self._step)

    def _run_analysis(self, myfile: str, frequency: int) -> None:
        """Perform the next iteration of the data analysis."""
        try:
            framedata = self.generator(myfile, frequency)
            self._draw_next_frame(framedata)
            match global_analysis_method:
                case AnalysisMethod.CONTINUOUS_SCAN:
                    pass
                case AnalysisMethod.FREQUENCY_MAP:
                    global_track.tracking(self.file, frequency)
        except StopIteration:
            return None
        ##########################################################################
        ### if the resize limit has been reached, resize and redraw the figure ###
        ##########################################################################
        if self.file == self.resize_limit:
            # Dont redraw if this is already the last file #
            if self.resize_limit < global_number_of_files_to_process:
                ###############################################################
                ### If this is the last frequency, move onto the next limit ###
                ###############################################################
                if self.count == self.frequency_limit:
                    self.resize_limit = self.resize_limit + self.resize_interval
                    ### If the resize limit is above the number of files (e.g.
                    ### going out of bounds for the last resize event) then
                    ### readjust the final interval to the number of files
                    if self.resize_limit >= global_number_of_files_to_process:
                        self.resize_limit = global_number_of_files_to_process
                ############################################################
                ### 'if' statement used to make sure the plots are not   ###
                ### erased when there are no more files to be visualized ###
                ############################################################
                self._redraw_figures()
                try:
                    self._redraw_figures()
                except Exception as exception:
                    print("_run_analysis: Could not redraw figure", exception)
        ##################################################################
        ### If the function has analyzed each frequency for this file, ###
        ### move on to the next file and reset the frequency index     ###
        ##################################################################
        if self.count == self.frequency_limit:
            ######################################################
            ### If there are multiple frequencies, perform     ###
            ### ratiometric analysis and visualize the data    ###
            ######################################################
            match global_analysis_method:
                case AnalysisMethod.CONTINUOUS_SCAN:
                    if len(global_frequency_list) > 1:
                        try:
                            framedata = self._ratiometric_generator()
                            self._draw_next_frame(framedata, fargs="ratiometric_analysis")
                        except StopIteration:
                            return None
                    global_track.tracking(self.file, None)
                case AnalysisMethod.FREQUENCY_MAP:
                    pass
            #########################################################################
            ### If the function has analyzed the final file, remove the callback ###
            #########################################################################
            if self.file == global_number_of_files_to_process:
                match global_analysis_method:
                    case AnalysisMethod.CONTINUOUS_SCAN:
                        global_post_analysis._analysis_finished()
                    case AnalysisMethod.FREQUENCY_MAP:
                        pass
            else:
                self.file += 1
                self.index += 1
                self.count = 0
                root.after(1, self._step)
        ##########################################################
        ### Elif the function has not analyzed each frequency  ###
        ### for this file, move onto the next frequency        ###
        ##########################################################
        elif self.count < self.frequency_limit:
            self.count += 1
            root.after(1, self._step)

    def _raw_generator(self, myfile: str, frequency: int) ->\
          Tuple[List[float], List[float], List[float], List[float], List[float], List[float]]:
        """Generate data for raw data visualization."""
        ################################################
        ### Polynomial Regression or Gauss Range (V) ###
        ################################################
        xstart = global_xstart[self.num][self.count]
        xend = global_xend[self.num][self.count]
        gauss_peak = global_gauss_peak[self.num][self.count]
        gauss_baseline = global_gauss_baseline[self.num][self.count]
        gauss_maxheight = global_gauss_maxheight[self.num][self.count]
        ###################################
        ### Retrieve data from the File ###
        ###################################
        potentials: List[float]
        currents: List[float]
        potential_to_current_map: Dict[float, float]
        smooth_currents: List[float]
        adjusted_currents: List[float]
        potentials, currents, potential_to_current_map = read_data(myfile, self.electrode)
        cut_value = 0
        for value in potentials:
            if value == 0:
                cut_value += 1
        if cut_value > 0:
            potentials = potentials[:-cut_value]
            currents = currents[:-cut_value]
        ################################################################
        ### Adjust the potentials depending on user-input parameters ###
        ################################################################
        adjusted_potentials = [value for value in potentials if xend <= value <= xstart]
        match global_peak_method:
            case PeakMethod.POLY:
                #########################################
                ### Savitzky-Golay Smoothing          ###
                #########################################
                min_potential = min(potentials)            # find the min potential
                sg_limit = global_savitzky_golay_window/1000                  # mV --> V
                # shift all values positive
                sg_potentials = [x - min_potential for x in potentials]
                # find how many points fit within the sg potential window
                # this will be how many points are included in the rolling average
                sg_range = len([x for x in sg_potentials if x <= sg_limit])
                #--- Savitzky-golay Window must be greater than the range ---#
                if sg_range <= global_savitzky_golay_degree:
                    sg_range = global_savitzky_golay_degree + 1
                #-- if the range is even, make it odd --#
                if sg_range % 2 == 0:
                    sg_range = sg_range + 1
                # Apply the smoothing function and create a dictionary pairing
                # each potential with its corresponding current
                try:
                    smooth_currents = list(savgol_filter(currents, sg_range, global_savitzky_golay_degree))
                    potential_to_current_map = dict(zip(potentials, smooth_currents))
                except ValueError:
                    smooth_currents = list(savgol_filter(currents, 15, global_savitzky_golay_degree))
                    potential_to_current_map = dict(zip(potentials, smooth_currents))
                #######################################
                ### adjust the smooth currents to   ###
                ### match the adjusted potentials   ###
                #######################################
                adjusted_currents = []
                for potential in adjusted_potentials:
                    adjusted_currents.append(potential_to_current_map[potential])
                ######################
                ### Polynomial fit ###
                ######################
                polynomial_coeffs = np.polyfit(adjusted_potentials, adjusted_currents, POLYFIT_DEGREE)
                #############################
                ### Polynomial Regression ###
                #############################
                eval_regress: List[float] = np.polyval(polynomial_coeffs, adjusted_potentials).tolist()
                center_regress: List[float] = eval_regress # Copy of polynomial fit because reserved for Gauss Method
                ###############################################
                ### Absolute Max/Min Peak Height Extraction ###
                ###############################################
                #-- If the user selects 'Absolute Max/Min' in the 'Peak Height Extraction Settings'
                #-- within the Settings toolbar this analysis method will be used for PHE
                fit_half = round(len(eval_regress)/2)
                min1 = min(eval_regress[:fit_half])
                min2 = min(eval_regress[fit_half:])
                max1 = max(eval_regress[:fit_half])
                max2 = max(eval_regress[fit_half:])
                data: float
                ################################################################
                ### If the user selected Peak Height Extraction, analyze PHE ###
                ################################################################
                peak_height = max(max1, max2) - min(min1, min2)
            case PeakMethod.GAUSS:
                ####################
                ### Gauss Method ###
                ####################
                #######################################
                ### adjust the currents to          ###
                ### match the adjusted potentials   ###
                #######################################
                adjusted_currents = []
                for potential in adjusted_potentials:
                    adjusted_currents.append(potential_to_current_map[potential])
                ####################################################
                ### Create variables and solver for Gauss Method ###
                ####################################################
                data_y = adjusted_currents
                data_x = adjusted_potentials
                length = len(data_x)
                paramVal = np.concatenate((data_y , data_x))
                ########################################
                ### Choose solver based on frequency ###
                ########################################
                solver = global_gauss_solver[self.count]
                ##############################################################
                ### Initial values for optimization parameters             ###
                ### (height,mean,1/(2*variance),baseline,noise,abs(noise)) ###
                ##############################################################
                initstd = np.array([0.05, 0.0913, 0.0913])
                initlambda = 1.0 / (2 * initstd ** 2)
                var0 = np.concatenate(([0.8, 0.3, 0.3], [gauss_peak, -0.7, 0.1], initlambda, [gauss_baseline], np.zeros(length), 1e-2 * np.ones(length)))
                lboundstd = np.array([0.0289, 0.0707, 0.0707])
                uboundstd = np.array([0.1, 0.1291, 0.1291])
                uboundlambda = 1.0 / (2 * lboundstd ** 2)
                lboundlambda = 1.0 / (2 * uboundstd ** 2)
                ##########################################################################
                ### Upper and lower bounds for optimization parameters and constraints ###
                ##########################################################################
                ubx = np.concatenate((
                    [gauss_maxheight, 2*gauss_maxheight, 2*gauss_maxheight], # heights
                    [var0[3] + 1 * np.sqrt(1 / (2 * var0[6])), np.min(data_x), np.max(data_x) + 1 * np.sqrt(1 / (2 * var0[8]))], # means
                    uboundlambda, # 1/(2*variances)
                    [var0[9]+0.5], # baseline
                    np.ones(length), # noise
                    np.ones(length) # abs(noise)
                ))
                lbx = np.concatenate((
                    [1e-5, 0, 0], # heights
                    [var0[3] - 1 * np.sqrt(1 / (2 * var0[6])), np.min(data_x) - 1 * np.sqrt(1 / (2 * var0[7])), np.max(data_x)], # means
                    lboundlambda, # 1/(2*variances)
                    [max(var0[9]-0.5,0)], # baseline
                    -1.0 * np.ones(length), # noise
                    -1.0 * np.ones(length) # abs(noise)
                ))
                ubg = np.concatenate((
                    np.zeros(2),
                    np.zeros(length), # noise - abs(noise)
                    np.inf * np.ones(length), # noise + abs(noise)
                    [1] # sum(abs(noise))
                ))
                lbg = np.concatenate((
                    -np.inf * np.ones(2),
                    -np.inf * np.ones(length), # noise - abs(noise)
                    np.zeros(length),  # noise + abs(noise)
                    [0] # sum(abs(noise))
                ))
                # Solve the problem
                solution = solver(
                    x0=var0,
                    p=paramVal,
                    ubx=ubx,
                    lbx=lbx,
                    ubg=ubg,
                    lbg=lbg
                )
                theta = solution['x'].full().flatten()
                # Calculate fit
                fit = (theta[9] + theta[0] * np.exp(-(data_x - theta[3])**2 * theta[6]) +
                    theta[1] * np.exp(-(data_x - theta[4])**2 * theta[7]) +
                    theta[2] * np.exp(-(data_x - theta[5])**2 * theta[8]))
                fit_center = theta[9] + theta[0] * np.exp(-(data_x - theta[3])**2 * theta[6])               
                # Regression output for Gauss Method 
                eval_regress: List[float] = fit.tolist()
                # Regression output for Gauss Method (Only Center)
                center_regress: List[float] = fit_center.tolist()
                # There is no smoothing in Gauss method
                smooth_currents: List[float] = currents
                ################################################################
                ### If the user selected Peak Height Extraction, analyze PHE ###
                ################################################################
                data: float
                peak_height=theta[0].item()
                peak_loc = theta[3].item()
                global_peak_list[self.num][self.count][self.index] = peak_loc

        match global_plot_summary_mode:
            case PlotSummaryMode.PHE:
                data = peak_height
            case PlotSummaryMode.AUC:
                ##################################
                ### Integrate Area Under the   ###
                ### Curve using a Riemann Sum  ###
                ##################################
                auc_index = 1
                auc: float = 0
                auc_potentials = adjusted_potentials
                match global_peak_method:
                    case PeakMethod.POLY: 
                        #--- Find the minimum value and normalize it to 0 ---#
                        auc_min = min(adjusted_currents)
                        auc_currents = [y - auc_min for y in adjusted_currents]
                    case PeakMethod.GAUSS:
                        #--- Find the minimum value and normalize it to 0 ---# 
                        auc_min = min(eval_regress)
                        auc_currents = [y - auc_min for y in eval_regress]
                #--- Midpoint Riemann Sum ---#
                while auc_index <= len(auc_currents) - 1:
                    auc_height = (auc_currents[auc_index] + auc_currents[auc_index - 1])/2
                    auc_width = auc_potentials[auc_index] - auc_potentials[auc_index - 1]
                    auc += (auc_height * auc_width)
                    auc_index += 1
                data = auc
        #######################################
        ### Save the data into global lists ###
        #######################################
        global_data_list[self.num][self.count][self.index] = data
        match global_analysis_method:
            case AnalysisMethod.CONTINUOUS_SCAN:
                global_data_normalization.normalize(self.file, data,\
                                                    self.num, self.count, self.index)
            case AnalysisMethod.FREQUENCY_MAP:
                local_frequency = global_frequency_list[self.count]
                self.frequency_axis.append(local_frequency)
                self.charge_axis.append(peak_height/local_frequency)
        #####################################################
        ### Return data to the animate function as 'args' ###
        #####################################################
        return potentials, adjusted_potentials, smooth_currents, adjusted_currents, eval_regress, center_regress

    def _continuous_func(self,\
                         framedata: Tuple[List[float],\
                                        List[float],\
                                        List[float],\
                                        List[float],\
                                        List[float],\
                                        List[float]],\
                        *args) -> List[matplotlib.artist.Artist]:
        """Generate plots for continuous scan visualization."""
        potentials: List[float]
        adjusted_potentials: List[float]
        smooth_currents: List[float]
        adjusted_currents: List[float]
        regression: List[float]
        center_regress: List[float]
        x_axis: Sequence[float]
        if global_key > 0:
            while True:
                potentials,\
                    adjusted_potentials,\
                    smooth_currents,\
                    adjusted_currents,\
                    regression,\
                    center_regress =\
                        framedata
                #############################################################
                ### Acquire the current frequency and get the xstart/xend ###
                ### parameters that will manipulate the visualized data   ###
                #############################################################
                ###################################
                ### Set the units of the X-axis ###
                ###################################
                match global_x_axis_mode:
                    case PlotTimeReportingMode.EXPERIMENT_TIME:
                        x_axis = self.sample_list
                    case PlotTimeReportingMode.FILE_NUMBER:
                        x_axis = self.file_list
                ################################################################
                ### Acquire the artists for this electrode at this frequency ###
                ### and get the data that will be visualized                 ###
                ################################################################
                # 'num' is the electrode index value
                # 'count' is the frequency index value
                #plots: List[matplotlib.artist.Artist] =\
                    #global_plot_list_continuous_scan[self.num][self.count]
                (smooth,\
                    regress,\
                    peak,\
                    peak_injection,\
                    normalization,\
                    norm_injection,\
                    poly,\
                    center) = global_plot_list_continuous_scan[self.num][self.count]
                ##########################
                ### Visualize the data ###
                ##########################
                #--- Peak Height ---#
                data = global_data_list[self.num][self.count][:len(self.file_list)]
                if global_frequency_list[self.count] == global_high_low_dictionary[HighLow.LOW]:
                    local_normalized_data_list =\
                          global_offset_normalized_data_list[self.num][:len(self.file_list)]
                else:
                    local_normalized_data_list =\
                          global_normalized_data_list[self.num][self.count][:len(self.file_list)]
                ####################################################
                ### Set the data of the artists to be visualized ###
                ####################################################
                if global_injection_point is None:
                    #Smooth current voltammogram:
                    smooth.set_data(potentials, smooth_currents)
                    regress.set_data(adjusted_potentials, regression)
                    center.set_data(adjusted_potentials,center_regress)
                    #Raw Data:
                    peak.set_data(x_axis, data)
                    #Normalized Data:
                    normalization.set_data(x_axis, local_normalized_data_list)
                ##########################################################
                ### If an Injection Point has been set, visualize the  ###
                ### data before and after the injection separately     ###
                ##########################################################
                else:
                    if self.file >= global_injection_point:
                        injection_index = global_injection_point - 1
                        ####################################################
                        ### Set the data of the artists to be visualized ###
                        ####################################################
                        #Smooth current voltammogram:
                        smooth.set_data(potentials, smooth_currents)
                        #Regression voltammogram:
                        regress.set_data(adjusted_potentials, regression)
                        center.set_data(adjusted_potentials,center_regress)
                        #Raw Data up until injection point:
                        peak.set_data(x_axis[:injection_index], data[:injection_index])
                        #Raw Data after injection point:
                        peak_injection.set_data(x_axis[injection_index:], data[injection_index:])
                        #Norm Data up until injection point:
                        normalization.set_data(x_axis[:injection_index],\
                                          local_normalized_data_list[:injection_index])
                        #Norm Data after injection point:
                        norm_injection.set_data(x_axis[injection_index:],\
                                          local_normalized_data_list[injection_index:])
                    else:
                        # Smooth current voltammogram:
                        smooth.set_data(potentials, smooth_currents)
                        regress.set_data(adjusted_potentials, regression)
                        center.set_data(adjusted_potentials,center_regress)
                        # Raw Data:
                        peak.set_data(x_axis, data)
                        # Clear the injection artist:
                        peak_injection.set_data([], [])
                        # Norm Data:
                        normalization.set_data(x_axis, local_normalized_data_list)
                        #Clear the injection artist:
                        norm_injection.set_data([], [])
                match global_plot_summary_mode:
                    case PlotSummaryMode.PHE:
                        pass
                    case PlotSummaryMode.AUC:
                        #--- Shaded region of Area Under the Curve ---#
                        vertices = [(adjusted_potentials[0], adjusted_currents[0]),\
                                        *zip(adjusted_potentials, adjusted_currents),\
                                            (adjusted_potentials[-1], adjusted_currents[-1])]
                        poly.set_xy(vertices)
                return [smooth,\
                    regress,\
                    peak,\
                    peak_injection,\
                    normalization,\
                    norm_injection,\
                    poly,\
                    center]
        else:
            time.sleep(0.1)
            print("Yielding Empty Plots in Animation")
            return []

    def _frequency_map_func(self,\
                            framedata: Tuple[List[float],\
                                            List[float],\
                                            List[float],\
                                            List[float],\
                                            List[float]],
                            *args) -> List[matplotlib.artist.Artist]:
        """Generate plots for frequency map visualization."""
        if global_key > 0:
            while True:
                potentials, adjusted_potentials, smooth_currents, _, regression = framedata
                ################################################################
                ### Acquire the artists for this electrode at this frequency ###
                ### and get the data that will be visualized                 ###
                ################################################################
                smooth: matplotlib.lines.Line2D
                regress: matplotlib.lines.Line2D
                charge: matplotlib.lines.Line2D
                (smooth, regress, charge) = global_plot_list_frequency_map[self.num]
                ##########################
                ### Visualize the data ###
                ##########################
                #--- Peak Height ---#
                ####################################################
                ### Set the data of the artists to be visualized ###
                ####################################################
                smooth.set_data(potentials, smooth_currents)     # Smooth current voltammogram
                regress.set_data(adjusted_potentials, regression) # Regression voltammogram
                charge.set_data(self.frequency_axis, self.charge_axis)
                return [smooth, regress, charge]
        else:
            time.sleep(0.1)
            print("Yielding Empty Plots in Animation")
            return []

    ############################
    ### Ratiometric Analysis ###
    ############################
    def _ratiometric_generator(self):
        """Generate data for ratiometric analysis."""
        local_high_frequency = global_high_low_dictionary[HighLow.HIGH]
        high_count = global_frequency_dict[local_high_frequency]
        high_point = global_normalized_data_list[self.num][high_count][self.index]
        low_point = global_offset_normalized_data_list[self.num][self.index]
        normalized_ratio = high_point/low_point
        match global_kdm_method:
            case KDMMethod.OLD:
                kdm = (high_point - low_point) + 1
            case KDMMethod.NEW:
                average = 0.5*(high_point + low_point)
                kdm = (high_point - low_point)/average + 1
        #-- save the data to global lists --#
        global_normalized_ratiometric_data_list[self.num].append(normalized_ratio)
        global_kdm_list[self.num].append(kdm)
        return normalized_ratio, kdm

    def _ratiometric_animation(self, framedata, *args) -> List[matplotlib.artist.Artist]:
        """Generate plots for ratiometric analysis."""
        plots = global_ratiometric_plots[self.num]
        x_axis: List[float]
        match global_x_axis_mode:
            case PlotTimeReportingMode.EXPERIMENT_TIME:
                x_axis = self.sample_list
            case PlotTimeReportingMode.FILE_NUMBER:
                x_axis = [float(i) for i in self.file_list]
        local_norm = [x*100 for x in global_normalized_ratiometric_data_list[self.num]]
        kdm = [x*100 for x in global_kdm_list[self.num]]
        if len(global_frequency_list) > 1:
            ratio_fig, ratio_axs = global_ratiometric_figures[self.num]
            x_left_bound, x_right_bound = compute_x_bounds(self.file)
            ratio_axs[0].set_xlim(x_left_bound, x_right_bound)
            ratio_axs[1].set_xlim(x_left_bound, x_right_bound)
            match global_y_norm_radiobutton:
                case YBound.AUTOMATIC:
                    local_norm_within_x_window = []
                    for i in range(len(x_axis)):
                        if x_left_bound <= x_axis[i] <= x_right_bound:
                            local_norm_within_x_window.append(local_norm[i])
                    norm_mean: float = float(np.mean(local_norm_within_x_window))
                    norm_std: float = float(np.std(local_norm_within_x_window))
                    ylimit_norm_low = norm_mean - 3*norm_std
                    ylimit_norm_high = norm_mean + 3*norm_std
                case YBound.MANUAL:
                    ylimit_norm_low = 100*global_min_norm
                    ylimit_norm_high = 100*global_max_norm
            match global_y_kdm_radiobutton:
                case YBound.AUTOMATIC:
                    local_kdm_within_x_window = []
                    for i in range(len(x_axis)):
                        if x_left_bound <= x_axis[i] <= x_right_bound:
                            local_kdm_within_x_window.append(kdm[i])
                    kdm_mean: float = float(np.mean(local_kdm_within_x_window))
                    kdm_std: float = float(np.std(local_kdm_within_x_window))
                    ylimit_kdm_low = kdm_mean - 3*kdm_std
                    ylimit_kdm_high = kdm_mean + 3*kdm_std
                case YBound.MANUAL:
                    ylimit_kdm_low = 100*global_min_kdm
                    ylimit_kdm_high = 100*global_max_kdm
            ratio_axs[0].set_ylim(ylimit_norm_low, ylimit_norm_high)
            ratio_axs[1].set_ylim(ylimit_kdm_low, ylimit_kdm_high)
            ratio_fig.canvas.draw()
        ##########################################
        ## If an injection point has not been   ##
        ## chosen, visualize the data as usual  ##
        ##########################################
        if global_injection_point is None:
            plots[0].set_data(x_axis, local_norm)
            plots[2].set_data(x_axis, kdm)
        ############################################
        ## If an injection point has been chosen  ##
        ## chosen, visualize the injection        ##
        ## points separately                      ##
        ############################################
        else:
            #-- list index value for the injection point --#
            injection_index = global_injection_point - 1
            if self.file >= global_injection_point:
                #-- if the injection point has already been --#
                #-- analyzed, separate the visualized data  --#
                plots[0].set_data(x_axis[:injection_index], local_norm[:injection_index])
                plots[1].set_data(x_axis[injection_index:], local_norm[injection_index:])
                plots[2].set_data(x_axis[:injection_index], kdm[:injection_index])
                plots[3].set_data(x_axis[injection_index:], kdm[injection_index:])
            else:
                #-- if the file is below the injectionpoint, wait until  --#
                #-- the point is reached to visualize the injection data --#
                plots[0].set_data(x_axis, local_norm)
                plots[1].set_data([], [])
                plots[2].set_data(x_axis, kdm)
                plots[3].set_data([], [])
        plots_as_artists: List[matplotlib.artist.Artist] = list(plots)
        return plots_as_artists
                                        ##############################
                                        ##############################
                                        ### END OF ANIMATION CLASS ###
                                        ##############################
                                        ##############################
#-------------------------------------------------------------------------------------------------#
class DataNormalization():
    """Normalize the data to the first file, or perhaps another file."""
    def __init__(self):
        pass

    def normalize(self, file: int, data: float, num: int, count: int, index: int) -> None:
        """Normalize the data to the given file."""
        global global_initialized_normalization
        sample: float = get_time(len(global_file_list))
        #######################################################
        ## Check the frequency and apply the baseline offset ##
        #######################################################
        frequency: int = global_frequency_list[count]
        offset: float
        if frequency == global_high_low_dictionary[HighLow.LOW]:
            match global_x_axis_mode:
                case PlotTimeReportingMode.EXPERIMENT_TIME:
                    offset = sample*global_low_frequency_slope + global_low_frequency_offset
                case PlotTimeReportingMode.FILE_NUMBER:
                    offset = file*global_low_frequency_slope + global_low_frequency_offset
        else:
            offset = 0
        normalization_index: int = global_normalization_point - 1
        #--- If the file being used as the standard has been analyzed,
        # normalize the data to that point ---#
        if file >= global_normalization_point:
            if global_normalization_point not in global_normalization_vault:
                global_normalization_vault.append(global_normalization_point)
            #-- if the software has still been normalizing to the first file,
            # start normalizing to the normalization point --#
            global_initialized_normalization = True
            ###########################################################
            ### If the rest of the data has already been normalized ###
            ### to this point, continue to normalize the data for   ###
            ### the current file to the normalization point         ###
            ###########################################################
            global_normalized_data_list[num][count][index] =\
                  data/global_data_list[num][count][normalization_index]
            ###########################################################################
            ### If this is a low frequency, apply the offset to the normalized data ###
            ###########################################################################
            if frequency == global_high_low_dictionary[HighLow.LOW]:
                global_offset_normalized_data_list[num][index] =\
                      global_normalized_data_list[num][count][index] + offset
        #######################################################################
        ### Elif the chosen normalization point is greater than the current ###
        ### file, continue to normalize to the previous normalization point ###
        #######################################################################
        elif global_initialized_normalization:
            ### Acquire the normalization point that was previously selected ###
            temp_normalization_point: int = global_normalization_vault[-1]
            temp_normalization_index: int = temp_normalization_point - 1
            global_normalized_data_list[num][count][index] =\
                  data/global_data_list[num][count][temp_normalization_index]
            ###########################################################################
            ### If this is a low frequency, apply the offset to the normalized data ###
            ###########################################################################
            if global_frequency_list[count] == global_high_low_dictionary[HighLow.LOW]:
                global_offset_normalized_data_list[num][index] =\
                      global_normalized_data_list[num][count][index] + offset
        #--- Else, if the initial normalization point has not yet been reached,
        # normalize to the first file ---#
        else:
            global_normalized_data_list[num][count][index] = data/global_data_list[num][count][0]
            ###########################################################################
            ### If this is a low frequency, apply the offset to the normalized data ###
            ###########################################################################
            if global_frequency_list[count] == global_high_low_dictionary[HighLow.LOW]:
                global_offset_normalized_data_list[num][index] =\
                      global_normalized_data_list[num][count][index] + offset

    def renormalize_data(self, file: int) -> None:
        """Renormalize the data to the new normalization point."""
        ##############################################################
        ## If the normalization point equals the current file,      ##
        ## normalize all of the data to the new normalization point ##
        #############################################################
        index: int
        normalization_index: int
        offset: float
        if file == global_normalization_point:
            index = file - 1
            normalization_index = global_normalization_point - 1
            for num in range(global_electrode_count):
                for count in range(len(global_frequency_list)):
                    global_normalized_data_list[num][count][:index] =\
                          [(idx/global_data_list[num][count][normalization_index])\
                           for idx in global_data_list[num][count][:index]]
                    ##################################################
                    ## If the frequency is below cutoff_frequency, ###
                    ## add the baseline Offset                     ###
                    ##################################################
                    if global_frequency_list[count] == global_high_low_dictionary[HighLow.LOW]:
                        for index in range(len(global_file_list)):
                            ##########################
                            ## Calculate the offset ##
                            ##########################
                            sample = global_sample_list[index]
                            file = global_file_list[index]
                            match global_x_axis_mode:
                                case PlotTimeReportingMode.EXPERIMENT_TIME:
                                    offset = sample*global_low_frequency_slope\
                                            + global_low_frequency_offset
                                case PlotTimeReportingMode.FILE_NUMBER:
                                    offset = file*global_low_frequency_slope\
                                            + global_low_frequency_offset
                            global_offset_normalized_data_list[num][index] =\
                                  global_normalized_data_list[num][count][index] + offset
            ################################################
            ### Analyze KDM using new normalization data ###
            ################################################
            if len(global_frequency_list) > 1:
                self.reset_ratiometric_data()
            ###############################
            ### GUI Normalization Label ###
            ###############################
            global_norm_warning.config(foreground="green")
            global_norm_warning.config(text=f"Normalized to file {global_normalization_point}")
            ########################################################################
            ### If .txt file export has been activated, update the exported data ###
            ########################################################################
            if global_text_file_export_activated:
                global_text_file_export.txt_file_normalization()
        #########################################################################
        ## If the Normalization Point has been changed and the current file is ##
        ## greater than the new point, renormalize the data to the new point   ##
        #########################################################################
        if global_normalization_waiting:
            index = file - 1
            normalization_index = global_normalization_point - 1
            for num in range(global_electrode_count):
                for count in range(len(global_frequency_list)):
                    ##########################
                    ## Renormalize the data ##
                    ##########################
                    global_normalized_data_list[num][count][:index] =\
                          [idx/global_data_list[num][count][normalization_index]\
                           for idx in global_data_list[num][count][:index]]
                    ##################################################
                    ## If the frequency is below cutoff_frequency,  ##
                    ## add the baseline Offset                      ##
                    ##################################################
                    if global_frequency_list[count] == global_high_low_dictionary[HighLow.LOW]:
                        for index in range(len(global_file_list)):
                            ##########################
                            ## Calculate the offset ##
                            ##########################
                            sample = global_sample_list[index]
                            file = index + 1
                            match global_x_axis_mode:
                                case PlotTimeReportingMode.EXPERIMENT_TIME:
                                    offset = sample*global_low_frequency_slope\
                                            + global_low_frequency_offset
                                case PlotTimeReportingMode.FILE_NUMBER:
                                    offset = file*global_low_frequency_slope\
                                            + global_low_frequency_offset
                            global_offset_normalized_data_list[num][index] =\
                                  global_normalized_data_list[num][count][index] + offset
            ################################################
            ## Using the newly normalized data, calculate ##
            ## the Normalized Ratio and KDM               ##
            ## for each file that has been analyzed       ##
            ################################################
            if len(global_frequency_list) > 1:
                self.reset_ratiometric_data()
            #--- Indicate that the data has been normalized to the new normalization_point ---#
            global_norm_warning.config(foreground="green")
            global_norm_warning.config(text=f"Normalized to file {global_normalization_point}")
            global_wait_time.normalization_proceed()
            #-- if .txt file export has been activated, update the exported data ---#
            if global_text_file_export_activated:
                global_text_file_export.txt_file_normalization()

    #############################################################
    ### Readjust the data to the new user-inputted parameters ###
    #############################################################
    def reset_ratiometric_data(self) -> None:
        ############################################
        ### Readjust Low Frequencies with Offset ###
        ############################################
        #-- Iterate through every frequency --#
        for frequency in global_frequency_list:
            #-- Only apply the offset if the frequency is below cutoff_frequency --#
            if frequency == global_high_low_dictionary[HighLow.LOW]:
                count = global_frequency_dict[frequency]
                #-- Apply the offset to every file --#
                for index in range(len(global_file_list)):
                    sample = global_sample_list[index]
                    file = global_file_list[index]
                    offset: float
                    match global_x_axis_mode:
                        case PlotTimeReportingMode.EXPERIMENT_TIME:
                            offset = sample*global_low_frequency_slope\
                                    + global_low_frequency_offset
                        case PlotTimeReportingMode.FILE_NUMBER:
                            offset = file*global_low_frequency_slope\
                                    + global_low_frequency_offset
                    for num in range(global_electrode_count):
                        global_offset_normalized_data_list[num][index] =\
                              global_normalized_data_list[num][count][index] + offset
        ####################################################
        ### Readjust KDM with newly adjusted frequencies ###
        ####################################################
        for file in global_file_list:
            index = file - 1
            for num in range(global_electrode_count):
                # grab the index value for the current high and low frequencies
                # used for ratiometric analysis #
                high_count: int = global_frequency_dict[global_high_frequency]
                high_point: float = global_normalized_data_list[num][high_count][index]
                low_point: float = global_offset_normalized_data_list[num][index]
                normalized_data_ratio: float = high_point/low_point
                global_normalized_ratiometric_data_list[num][index] = normalized_data_ratio
                #-- KDM ---#
                match global_kdm_method:
                    case KDMMethod.OLD:
                        kdm: float = (high_point - low_point) + 1
                    case KDMMethod.NEW:
                        average: float = 0.5*(high_point + low_point)
                        kdm: float = (high_point - low_point)/average + 1
                global_kdm_list[num][index] = kdm
        #-- if .txt file export has been activated, update the exported data ---#
        if global_text_file_export_activated:
            if not global_analysis_complete:
                global_text_file_export.txt_file_normalization()

class PostAnalysis(ttk.Frame):
    """Post Analysis Module for data manipulation after the completion of data analysis."""
    def __init__(self, parent: ttk.Frame, container):
        global global_analysis_complete
        ############################
        ### Class-wide variables ###
        ############################
        ttk.Frame.__init__(self, parent)
        self.parent = parent
        self.container = container
        #-- global boolean to control activation of this class --#
        global_analysis_complete = False
        self.export_top_level_exists = False
        #-- once completion value == electrode_count, analysis_complete --#
        #-- will be changed from False to True                          --#
        self.completion_value = 0
        self.high = global_frequency_list[-1] > CUTOFF_FREQUENCY
        self.low = global_frequency_list[0] <= CUTOFF_FREQUENCY
        ##########################################
        ### Initialize the Post Analysis Frame ###
        ##########################################
        self._initialize_frame()
        #declarations for fields that will be initialized in other methods
        self.frequency_frame: ttk.Frame
        self.set_point_norm: ttk.Entry
        self.norm_warning: ttk.Label
        self.high_frequency_entry: ttk.Entry
        self.low_frequency_entry: ttk.Entry
        self.low_frequency_offset: ttk.Entry
        self.low_frequency_slope: ttk.Entry
        self.win: tk.Toplevel
        self.select_file_path: ttk.Button
        self.no_selected_path: ttk.Label
        self.path_warning_exists: bool
        self.filehandle: ttk.Entry
        self.electrode_label: ttk.Label
        self.electrode_count: tk.Listbox
        self.frequency_label: ttk.Label
        self.frequency_list: tk.Listbox
        self.electrode_list_exists: bool
        self.frequency_list_exists: bool

    def _initialize_frame(self) -> None:
        """Initialize the Post Analysis Frame."""
        ttk.Frame.__init__(self, self.parent)
        ttk.Label(self, text="Post Analysis", font=HUGE_FONT).grid(row=0, column=0, columnspan=2)
        data_adjustment_frame = ttk.Frame(self, relief="groove")
        data_adjustment_frame.grid(row=1, column=0, columnspan=2,\
                                   pady=5, ipadx=50, padx=2, sticky="ns")
        normalization_frame = ttk.Frame(data_adjustment_frame)
        normalization_frame.grid(row=1, column=0, pady=5)
        #--- Real-time Normalization Variable ---#
        ttk.Label(normalization_frame, text="Set Normalization Point", font=MEDIUM_FONT).\
            grid(row=0, column=0, pady=5)
        normalization_var = tk.StringVar()
        norm_string = str(global_normalization_point)
        normalization_var.set(norm_string)
        self.set_point_norm = ttk.Entry(normalization_frame, textvariable=normalization_var,\
                                         width=8)
        self.set_point_norm.grid(row=1, column=0, pady=5)
        #--- Button to apply any changes to the normalization variable ---#
        ttk.Button(normalization_frame, text="Apply Norm",\
                   command=self.post_analysis_normalization, width=10).\
                    grid(row=2, column=0)
        self.norm_warning = ttk.Label(normalization_frame, text="",\
                                      foreground="red", font=MEDIUM_FONT)
        if len(global_frequency_list) > 1:
            self.frequency_frame = ttk.Frame(data_adjustment_frame, relief="groove")
            self.frequency_frame.grid(row=2, column=0, pady=10, padx=3, ipady=2)
            #--- Drift Correction Title ---#
            ttk.Label(self.frequency_frame, text="Drift Correction", font=LARGE_FONT).\
                grid(row=0, column=0, columnspan=3, pady=1, padx=5)
            #--- High Frequency Selection for KDM and Ratiometric Analysis ---#
            ttk.Label(self.frequency_frame, text="High Frequency", font=MEDIUM_FONT).\
                        grid(row=1, column=1, pady=5, padx=5)
            self.high_frequency_entry = ttk.Entry(self.frequency_frame, width=7)
            self.high_frequency_entry.insert(tk.END, str(global_high_frequency))
            self.high_frequency_entry.grid(row=2, column=1, padx=5)
            #--- Low Frequency Selection for KDM and Ratiometric Analysis ---#
            ttk.Label(self.frequency_frame, text="Low Frequency", font=MEDIUM_FONT).\
                grid(row=1, column=0, pady=5, padx=5)
            self.low_frequency_entry = ttk.Entry(self.frequency_frame, width=7)
            self.low_frequency_entry.insert(tk.END, str(global_low_frequency))
            self.low_frequency_entry.grid(row=2, column=0, padx=5)
            tk.Label(self.frequency_frame, text="Low Frequency\n Offset",\
                     font=MEDIUM_FONT).\
                        grid(row=3, column=0, pady=2, padx=2)
            self.low_frequency_offset = ttk.Entry(self.frequency_frame, width=7)
            self.low_frequency_offset.insert(tk.END, str(global_low_frequency_offset))
            self.low_frequency_offset.grid(row=4, column=0, padx=2, pady=2)
            ttk.Label(self.frequency_frame, text="Low Frequency\n Slope Manipulation",\
                      font=MEDIUM_FONT).\
                        grid(row=3, column=1, pady=2, padx=2)
            self.low_frequency_slope = ttk.Entry(self.frequency_frame, width=7)
            self.low_frequency_slope.insert(tk.END, str(global_low_frequency_slope))
            self.low_frequency_slope.grid(row=4, column=1, padx=2, pady=2)
            ttk.Button(self.frequency_frame, text="Apply Frequencies",\
                       command=self.post_analysis_kdm).\
                        grid(row=5, column=0, columnspan=2, pady=5, padx=5)
        ttk.Button(data_adjustment_frame, text="Redraw Figures",\
                   command=self._draw, width=12).\
                    grid(row=3, column=0, pady=7)
        data_adjustment_frame.columnconfigure(0, weight=1)
        row_value = 3
        self.data_export_frame = ttk.Frame(self, relief="groove")
        self.data_export_frame.grid(row=row_value, column=0, pady=5, ipady=5)
        ttk.Button(self.data_export_frame, text="Data Export Settings",\
                  command=lambda: self.data_export_top_level)
        #---Buttons to switch between electrode frames---#
        frame_value = 0
        column_value = 0
        for _ in global_plot_values:
            ttk.Button(self, text=global_frame_list[frame_value],\
                       command=lambda frame_value=frame_value:\
                        self.show_plot(global_plot_values[frame_value])).\
                            grid(row=row_value, column=column_value, pady=2, padx=5)
            ## allows .grid() to alternate between
            ## packing into column 1 and column 2
            if column_value == 1:
                column_value = 0
                row_value += 1
            ## if gridding into the 1st column,
            ## grid the next into the 2nd column
            else:
                column_value += 1
            frame_value += 1
        row_value += 1
        export_settings = ttk.Frame(self)
        export_settings.grid(row=row_value, column=0, columnspan=2, pady=5, ipady=10)
        ttk.Button(export_settings, text="Post-analysis data export",\
                  command=self.data_export_top_level).\
                    grid(row=0, column=0, padx=5)
        export_settings.columnconfigure(1, weight=1)
        export_settings.rowconfigure(1, weight=1)
        row_value += 1
        #--- Reset ---#
        ttk.Button(self, text="Reset", style="Fun.TButton", command=self._reset).\
            grid(row=row_value, column=1, pady=5, padx=5)
        #--- Quit ---#
        ttk.Button(self, text="Quit", command=lambda: os._exit(0)).\
            grid(row=row_value, column=0, pady=5)
        for row in range(row_value):
            self.rowconfigure(row, weight=1)
        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)

    def _analysis_finished(self) -> None:
        global global_analysis_complete
        self.completion_value += 1
        if self.completion_value == global_electrode_count:
            global_analysis_complete = True
            #####################################
            ### Raise the Post Analysis Frame ###
            #####################################
            global_show_frames[FRAME_POST_ANALYSIS].tkraise()

    def _adjust_data(self) -> None:
        ###################################
        ### Renormalize all of the data ###
        ###################################
        normalization_index = global_normalization_point - 1
        if global_normalization_point <= global_number_of_files_to_process:
            for num in range(global_electrode_count):
                for count in range(len(global_frequency_list)):
                    global_normalized_data_list[num][count] =\
                        [(idx/global_data_list[num][count][normalization_index])\
                         for idx in global_data_list[num][count]]
                    ##################################################
                    ## If the frequency is below cutoff_frequency,  ##
                    ## add the baseline Offset                      ##
                    ##################################################
                    if global_frequency_list[count] == global_high_low_dictionary[HighLow.LOW]:
                        for index in range(global_number_of_files_to_process):
                            ##########################
                            ## Calculate the offset ##
                            ##########################
                            sample = global_sample_list[index]
                            file = global_file_list[index]
                            match global_x_axis_mode:
                                case PlotTimeReportingMode.EXPERIMENT_TIME:
                                    offset = (sample*global_low_frequency_slope)\
                                            + global_low_frequency_offset
                                case PlotTimeReportingMode.FILE_NUMBER:
                                    offset = (file*global_low_frequency_slope)\
                                            + global_low_frequency_offset
                            global_offset_normalized_data_list[num][index] =\
                                global_normalized_data_list[num][count][index] + offset
        global_data_normalization.reset_ratiometric_data()
        self.norm_warning.config(foreground="green")
        self.norm_warning.config(text=f"Normalized to File %{global_normalization_point}")
        if global_text_file_export_activated:
            global_text_file_export.txt_file_normalization()
        ### Draw the readjusted data
        self._draw()

    def _draw(self) -> None:
        x_axis: Sequence[float]
        match global_x_axis_mode:
            case PlotTimeReportingMode.EXPERIMENT_TIME:
                x_axis = global_sample_list
            case PlotTimeReportingMode.FILE_NUMBER:
                x_axis = global_file_list
        for num in range(global_electrode_count):
            ## get the figure for the electrode ##
            fig, axs = global_figures[num]
            for count in range(len(global_frequency_list)):
                ###################################
                ### Set the units of the X-axis ###
                ###################################
                ################################################################
                ### Acquire the artists for this electrode at this frequency ###
                ### and get the data that will be visualized                 ###
                ################################################################
                ##########################
                ### Visualize the data ###
                ##########################
                #--- Peak Height ---#
                if global_frequency_list[count] == global_high_low_dictionary[HighLow.LOW]:
                    local_normalized_data_list = global_offset_normalized_data_list[num]
                else:
                    local_normalized_data_list = global_normalized_data_list[num][count]
                ### Draw new data ###
                axs[1][count].clear()
                axs[2][count].clear()
                supposedly_global_norm = axs[2][count].\
                    plot(x_axis, local_normalized_data_list, "ko", markersize=1)
                #####################
                ## Set the Y Label ##
                #####################
                axs[0][0].set_ylabel("Current\n(µA)", fontweight="bold")
                match global_plot_summary_mode:
                    case PlotSummaryMode.PHE:
                        axs[1][0].set_ylabel("Peak Height\n(µA)", fontweight="bold")
                    case PlotSummaryMode.AUC:
                        axs[1][0].set_ylabel("AUC (a.u.)", fontweight="bold")
                axs[2][0].set_ylabel("Normalized", fontweight="bold")
            fig.canvas.draw_idle()
            ### If necessary, redraw ratiometric data ###
            if len(global_frequency_list) > 1:
                ratio_fig, ratio_axs = global_ratiometric_figures[num]
                supposedly_global_norm = [x*100\
                                            for x in global_normalized_ratiometric_data_list[num]]
                kdm = [x*100 for x in global_kdm_list[num]]
                #-- Clear the Plots --#
                ratio_axs[0].clear()
                ratio_axs[1].clear()
                #-- Redraw the titles --#
                ratio_axs[0].set_title("Normalized Ratio")
                ratio_axs[1].set_title("KDM")
                ratio_axs[0].set_ylabel("% Signal", fontweight="bold")
                ratio_axs[1].set_ylabel("% Signal", fontweight="bold")
                #-- Plot the Data --#
                ratio_axs[0].plot(x_axis, supposedly_global_norm, "ro", markersize=1)
                # normalized ratio of high and low freq'
                ratio_axs[1].plot(x_axis, kdm, "ro", markersize=1)
                ratio_fig.canvas.draw_idle()

    def post_analysis_kdm(self) -> None:
        """Adjust the High and Low frequencies used for KDM and ratiometric analysis."""
        global global_high_frequency,\
            global_low_frequency_offset,\
            global_low_frequency_slope,\
            global_low_frequency,\
            global_frequency_warning_label_exists,\
            global_wrong_frequency_label,\
            global_ratiometric_check
        global_high_frequency = int(self.high_frequency_entry.get())
        global_low_frequency = int(self.low_frequency_entry.get())
        global_low_frequency_offset = float(self.low_frequency_offset.get())
        global_low_frequency_slope = float(self.low_frequency_slope.get())
        #--- Reset the variable for the Warning Label (WrongFrequencyLabel) ---#
        if global_high_frequency not in global_frequency_list\
            and global_low_frequency in global_frequency_list:
            if global_frequency_warning_label_exists:
                global_wrong_frequency_label.grid_forget()
            global_wrong_frequency_label = ttk.Label(self.frequency_frame,\
                                                    text="High Frequency Does Not Exist",\
                                                        foreground="red")
            global_wrong_frequency_label.grid(row=6, column=0, columnspan=4)
            if not global_frequency_warning_label_exists:
                global_frequency_warning_label_exists = True
        elif global_high_frequency in global_frequency_list\
            and global_low_frequency not in global_frequency_list:
            if global_frequency_warning_label_exists:
                global_wrong_frequency_label.grid_forget()
            global_wrong_frequency_label = ttk.Label(self.frequency_frame,\
                                            text="Low Frequency Does Not Exist",\
                                                foreground="red")
            global_wrong_frequency_label.grid(row=6, column=0, columnspan=4)
            if not global_frequency_warning_label_exists:
                global_frequency_warning_label_exists = True
        elif global_high_frequency not in global_frequency_list\
            and global_low_frequency not in global_frequency_list:
            if global_frequency_warning_label_exists:
                global_wrong_frequency_label.grid_forget()
            global_wrong_frequency_label = ttk.Label(self.frequency_frame,\
                                            text="High and Low Frequencies Do Not Exist",\
                                                foreground="red")
            global_wrong_frequency_label.grid(row=6, column=0, columnspan=4)
            if not global_frequency_warning_label_exists:
                global_frequency_warning_label_exists = True
        else:
            global_high_low_dictionary[HighLow.HIGH] = global_high_frequency
            global_high_low_dictionary[HighLow.LOW] = global_low_frequency
            global_data_normalization.reset_ratiometric_data()
            #--- if a warning label exists, forget it ---#
            if global_frequency_warning_label_exists:
                global_wrong_frequency_label.grid_forget()
            #--- Tells RawVoltammogramVisualization to revisualize data
            # for new High and Low frequencies ---#
            global_ratiometric_check = True
            self._adjust_data()

    def post_analysis_normalization(self) -> None:
        """Adjust the normalization point for real-time normalization."""
        global global_normalization_point
        global_normalization_point = int(self.set_point_norm.get())
        file = int(global_file_label.cget("text"))
        if file >= global_normalization_point:
            global_wait_time.normalization_wait_time()
        else:
            global_norm_warning.config(foreground="red")
            global_norm_warning.config(text=f"File {global_normalization_point}" +\
                                       f"has not been analyzed yet")
        if global_analysis_complete:
            global_post_analysis._adjust_data()

    def data_export_top_level(self) -> None:
        """Create a TopLevel window for data export."""
        self.win = tk.Toplevel()
        self.win.wm_title("Post Analysis Data Export")
        self.export_top_level_exists = True
        #--- File Path ---#
        self.select_file_path = ttk.Button(self.win, style="On.TButton",\
                                           text=global_data_directory,\
                                            command=lambda: self.find_file(self.parent))
        self.select_file_path.grid(row=0, column=0, columnspan=2)
        self.no_selected_path = ttk.Label(self.win, text="No File Path Selected",\
                                        font=MEDIUM_FONT, foreground="red")
        self.path_warning_exists = False
        #--- File Handle Input ---#
        ttk.Label(self.win, text="Exported File Handle:", font=LARGE_FONT).\
            grid(row=4, column=0, columnspan=2)
        self.filehandle = ttk.Entry(self.win)
        self.filehandle.insert(tk.END, global_file_handle)
        self.filehandle.grid(row=5, column=0, columnspan=2, pady=5)
        self.electrode_label = ttk.Label(self.win, text="Select Electrodes:", font=LARGE_FONT)
        self.electrode_label.grid(row=10, column=0, sticky="nswe")
        self.electrode_count = tk.Listbox(self.win, relief="groove",\
                                          exportselection=0, width=10,\
                                          font=LARGE_FONT, height=6,\
                                            selectmode=tk.MULTIPLE, bd=3)
        self.electrode_count.bind("<<ListboxSelect>>", self.electrode_cur_select)
        self.electrode_count.grid(row=11, column=0, padx=10, sticky="nswe")
        for electrode in global_electrode_list:
            self.electrode_count.insert(tk.END, electrode)
        self.frequency_label = ttk.Label(self.win, text="Select Frequencies", font=LARGE_FONT)
        self.frequency_label.grid(row=10, column=1, padx=10)
        self.frequency_list = tk.Listbox(self.win, relief="groove",\
                                         exportselection=0, width=5,\
                                         font=LARGE_FONT, height=5,\
                                            selectmode=tk.MULTIPLE, bd=3)
        self.frequency_list.bind("<<ListboxSelect>>", self.frequency_cur_select)
        self.frequency_list.grid(row=11, column=1, padx=10, sticky="nswe")
        for frequency in global_frequency_list:
            self.frequency_list.insert(tk.END, frequency)
        ttk.Button(self.win, text="Export Data", command=self.post_analysis_data_export).\
            grid(row=15, column=0, columnspan=2)
        ttk.Button(self.win, text="Close", command=self.win.destroy).\
            grid(row=16, column=0, columnspan=2, pady=10)

    def electrode_cur_select(self, event) -> None:
        """Select the electrodes to be exported."""
        if global_electrode_count == 0:
            self.electrode_list_exists = False
            self.electrode_label.config(foreground="red")
        else:
            self.electrode_list_exists = True
            self.electrode_label.config(foreground="black")

    #--- Frequency Selection ---#
    def frequency_cur_select(self, event) -> None:
        """Select the frequencies to be exported."""
        if len(global_frequency_list) == 0:
            self.frequency_list_exists = False
            self.frequency_label.config(foreground="red")
        else:
            self.frequency_list_exists = True
            self.frequency_label.config(foreground="black")
            for var in global_frequency_list:
                var = int(var)

    def find_file(self, parent) -> None:
        """Let the user set the file path for data export."""
        global global_file_path, global_export_path, global_found_file_path
        try:
            ### prompt the user to select a  ###
            ### directory for  data analysis ###
            global_file_path = filedialog.askdirectory(parent=parent)
            global_file_path = global_file_path + "/"
            ### Path for directory in which the    ###
            ### exported .txt file will be placed  ###
            local_export_path: List[str] = global_file_path.split("/")
            #-- change the text of the find file button to the folder the user chose --#
            local_data_folder = f"{local_export_path[-3]}/{local_export_path[-2]}"
            self.select_file_path.config(style="On.TButton")
            self.select_file_path.config(text=local_data_folder)
            del local_export_path[-1]
            global_export_path = "/".join(local_export_path) + "/"
            ## Indicates that the user has selected a File Path ###
            global_found_file_path = True
            if self.path_warning_exists:
                self.no_selected_path.config(text="")
                self.no_selected_path.grid_forget()
        except FileNotFoundError:
            global_found_file_path = False
            self.no_selected_path.grid(row=1, column=0, columnspan=4)
            self.path_warning_exists = True
            self.select_file_path.config(style="Off.TButton")
            self.select_file_path.config(text="")

    def post_analysis_data_export(self) -> None:
        """Export the data after post-analysis."""
        global global_export_file_path
        local_file_handle = str(self.filehandle.get())
        global_export_file_path = global_export_path + local_file_handle
        post_analysis_export = TextFileExport().initialize()
        post_analysis_export.txt_file_normalization()

    def _reset(self) -> None:
        """Reset the experiment."""
        global global_high_already_reset,\
            global_low_already_reset,\
            global_already_reset,\
            global_analysis_already_initiated,\
            global_poison_pill,\
            global_key,\
            global_analysis_complete
        self.completion_value = 0
        global_analysis_complete = False
        if self.export_top_level_exists is True:
            self.win.destroy()
        global_key = 0
        global_poison_pill = True
        global_analysis_already_initiated = False # reset the start variable
        global_already_reset = True
        # Raise the initial user input frame
        self.show_frame(FRAME_INPUT)
        self.close_frame(FRAME_POST_ANALYSIS)
        ## Take resize weight away from the Visualization Canvas
        global_container.columnconfigure(1, weight=0)

    #--- Function to switch between visualization frames ---#
    def show_plot(self, frame: ttk.Frame) -> None:
        """Raise the given frame."""
        frame.tkraise()

    def show_frame(self, cont: str) -> None:
        """Raise the appropriate frame."""
        frame = global_show_frames[cont]
        frame.tkraise()

    def close_frame(self, cont: str) -> None:
        """Destroy the rames on reset."""
        frame = global_show_frames[cont]
        frame.grid_forget()
        # close all matplotlib figures
        plt.close("all")
        # destroy the frames holding the figures
        for frame in global_plot_values:
            frame.destroy()
        # destroy the container holding those frames
        global_plot_container.destroy()
#-----------------------------------------------------------------------------------------------#
                ###############################################################################
                ###############################################################################
                ###### Classes and Functions for Real-Time Tracking and Text File Export ######
                ###############################################################################
                ###############################################################################
class WaitTime():
    """Class for normalization signalling."""
    def __init__(self):
        global global_normalization_waiting
        global_normalization_waiting = False

    def normalization_wait_time(self) -> None:
        global global_normalization_waiting
        global_normalization_waiting = True

    def normalization_proceed(self) -> None:
        global global_normalization_waiting
        global_normalization_waiting = False

class Track():
    """Class for tracking the number of files analyzed."""
    def __init__(self):
        self.track_list: List[int] = [1]*global_number_of_files_to_process

    def tracking(self, file: int, frequency: Optional[int]) -> None:
        global global_ratiometric_check
        index: int = file - 1
        if self.track_list[index] == global_electrode_count:
            match global_analysis_method:
                case AnalysisMethod.CONTINUOUS_SCAN:
                    _update_global_lists(file)
                    global_data_normalization.renormalize_data(file)
                    if global_text_file_export_activated:
                        global_text_file_export.continuous_scan_export(file)
                    #--- if the high and low frequencies have been changed, adjust the data ---#
                    if global_ratiometric_check:
                        global_data_normalization.reset_ratiometric_data()
                        #-- if the data is being exported, reset the exported data file --#
                        if global_text_file_export_activated:
                            global_text_file_export.txt_file_normalization()
                        global_ratiometric_check = False
                case AnalysisMethod.FREQUENCY_MAP:
                    if self.track_list[index] == global_electrode_count:
                        if global_text_file_export_activated:
                            if frequency is None:
                                internal_error("tracking: frequency is None")
                            else:
                                global_text_file_export.frequency_map_export(file, frequency)
            self.track_list[index] = 1
        else:
            self.track_list[index] += 1
#-------------------------------------------------------------------------------------------------#
class TextFileExport():
    """Class for exporting data to a .txt file."""
    def __init__(self,\
                 electrodes: Optional[List[int]] = None,\
                 frequencies: Optional[List[int]] = None):
        #declarations for fields that will be initialized in other methods
        self.electrode_list: List[int]
        self.electrode_count: int
        self.frequency_list: List[int]
        self.text_file_handle: str

    def initialize(self,\
                 electrodes: Optional[List[int]] = None,\
                 frequencies: Optional[List[int]] = None):
        if electrodes is None:
            self.electrode_list = global_electrode_list
        else:
            self.electrode_list = electrodes
        self.electrode_count = len(self.electrode_list)
        if frequencies is None:
            self.frequency_list = global_frequency_list
        else:
            self.frequency_list = frequencies
        self.text_file_handle = global_export_file_path
        match global_analysis_method:
            case AnalysisMethod.CONTINUOUS_SCAN:
                txt_list: List[str] = []
                match global_plot_summary_mode:
                    case PlotSummaryMode.PHE:
                        match global_peak_method:
                            case PeakMethod.POLY:
                                txt_list.append("Peak Method: Poly Fit")
                            case PeakMethod.GAUSS:
                                txt_list.append("Peak Method: Gauss")
                        with open(self.text_file_handle, "w+", encoding="utf-8", newline="") as output:
                            writer = csv.writer(output, delimiter=" ")
                            writer.writerow(txt_list)
                txt_list: List[str] = []
                txt_list.append("File")
                txt_list.append("Time(Hrs)")
                for frequency in self.frequency_list:
                    for electrode in self.electrode_list:
                        match global_plot_summary_mode:
                            case PlotSummaryMode.PHE:
                                txt_list.append("PeakHeight_E%d_%dHz" % (electrode, frequency))
                                match global_peak_method:
                                    case PeakMethod.GAUSS:
                                        txt_list.append("PeakLocation_E%d_%dHz" % (electrode, frequency))
                            case PlotSummaryMode.AUC:
                                txt_list.append("AUC_E%d_%dHz" % (electrode, frequency))
                if self.electrode_count > 1:
                    for frequency in self.frequency_list:
                        match global_plot_summary_mode:
                            case PlotSummaryMode.PHE:
                                txt_list.append("Avg_PeakHeight_%dHz" % frequency)
                            case PlotSummaryMode.AUC:
                                txt_list.append("Avg_AUC_%dHz" % frequency)
                for frequency in self.frequency_list:
                    for electrode in self.electrode_list:
                        txt_list.append(f"Norm_E{electrode}_{frequency}Hz")
                if self.electrode_count > 1:
                    for frequency in self.frequency_list:
                        txt_list.append(f"Average_Norm_{frequency}Hz")
                    for frequency in self.frequency_list:
                        txt_list.append(f"SD_Norm_{frequency}Hz")
                if len(self.frequency_list) > 1:
                    for electrode in self.electrode_list:
                        txt_list.append(f"NormalizedRatio_E{electrode}")
                    if self.electrode_count > 1:
                        txt_list.append("NormalizedRatioAvg")
                        txt_list.append("NormalizedRatioSTD")
                    for electrode in self.electrode_list:
                        txt_list.append(f"KDM_E{electrode}")
                    if self.electrode_count > 1:
                        txt_list.append("AvgKDM")
                        txt_list.append("KDM_STD")
                with open(self.text_file_handle, "a", encoding="utf-8", newline="") as output:
                    writer = csv.writer(output, delimiter=" ")
                    writer.writerow(txt_list)
            case AnalysisMethod.FREQUENCY_MAP:
                txt_list = []
                txt_list.append("Frequency(Hz)")
                e_count = 1
                for _ in global_frame_list:
                    match global_plot_summary_mode:
                        case PlotSummaryMode.PHE:
                            txt_list.append(f"PeakHeight_E{e_count}(µA)")
                        case PlotSummaryMode.AUC:
                            txt_list.append(f"AUC_E{e_count}")
                    txt_list.append(f"Charge_E{e_count}(µC)")
                    e_count += 1
                if self.electrode_count > 1:
                    txt_list.append("Avg.PeakHeight(µA)")
                    txt_list.append("Standard_Deviation(µA)")
                    txt_list.append("Avg.Charge(µC)")
                    txt_list.append("Standard_Deviation(µC)")
                with open(self.text_file_handle, "w+", encoding="utf-8", newline="") as output:
                    writer = csv.writer(output, delimiter=" ")
                    writer.writerow(txt_list)
        return self

    def continuous_scan_export(self, file: int) -> None:
        """Export the data from the current file."""
        normalized_frequency_currents: List[float]
        norm_list: List[float]
        kdm_list: List[float]
        running_sum: float
        average: float
        average_norm: float
        index: int = file - 1
        output_list: List[str] = []
        output_list.append(str(file))
        output_list.append(str(global_sample_list[index]))
        #--- Peak Height ---#
        for count in range(len(global_frequency_list)):
            for num in range(global_electrode_count):
                output_list.append(str(global_data_list[num][count][index]))
                match global_peak_method:
                    case PeakMethod.GAUSS:
                        output_list.append(str(global_peak_list[num][count][index]))
        #--- Avg. Peak Height ---#
        if self.electrode_count > 1:
            for count in range(len(global_frequency_list)):
                running_sum = 0
                for num in range(global_electrode_count):
                    running_sum += global_data_list[num][count][index]
                average = running_sum/global_electrode_count
                output_list.append(str(average))
        #--- Peak Height/AUC Data Normalization ---#
        for count in range(len(global_frequency_list)):
            for num in range(global_electrode_count):
                if global_frequency_list[count] == global_high_low_dictionary[HighLow.LOW]:
                    output_list.append(str(global_offset_normalized_data_list[num][index]))
                else:
                    output_list.append(str(global_normalized_data_list[num][count][index]))
        #--- Average normalized data across all electrodes for each frequency ---#
        if self.electrode_count > 1:
            for count in range(len(global_frequency_list)):
                normalized_frequency_currents = []
                for num in range(global_electrode_count):
                    if global_frequency_list[count] == global_high_low_dictionary[HighLow.LOW]:
                        normalized_frequency_currents.\
                            append(global_offset_normalized_data_list[num][index])
                    else:
                        normalized_frequency_currents.\
                            append(global_normalized_data_list[num][count][index])
                average_norm = sum(normalized_frequency_currents)/global_electrode_count
                output_list.append(str(average_norm))
        #--- Standard Deviation ---#
        if self.electrode_count > 1:
            for count in range(len(global_frequency_list)):
                normalized_frequency_currents = []
                for num in range(global_electrode_count):
                    normalized_frequency_currents.\
                        append(global_normalized_data_list[num][count][index])
                average_norm = sum(normalized_frequency_currents)/global_electrode_count
                std_list = [(x - average_norm)**2 for x in normalized_frequency_currents]
                standard_deviation = float(sqrt(sum(std_list)/(global_electrode_count - 1)))
                output_list.append(str(standard_deviation))
        if len(global_frequency_list) > 1:
            #--- Append Normalized Ratiometric Data ---#
            norm_list = []
            for num in range(global_electrode_count):
                output_list.append(str(global_normalized_ratiometric_data_list[num][index]))
                norm_list.append(global_normalized_ratiometric_data_list[num][index])
            if self.electrode_count > 1:
                norm_average: float = sum(norm_list)/global_electrode_count
                output_list.append(str(norm_average))
                norm_std_list: List[float] = [(x - norm_average)**2 for x in norm_list]
                norm_standard_deviation = sqrt(sum(norm_std_list)/(global_electrode_count - 1))
                output_list.append(str(norm_standard_deviation))
            #--- Append KDM ---#
            kdm_list = []
            for num in range(global_electrode_count):
                output_list.append(str(global_kdm_list[num][index]))
                kdm_list.append(global_kdm_list[num][index])
            if self.electrode_count > 1:
                kdm_average: float = sum(kdm_list)/global_electrode_count
                output_list.append(str(kdm_average))
                kdm_std_list: List[float] = [(x - kdm_average)**2 for x in kdm_list]
                kdm_std: float = sqrt(sum(kdm_std_list)/(global_electrode_count - 1))
                output_list.append(str(kdm_std))
        #--- Write the data into the .txt file ---#
        with open(self.text_file_handle, "a", encoding="utf-8", newline="") as text_io_wrapper:
            writer = csv.writer(text_io_wrapper, delimiter=" ")
            writer.writerow(output_list)
        with open(self.text_file_handle, "r", encoding="utf-8", newline="") as filecontents:
            filedata = filecontents.read()
        filedata = filedata.replace("[", "")
        filedata = filedata.replace("\"", "")
        filedata = filedata.replace("]", "")
        filedata = filedata.replace(",", "")
        filedata = filedata.replace("'", "")
        with open(self.text_file_handle, "w", encoding="utf-8", newline="") as output:
            output.write(filedata)

    def frequency_map_export(self, file: int, frequency: int) -> None:
        """Export the data from the current file."""
        output_list: List[str] = []
        index: int = file - 1
        running_sum: float
        average: float
        try:
            output_list.append(str(frequency))
            count = global_frequency_dict[frequency]
            # Peak Height / AUC
            for num in range(global_electrode_count):
                output_list.append(str(global_data_list[num][count][index]))
                output_list.append(str(global_data_list[num][count][index]/frequency))
            # Average Peak Height / AUC
            if self.electrode_count > 1:
                running_sum = 0
                for num in range(global_electrode_count):
                    running_sum += global_data_list[num][count][index]
                average = running_sum/global_electrode_count
                output_list.append(str(average))
                # Standard Deviation of a Sample across all electrodes
                # for Peak Height/AUC
                std_list: List[float] = []
                for num in range(global_electrode_count):
                    std_list.append(global_data_list[num][count][index])
                std_list = [(value - average)**2 for value in std_list]
                standard_deviation = sqrt(sum(std_list)/(global_electrode_count - 1))
                output_list.append(str(standard_deviation))
                #-- Average Charge --#
                avg_charge = average/frequency
                output_list.append(str(avg_charge))
                #-- Charge STD --#
                std_list = []
                for num in range(global_electrode_count):
                    std_list.append(global_data_list[num][count][index])
                std_list = [x/frequency for x in std_list]
                std_list = [(value - avg_charge)**2 for value in std_list]
                charge_standard_deviation = sqrt(sum(std_list)/(global_electrode_count - 1))
                output_list.append(str(charge_standard_deviation))
            #--- Write the data into the .txt file ---#
            with open(self.text_file_handle, "a", encoding="utf-8", newline="") as output:
                writer = csv.writer(output, delimiter=" ")
                writer.writerow(output_list)
            with open(self.text_file_handle, "r", encoding="utf-8", newline="") as filecontents:
                filedata = filecontents.read()
            filedata = filedata.replace("[", "")
            filedata = filedata.replace("\"", "")
            filedata = filedata.replace("]", "")
            filedata = filedata.replace(",", "")
            filedata = filedata.replace("'", "")
            with open(self.text_file_handle, "w", encoding="utf-8", newline="") as output:
                output.write(filedata)
        except Exception as exception:
            print("Error in frequency_map_export", exception)
            time.sleep(3)

    def txt_file_normalization(self) -> None:
        """Normalize the data within the text file."""
        try:
            #--- reinitialize the .txt file ---#
            self.initialize(electrodes=self.electrode_list, frequencies=self.frequency_list)
            #--- rewrite the data for the files that have already been analyzed and normalize
            # them to the new standard---#
            if global_analysis_complete:
                analysis_range = len(global_file_list)
            else:
                analysis_range = len(global_file_list) - 1
            for index in range(analysis_range):
                file: int = index + 1
                txt_list: List[str] = []
                #AvgList = [] #not used
                txt_list.append(str(file))
                txt_list.append(str(global_sample_list[index]))
                #--- peak height ---#
                for frequency in self.frequency_list:
                    count = global_frequency_dict[frequency]
                    for electrode in self.electrode_list:
                        num = global_electrode_dict[electrode]
                        txt_list.append(str(global_data_list[num][count][index]))
                        match global_peak_method:
                            case PeakMethod.GAUSS:
                                txt_list.append(str(global_peak_list[num][count][index]))
                #--- Avg. Peak Height ---#
                if self.electrode_count > 1:
                    for frequency in self.frequency_list:
                        count = global_frequency_dict[frequency]
                        running_sum: float = 0
                        for electrode in self.electrode_list:
                            num = global_electrode_dict[electrode]
                            running_sum += global_data_list[num][count][index]
                        average = running_sum/self.electrode_count
                        txt_list.append(str(average))
                #--- Data Normalization ---#
                for frequency in self.frequency_list:
                    count = global_frequency_dict[frequency]
                    normalized_frequency_currents: List[float] = []
                    for electrode in self.electrode_list:
                        num = global_electrode_dict[electrode]
                        if frequency == global_high_low_dictionary[HighLow.LOW]:
                            txt_list.append(str(global_offset_normalized_data_list[num][index]))
                        else:
                            txt_list.append(str(global_normalized_data_list[num][count][index]))
                #--- Average Data Normalization ---#
                if self.electrode_count > 1:
                    for frequency in self.frequency_list:
                        count = global_frequency_dict[frequency]
                        normalized_frequency_currents = []
                        for electrode in self.electrode_list:
                            num = global_electrode_dict[electrode]
                            if frequency == global_high_low_dictionary[HighLow.LOW]:
                                normalized_frequency_currents.\
                                    append(global_offset_normalized_data_list[num][index])
                            else:
                                normalized_frequency_currents.\
                                    append(global_normalized_data_list[num][count][index])
                        average_norm = sum(normalized_frequency_currents)/self.electrode_count
                        txt_list.append(str(average_norm))
                    #--- Standard Deviation between electrodes ---#
                    for frequency in self.frequency_list:
                        count = global_frequency_dict[frequency]
                        normalized_frequency_currents = []
                        for electrode in self.electrode_list:
                            num = global_electrode_dict[electrode]
                            if frequency == global_high_low_dictionary[HighLow.LOW]:
                                normalized_frequency_currents.\
                                    append(global_offset_normalized_data_list[num][index])
                            else:
                                normalized_frequency_currents.\
                                    append(global_normalized_data_list[num][count][index])
                        average_norm = sum(normalized_frequency_currents)/self.electrode_count
                        std_list = [(x - average_norm)**2 for x in normalized_frequency_currents]
                        standard_deviation = sqrt(sum(std_list)/(self.electrode_count - 1))
                        txt_list.append(str(standard_deviation))
                if len(self.frequency_list) > 1:
                    #--- Append Normalized Ratiometric Data ---#
                    norm_list: List[float] = []
                    for electrode in self.electrode_list:
                        num = global_electrode_dict[electrode]
                        txt_list.append(str(global_normalized_ratiometric_data_list[num][index]))
                        norm_list.append(global_normalized_ratiometric_data_list[num][index])
                    if self.electrode_count > 1:
                        norm_average = sum(norm_list)/self.electrode_count
                        txt_list.append(str(norm_average))
                        norm_std_list = [(x - norm_average)**2 for x in norm_list]
                        norm_standard_deviation = sqrt(sum(norm_std_list)/(self.electrode_count-1))
                        txt_list.append(str(norm_standard_deviation))
                    #--- Append KDM ---#
                    kdm_list = []
                    for electrode in self.electrode_list:
                        num = global_electrode_dict[electrode]
                        txt_list.append(str(global_kdm_list[num][index]))
                        kdm_list.append(global_kdm_list[num][index])
                    if self.electrode_count > 1:
                        kdm_average: float = sum(kdm_list)/self.electrode_count
                        txt_list.append(str(kdm_average))
                        kdm_std_list: List[float] = [(x - kdm_average)**2 for x in kdm_list]
                        kdm_std: float = sqrt(sum(kdm_std_list)/(self.electrode_count - 1))
                        txt_list.append(str(kdm_std))
                #--- Write the data into the .txt file ---#
                with open(self.text_file_handle, "a", encoding="utf-8", newline="") as output:
                    writer = csv.writer(output, delimiter=" ")
                    writer.writerow(txt_list)
                with open(self.text_file_handle, "r", encoding="utf-8", newline="") as filecontents:
                    filedata = filecontents.read()
                filedata = filedata.replace("[", "")
                filedata = filedata.replace("\"", "")
                filedata = filedata.replace("]", "")
                filedata = filedata.replace(",", "")
                filedata = filedata.replace("'", "")
                with open(self.text_file_handle, "w", encoding="utf-8", newline="") as output:
                    output.write(filedata)
        except Exception as exception:
            print("Error in txt_file_normalization", exception)
            time.sleep(0.1)

#--------------------------------------------------------------------------------------- #
class AnalysisMethod (Enum):
    CONTINUOUS_SCAN = 0
    FREQUENCY_MAP = 1
    def __str__(self) -> str:
        match self:
            case AnalysisMethod.CONTINUOUS_SCAN:
                return "Continuous Scan"
            case AnalysisMethod.FREQUENCY_MAP:
                return "Frequency Map"
class KDMMethod (Enum):
    OLD = 0
    NEW = 1
    def __str__(self) -> str:
        match self:
            case KDMMethod.OLD:
                return "Old KDM"
            case KDMMethod.NEW:
                return "New KDM"
class PeakMethod (Enum):
    POLY = 0
    GAUSS = 1
    def __str__(self) -> str:
        match self:
            case PeakMethod.POLY:
                return "Poly Fit"
            case PeakMethod.GAUSS:
                return "Gauss"
global_analysis_method: AnalysisMethod = AnalysisMethod.CONTINUOUS_SCAN
global_kdm_method: KDMMethod = KDMMethod.OLD
global_peak_method: PeakMethod = PeakMethod.POLY
global_container: ttk.Frame
global_number_of_files_to_process: int
global_injection_point: Optional[int]
global_initialized_normalization: bool
global_ratiometric_check: bool
global_normalization_point: int
global_text_file_export_activated: bool
global_injection_selected: bool
global_search_interval: int
global_high_frequency: int
global_low_frequency: int
class HighLow (Enum):
    HIGH = 0
    LOW = 1
global_high_low_dictionary: Dict[HighLow, int]
global_normalization_vault: List[int]
global_min_kdm: float
global_max_kdm: float
global_min_norm: float
global_max_norm: float
global_min_raw: float
global_max_raw: float
global_min_data: float
global_max_data: float
global_sample_rate: float
global_queue: Queue
global_export_path: str = ""
global_data_directory: str
global_file_handle: str
global_export_file_path: str
global_wait_time: WaitTime
global_track: Track
global_data_normalization: DataNormalization
global_post_analysis: PostAnalysis
global_file_label: ttk.Label
global_real_time_sample_label: ttk.Label
global_set_point_norm: ttk.Entry
global_norm_warning: ttk.Label
global_low_frequency_entry: ttk.Entry
global_show_frames: Dict[str, ttk.Frame]
global_plot_values: List[ttk.Frame]
global_plot_frames: Dict[str, ttk.Frame]
global_empty_plots: List[matplotlib.lines.Line2D]
global_x_axis_mode: PlotTimeReportingMode = PlotTimeReportingMode.EXPERIMENT_TIME
global_frame_reference: Union[ContinuousScanVisualizationFrame, FrequencyMapVisualizationFrame]
global_plot_container: ttk.Frame
global_kdm_list: List[List[float]]
global_normalization_waiting: bool
global_figures: List[Tuple[matplotlib.figure.Figure, List[List[matplotlib.axes.Axes]]]]
global_electrode_count: int
global_electrode_list: List[int]
global_electrode_dict: Dict[int, int]
global_frequency_list: List[int]
global_frequency_dict: Dict[int, int]
global_resize_interval: int
global_high_xstart_entry: ttk.Entry # delete once freq. map is fixed
global_low_xstart_entry: ttk.Entry # delete once freq. map is fixed
global_high_xend_entry: ttk.Entry # delete once freq. map is fixed
global_low_xend_entry: ttk.Entry # delete once freq. map is fixed
global_xstart_entry: List[List[ttk.Entry]]
global_xend_entry: List[List[ttk.Entry]]
global_high_frequency_entry: ttk.Entry
global_normalization_var: tk.StringVar
global_high_xstart: float # delete once freq. map is fixed
global_low_xstart: float # delete once freq. map is fixed
global_high_xend: float # delete once freq. map is fixed
global_low_xend: float # delete once freq. map is fixed
global_xstart: List[List[float]]
global_xend: List[List[float]]
global_wrong_frequency_label: ttk.Label
global_plot_list_continuous_scan: List[List[Tuple[matplotlib.lines.Line2D,\
                                                    matplotlib.lines.Line2D,\
                                                    matplotlib.lines.Line2D,\
                                                    matplotlib.lines.Line2D,\
                                                    matplotlib.lines.Line2D,\
                                                    matplotlib.lines.Line2D,\
                                                    Polygon,\
                                                    matplotlib.lines.Line2D]]]
global_plot_list_frequency_map: List[Tuple[matplotlib.lines.Line2D,\
                                           matplotlib.lines.Line2D,\
                                            matplotlib.lines.Line2D]]
global_frame_list: List[str]
global_ratiometric_plots: List[List[matplotlib.lines.Line2D]]
global_text_file_export: TextFileExport
global_file_list: List[int]
global_animations: List[ElectrochemicalAnimation]
global_data_list: List[List[List[float]]]
global_peak_list: List[List[List[float]]]
global_sample_list: List[float]
global_normalized_data_list: List[List[List[float]]]
global_offset_normalized_data_list: List[List[float]]
global_normalized_ratiometric_data_list: List[List[float]]
global_ratiometric_figures: List[Tuple[matplotlib.figure.Figure, List[matplotlib.axes.Axes]]]
                    ############################################
                    ### Initialize GUI to start the program  ###
                    ############################################
if __name__ == "__main__":
    root = tk.Tk()
    app = MainWindow(root)
    #root.title("Electrochemical Analysis")
    style = ttk.Style()
    style.configure("On.TButton", foreground="blue", font=MEDIUM_FONT, relief="raised", border=100)
    style.configure("Off.TButton", foreground="black", font=MEDIUM_FONT, relief="sunken", border=5)
    style.configure("TFrame", borderwidth=2)
    while True:
        #--- initiate the main loop ---#
        try:
            root.mainloop()
        #--- escape scrolling error ---#
        except UnicodeDecodeError:
            pass
