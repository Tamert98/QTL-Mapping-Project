from customtkinter import *
from PIL import Image
import os

from gui_frames.file_selection_frame import FileSelectionFrame
from gui_frames.loading_data_frame import LoadingDataFrame
from gui_frames.smt_results_frame import SMTResultsFrame
from gui_frames.file_format_window import open_format_window
from styles import get_styles

# === Theme & Styles ===
current_mode = "dark"
red_button_style = {}
white_button_style = {}
bg_image = None

def apply_styles_for_mode():
    global red_button_style, white_button_style
    styles = get_styles(current_mode)
    red_button_style = styles["red"]
    white_button_style = styles["white"]

def toggle_theme():
    global current_mode
    current_mode = "light" if current_mode == "dark" else "dark"
    set_appearance_mode(current_mode)
    clear_main_widgets()
    setup_main_page()

# === App Setup ===
set_appearance_mode("dark")
set_default_color_theme("blue")

app = CTk()
app.title("QTL Mapping Project")
width, height = 1100, 650
app.update_idletasks()
screen_width = app.winfo_screenwidth()
screen_height = app.winfo_screenheight()
x = (screen_width // 2) - (width // 2)
y = (screen_height // 2) - (height // 2) - 100
app.geometry(f"{width}x{height}+{x}+{y}")

try:
    app.iconbitmap("AppIcon.ico")
except:
    print("Icon not found — skipping.")

def clear_main_widgets():
    for widget in app.winfo_children():
        widget.destroy()

def launch_file_selection():
    clear_main_widgets()
    FileSelectionFrame(
        master=app,
        on_back=setup_main_page,
        on_next=go_to_next_step,
        styles={"red": red_button_style, "white": white_button_style}
    )

def go_to_next_step(vcf_path, trait_path):
    def after_loading_done(vcf_data, sample_names, traits, unfiltered_maps, filtered_maps):
        app.vcf_data = vcf_data
        app.sample_names = sample_names
        app.traits = traits
        app.unfiltered_maps = unfiltered_maps
        app.filtered_maps = filtered_maps
        launch_smt_results_frame()

    clear_main_widgets()
    LoadingDataFrame(
        master=app,
        vcf_path=vcf_path,
        trait_path=trait_path,
        on_done=after_loading_done,
        styles={"red": red_button_style, "white": white_button_style}
    )

def launch_smt_results_frame():
    clear_main_widgets()
    SMTResultsFrame(
        master=app,
        vcf_data=app.vcf_data,
        traits=app.traits,
        sample_names=app.sample_names,
        genetic_maps_unfiltered=app.unfiltered_maps,
        styles={"red": red_button_style, "white": white_button_style}
    )

def setup_main_page():
    global bg_image, title_label, text_box, theme_toggle_button

    width, height = 1100, 650
    app.update_idletasks()
    screen_width = app.winfo_screenwidth()
    screen_height = app.winfo_screenheight()
    x = (screen_width // 2) - (width // 2)
    y = (screen_height // 2) - (height // 2) - 100
    app.geometry(f"{width}x{height}+{x}+{y}")
    clear_main_widgets()

    try:
        bg_image = CTkImage(light_image=Image.open("sideimg.png"), size=(250, 650))
        CTkLabel(app, image=bg_image, text="").place(x=0, y=0)
    except:
        print("Background image not found.")

    theme_toggle_button = CTkButton(
        app,
        text="Apply Light Mode" if current_mode == "dark" else "Apply Dark Mode",
        command=toggle_theme,
        font=("Segoe UI", 12, "bold"),
        width=140,
        fg_color="#222222" if current_mode == "dark" else "#dddddd",
        hover_color="#333333" if current_mode == "dark" else "#cccccc",
        text_color="#ffffff" if current_mode == "dark" else "#000000"
    )
    theme_toggle_button.place(relx=1.0, x=-20, y=20, anchor="ne")

    title_label = CTkLabel(
        app,
        text="Welcome to the QTL Mapping Project",
        font=("Segoe UI", 22, "bold"),
        text_color="#ffffff" if current_mode == "dark" else "#000000"
    )
    title_label.pack(pady=(20, 10), padx=(270, 20))

    welcome_message = (
        "This project was developed by Tamer Talhami and Eyas Rizik.\n\n"
        "It focuses on identifying the DNA regions responsible for specific trait values "
        "using two powerful QTL mapping algorithms:\n"
        "• SMT (Single Marker Test): a simpler, fast approach.\n"
        "• SIM (Single Interval Mapping): a more advanced method.\n\n"
        "To begin:\n"
        "1. Prepare your VCF and trait files.\n"
        "2. Ensure the formats are correct.\n"
        "3. Click 'File Format Requirements' to learn more.\n"
        "4. Then click 'Next' to proceed."
    )

    text_box = CTkTextbox(
        app,
        width=750,
        height=350,
        font=("Segoe UI", 14),
        wrap="word",
        border_width=1,
        border_color="#333333"
    )
    text_box.insert("0.0", welcome_message)
    text_box.configure(state="disabled")
    text_box.pack(pady=10, padx=(270, 20))

    button_frame = CTkFrame(app, fg_color="transparent")
    button_frame.place(relx=0.72, rely=0.75, anchor="n")

    format_button = CTkButton(
        button_frame,
        text="File Format Requirements",
        font=("Segoe UI", 14, "bold"),
        width=240,
        height=50,
        command=lambda: open_format_window({"red": red_button_style, "white": white_button_style}, app),
        **white_button_style
    )

    next_button = CTkButton(
        button_frame,
        text="Next",
        font=("Segoe UI", 14, "bold"),
        width=200,
        height=50,
        command=launch_file_selection,
        **red_button_style
    )

    format_button.grid(row=0, column=0, padx=(0, 150))
    next_button.grid(row=0, column=1, padx=(0, 220))

# === START ===
apply_styles_for_mode()
setup_main_page()
app.mainloop()
