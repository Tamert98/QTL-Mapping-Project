# ================================ #
#  LoadingDataFrame GUI Component  #
#  For Genetic Map Preparation     #
# ================================ #

# --- External Libraries ---
from customtkinter import *
from PIL import Image
import os
import threading
import time

# --- Internal Modules: Data Handling ---
from vcf_data_handler import (
    parse_vcf_file,
    parse_trait_file,
    generate_combined_genetic_maps,
    generate_combined_genetic_maps_filtered,
    select_evenly_spaced_markers,
)

# --- Internal Modules: Plotting ---
from plot_utils import (
    compute_global_physical_max,
    generate_genetic_map_images_and_pdf,
    generate_comparative_genetic_map_images_and_pdf,
    print_heatmap_pdf_path,
)

# --- Internal Modules: Frame Navigation ---
from gui_frames.smt_sim_runframe import SMT_SIM_RunFrame


class LoadingDataFrame(CTkFrame):
    """
    Frame for loading and preprocessing genotype and trait data.
    Generates genetic maps and heatmaps, then transitions to the SMT+SIM stage.
    """

    def __init__(self, master, vcf_path, trait_path, on_back, styles, on_sim_done):
        """
        Initialize the LoadingDataFrame.

        Args:
            master (Tk): Parent window.
            vcf_path (str): Path to the input VCF file.
            trait_path (str): Path to the trait data file.
            on_back (func): Callback to go back to previous frame.
            styles (dict): Styling dictionary for buttons.
            on_sim_done (func): Callback when SIM stage finishes.
        """
        super().__init__(master, width=1100, height=650)
        self.master = master
        self.vcf_path = vcf_path
        self.trait_path = trait_path
        self.on_back = on_back
        self.on_sim_done = on_sim_done
        self.styles = styles

        # === Data containers ===
        self.vcf_data = None
        self.sample_names = None
        self.traits = None
        self.unfiltered_maps = None
        self.filtered_maps = None
        self.selected_markers = None

        # === State flags and widgets ===
        self.data_ready = False
        self.maps_saved_label = None

        self.pack(fill="both", expand=True)
        self.build_ui()

        # Start loading and generation in background
        threading.Thread(target=self.load_and_generate, daemon=True).start()

    def build_ui(self):
        """Construct the GUI layout."""
        width, height = 1100, 650
        self.master.geometry(f"{width}x{height}")
        self.master.update_idletasks()
        x = (self.master.winfo_screenwidth() // 2) - (width // 2)
        y = (self.master.winfo_screenheight() // 2) - (height // 2) - 100
        self.master.geometry(f"{width}x{height}+{x}+{y}")

        # === Step 1 Image (left panel) ===
        try:
            step1_image = CTkImage(light_image=Image.open("gui_frames/Step1.png"), size=(250, 650))
            CTkLabel(self, image=step1_image, text="").place(x=0, y=0)
            self.step1_image = step1_image
        except Exception as e:
            print(f"Step1.png not found: {e}")

        # === Right-side content panel ===
        content_frame = CTkFrame(self, fg_color="transparent")
        content_frame.place(relx=0.62, rely=0.5, anchor="center")

        self.title_label = CTkLabel(content_frame, text="Loading Data", font=("Segoe UI", 24, "bold"))
        self.title_label.pack(pady=(10, 10))

        self.intro_label = CTkLabel(
            content_frame,
            text="This process will prepare genetic maps and heatmaps used in the QTL analysis pipeline.",
            font=("Segoe UI", 16, "bold"),
            justify="left",
            wraplength=650,
            text_color="#ffffff"
        )
        self.intro_label.pack(pady=5)

        # === Status Labels ===
        self.vcf_label = CTkLabel(content_frame, text="Loading VCF data...", font=("Segoe UI", 20, "bold"), text_color="#ffffff")
        self.vcf_label.pack(pady=(20, 5), anchor="w")

        self.trait_label = CTkLabel(content_frame, text="Waiting to load trait data...", font=("Segoe UI", 20, "bold"), text_color="#ffffff")
        self.trait_label.pack(pady=(10, 5), anchor="w")

        self.maps_saved_label = CTkLabel(content_frame, text="", font=("Segoe UI", 16, "bold"), text_color="#008000")
        self.maps_saved_label.pack(pady=(10, 10), anchor="w")

        self.output_label = CTkLabel(content_frame, text="", font=("Segoe UI", 14, "bold"), justify="left", text_color="#008000")
        self.output_label.pack(pady=(5, 20), anchor="w")

        # === Navigation Buttons ===
        nav_frame = CTkFrame(content_frame, fg_color="transparent")
        nav_frame.pack(pady=10)

        self.back_btn = CTkButton(nav_frame, text="Back", font=("Segoe UI", 14, "bold"),
                                  width=120, command=self.on_back, **self.styles["white"])
        self.back_btn.pack(side="left", padx=(0, 20))

    def load_and_generate(self):
        """Load VCF and trait data, generate maps, and continue to SMT+SIM."""
        # --- Step 1: Parse input files ---
        self.vcf_data, self.sample_names = parse_vcf_file(self.vcf_path)
        self.vcf_label.configure(text="✔ VCF data has been loaded", text_color="#008000")
        time.sleep(0.5)

        _, self.traits = parse_trait_file(self.trait_path)
        self.trait_label.configure(text="✔ Trait data has been loaded", text_color="#008000")
        time.sleep(1)

        # --- Step 2: Start analysis and update UI ---
        self.vcf_label.configure(text="Generating genetic maps and heatmaps...", text_color="#ffffff")
        self.trait_label.configure(text="")
        self.title_label.configure(text="Analysing Genotypes")

        try:
            step2_image = CTkImage(light_image=Image.open("gui_frames/Step2.png"), size=(250, 650))
            CTkLabel(self, image=step2_image, text="").place(x=0, y=0)
            self.step2_image = step2_image
        except Exception as e:
            print(f"Step2.png not found: {e}")
        time.sleep(1.2)

        # --- Step 3: Genetic map and heatmap generation ---
        self.unfiltered_maps, comp_unfiltered = generate_combined_genetic_maps(self.vcf_data, self.sample_names)
        self.filtered_maps, comp_filtered = generate_combined_genetic_maps_filtered(self.vcf_data, self.sample_names)

        global_xmax = compute_global_physical_max(self.unfiltered_maps)
        generate_genetic_map_images_and_pdf(self.filtered_maps, global_xmax)
        generate_comparative_genetic_map_images_and_pdf(comp_unfiltered, comp_filtered, global_xmax)
        print_heatmap_pdf_path(output_dir="Results/Genetic_Maps")

        # --- Step 4: Marker selection ---
        self.selected_markers = select_evenly_spaced_markers(self.filtered_maps, self.vcf_data)

        # --- Step 5: Final status update ---
        self.vcf_label.configure(text="✔ Genetic maps and heatmaps have been generated", text_color="#008000")
        self.maps_saved_label.configure(text="✔ Genetic maps, heatmaps, and comparison maps saved to Results folder.", text_color="#008000")
        self.output_label.configure(text=(
            "Genetic map PDF: Results/Genetic_Maps/final_genetic_maps.pdf\n"
            "Unfiltered comparison: Results/CompareOfDistances_unfiltered/comparison_genetic_maps.pdf\n"
            "Filtered comparison: Results/CompareOfDistances_filtered/comparison_filtered_genetic_maps.pdf\n"
            "Heatmap PDF: Results/Genetic_Maps/final_heatmaps.pdf"
        ), text_color="#ffffff")
        self.title_label.configure(text="Analysed Genotypes ✔", text_color="#008000")

        self.data_ready = True
        self.after(5000, self.move_to_smt_sim)

    def move_to_smt_sim(self):
        """Destroy this frame and move to SMT+SIM analysis stage."""
        self.destroy()
        SMT_SIM_RunFrame(
            master=self.master,
            vcf_data=self.vcf_data,
            traits=self.traits,
            sample_names=self.sample_names,
            selected_markers=self.selected_markers,
            styles=self.styles,
            on_done=self.on_sim_done,
            genetic_maps_unfiltered=self.unfiltered_maps
        )
