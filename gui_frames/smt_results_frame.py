# ======================================================= #
#  SMT Results Viewer Frames                              #
#  - SMTConcatenatedResultsFrame: Genome-wide view        #
#  - SMTPerChromosomeResultsFrame: Chromosome-specific    #
# ======================================================= #

# --- External Libraries ---
from customtkinter import *
from PIL import Image
import os

# --- Internal Modules ---
from single_marker_test import get_best_marker_info


# ======================================================= #
#        SMTConcatenatedResultsFrame                      #
# ======================================================= #

class SMTConcatenatedResultsFrame(CTkFrame):
    """
    Displays genome-wide SMT -log10(p) plots for each trait.
    Shows the best QTL (highest significance) across the genome.
    """

    def __init__(self, master, traits, all_smt_results, genetic_maps_unfiltered, styles, on_back):
        """
        Initialize the concatenated SMT results viewer.

        Args:
            master (Tk): Parent window.
            traits (dict): Trait dictionary.
            all_smt_results (dict): SMT result dict per trait.
            genetic_maps_unfiltered (dict): Map of physical positions.
            styles (dict): UI style dict.
            on_back (func): Callback to return to previous frame.
        """
        super().__init__(master, width=1100, height=650)
        self.master = master
        self.traits = traits
        self.all_smt_results = all_smt_results
        self.genetic_maps_unfiltered = genetic_maps_unfiltered
        self.styles = styles
        self.on_back = on_back

        self.image_label = None
        self.qtl_label = None

        self.pack(fill="both", expand=True)
        self.build_ui()

    def build_ui(self):
        """Construct layout and trait dropdown."""
        CTkLabel(self, text="Overall SMT Results - Locate Best QTL", font=("Segoe UI", 22, "bold")).pack(pady=(10, 5))

        CTkLabel(self, text="Pick a trait from the list:", font=("Segoe UI", 14)).pack(pady=(0, 2))
        self.trait_dropdown = CTkOptionMenu(
            self,
            values=list(self.traits.keys()),
            width=300,
            fg_color="#ffffff",
            button_color="#ffffff",
            text_color="#000000",
            font=("Segoe UI", 16, "bold"),
            command=self.show_overall_plot
        )
        self.trait_dropdown.pack(pady=(0, 5))

        self.qtl_label = CTkLabel(self, text="", font=("Segoe UI", 14, "bold"))
        self.qtl_label.pack(pady=(5, 5))

        self.image_label = CTkLabel(self, text="")
        self.image_label.pack(padx=20, pady=(5, 5), fill="both", expand=True)

        CTkButton(
            self,
            text="Back",
            font=("Segoe UI", 16, "bold"),
            width=160,
            height=46,
            command=self.on_back,
            **self.styles["white"]
        ).pack(pady=(5, 10))

    def show_overall_plot(self, trait):
        """
        Display genome-wide SMT result image for the selected trait,
        and highlight the most significant QTL location.
        """
        image_path = f"Results/SMT/Plots/smt_concatenated_{trait}.jpg"
        self.display_image(image_path)

        # Get best marker info: (chr, cm, bp, -log10(p))
        chr, cm, bp, logp = get_best_marker_info(trait, self.all_smt_results, self.genetic_maps_unfiltered)
        self.qtl_label.configure(
            text=f"Best QTL for '{trait}': Chromosome = {chr}, cM = {cm}, BP = {bp}, -log10(p) = {logp:.2f}"
        )

    def display_image(self, path):
        """
        Load and display image if path exists, otherwise show warning.

        Args:
            path (str): Path to image file.
        """
        if os.path.exists(path):
            img = Image.open(path)
            resized = img.resize((880, 400))  # Scaled for clarity
            tk_img = CTkImage(light_image=resized, size=resized.size)
            self.image_label.configure(image=tk_img, text="")
            self.image_label.image = tk_img  # Prevent GC
        else:
            self.image_label.configure(image=None, text="Image not found")


# ======================================================= #
#        SMTPerChromosomeResultsFrame                     #
# ======================================================= #

class SMTPerChromosomeResultsFrame(CTkFrame):
    """
    Displays SMT -log10(p) plots for each chromosome and trait.
    Enables selection of trait and chromosome to view the per-chromosome p-value curve.
    """

    def __init__(self, master, traits, styles, on_back):
        """
        Initialize the per-chromosome SMT results viewer.

        Args:
            master (Tk): Parent window.
            traits (dict): Dictionary of available traits.
            styles (dict): Button styles.
            on_back (func): Callback to return to previous frame.
        """
        super().__init__(master, width=1100, height=650)
        self.master = master
        self.traits = traits
        self.styles = styles
        self.on_back = on_back

        self.trait_dropdown = None
        self.chrom_dropdown = None
        self.image_label = None

        self.pack(fill="both", expand=True)
        self.build_ui()

    def build_ui(self):
        """Construct the layout for trait/chromosome selection and image display."""
        CTkLabel(self, text="SMT p-value plots per chromosome", font=("Segoe UI", 22, "bold")).pack(pady=(20, 10))

        # === Dropdowns section ===
        frame = CTkFrame(self, fg_color="transparent")
        frame.pack(pady=5)

        # --- Trait Dropdown ---
        left_frame = CTkFrame(frame, fg_color="transparent")
        left_frame.pack(side="left", padx=10)

        CTkLabel(left_frame, text="Pick a trait from the list:", font=("Segoe UI", 14)).pack()
        self.trait_dropdown = CTkOptionMenu(
            left_frame,
            values=list(self.traits.keys()),
            width=300,
            fg_color="#ffffff",
            button_color="#ffffff",
            text_color="#000000",
            font=("Segoe UI", 16, "bold"),
            command=self.update_plot
        )
        self.trait_dropdown.pack()

        # --- Chromosome Dropdown ---
        right_frame = CTkFrame(frame, fg_color="transparent")
        right_frame.pack(side="left", padx=10)

        CTkLabel(right_frame, text="Pick a chromosome from the list:", font=("Segoe UI", 14)).pack()
        self.chrom_dropdown = CTkOptionMenu(
            right_frame,
            values=[f"chr{str(i).zfill(2)}" for i in range(1, 27)],
            width=200,
            fg_color="#ffffff",
            button_color="#ffffff",
            text_color="#000000",
            font=("Segoe UI", 16, "bold"),
            command=self.update_plot
        )
        self.chrom_dropdown.pack()

        # === Image Output ===
        self.image_label = CTkLabel(self, text="")
        self.image_label.pack(expand=True, fill="both", pady=(10, 10))

        CTkButton(
            self,
            text="Back",
            font=("Segoe UI", 16, "bold"),
            width=160,
            height=46,
            command=self.on_back,
            **self.styles["white"]
        ).pack(pady=(10, 10))

    def update_plot(self, *_):
        """
        Display SMT result image for the selected trait and chromosome.
        """
        trait = self.trait_dropdown.get()
        chrom = self.chrom_dropdown.get().replace("chr", "")
        image_path = f"Results/SMT/Pvalue_graphs/{trait}/smt-chr{chrom}.jpg"
        self.display_image(image_path)

    def display_image(self, path):
        """
        Load and display image from the specified path.

        Args:
            path (str): Image file path to be displayed.
        """
        if os.path.exists(path):
            img = Image.open(path)
            resized = img.resize((880, 400))  # Match width for consistency
            tk_img = CTkImage(light_image=resized, size=(880, 400))
            self.image_label.configure(image=tk_img, text="")
            self.image_label.image = tk_img
        else:
            self.image_label.configure(image=None, text=f"Image not found: {os.path.basename(path)}")
