# =================================================== #
#  SMT_SIM_CompareFrame                               #
#  GUI to compare SMT and SIM results side-by-side   #
# =================================================== #

# --- External Libraries ---
import os
from customtkinter import *
from PIL import Image


# --- Utility Function ---
def load_and_resize_image(path, size=(500, 380)):
    """
    Load an image from the given path and resize it.

    Args:
        path (str): Path to the image file.
        size (tuple): Desired (width, height) of the output image.

    Returns:
        CTkImage: Resized image if successful, None otherwise.
    """
    try:
        image = Image.open(path)
        image = image.resize(size, Image.Resampling.LANCZOS)
        return CTkImage(light_image=image, size=size)
    except Exception as e:
        print(f"Failed to load image {path}: {e}")
        return None


# --- Main Frame Class ---
class SMT_SIM_CompareFrame(CTkFrame):
    """
    Frame to compare SMT vs SIM results visually for a selected trait and chromosome.
    Displays both p-value and LOD score images side-by-side.
    """

    def __init__(self, master, traits, on_back):
        """
        Initialize the comparison frame.

        Args:
            master (Tk): Parent window.
            traits (dict): Available traits.
            on_back (function): Callback to return to previous frame.
        """
        super().__init__(master, fg_color="#1e1e1e")
        self.master = master
        self.traits = traits
        self.on_back = on_back

        self.selected_trait = None
        self.selected_chrom = None

        self.pack(fill="both", expand=True)
        self.center_window(1100, 700)
        self.build_ui()

    def center_window(self, width, height):
        """Center the window on screen."""
        screen_width = self.master.winfo_screenwidth()
        screen_height = self.master.winfo_screenheight()
        x = int((screen_width - width) / 2)
        y = int((screen_height - height) / 2) - 60  # Slight upward offset
        self.master.geometry(f"{width}x{height}+{x}+{y}")

    def build_ui(self):
        """Construct the layout with dropdowns, image zones, and back button."""
        CTkLabel(
            self,
            text="Comparison of SMT vs SIM outputs",
            font=("Arial", 24, "bold"),
            text_color="white"
        ).place(relx=0.5, rely=0.07, anchor="center")

        # === Dropdown Layout ===
        spacing = 20
        dropdown_width = 200
        label_width = 210
        total_width = label_width + dropdown_width + spacing + label_width + dropdown_width
        start_x = (1100 - total_width) // 2

        # Trait Selector
        CTkLabel(self, text="Pick a trait from list:", font=("Arial", 16), text_color="white")\
            .place(x=start_x + 50, y=90)
        self.trait_dropdown = CTkOptionMenu(
            self,
            values=list(self.traits.keys()),
            command=self.update_trait,
            width=dropdown_width
        )
        self.trait_dropdown.place(x=start_x + label_width, y=90)

        # Chromosome Selector
        CTkLabel(self, text="Pick a chromosome from list:", font=("Arial", 16), text_color="white")\
            .place(x=start_x + label_width + dropdown_width + spacing, y=90)
        self.chrom_dropdown = CTkOptionMenu(
            self,
            values=[f"chr{str(i).zfill(2)}" for i in range(1, 11)],  # Limit to chr01–chr10 for now
            command=self.update_chrom,
            width=dropdown_width
        )
        self.chrom_dropdown.place(x=start_x + label_width * 2 + dropdown_width + spacing, y=90)

        # === Image Display Areas ===
        self.left_image_label = CTkLabel(self, text="")  # SMT image
        self.left_image_label.place(x=40, y=180)

        self.right_image_label = CTkLabel(self, text="")  # SIM image
        self.right_image_label.place(x=570, y=180)

        # === Back Button ===
        CTkButton(
            self,
            text="Back",
            command=self.go_back,
            font=("Arial", 14, "bold"),
            fg_color="white",
            text_color="black",
            hover_color="#dddddd",
            width=100
        ).place(relx=0.5, y=670, anchor="center")

    def update_trait(self, value):
        """Called when a new trait is selected."""
        self.selected_trait = value
        self.try_update_images()

    def update_chrom(self, value):
        """Called when a new chromosome is selected."""
        self.selected_chrom = value
        self.try_update_images()

    def try_update_images(self):
        """Display SMT and SIM images if both trait and chromosome are selected."""
        if self.selected_trait and self.selected_chrom:
            smt_path = f"Results/SMT/Pvalue_graphs/{self.selected_trait}/smt-{self.selected_chrom}.jpg"
            sim_path = f"Results/SIM/Plots/{self.selected_trait}/{self.selected_chrom}_LOD_curve.png"

            smt_img = load_and_resize_image(smt_path)
            if smt_img:
                self.left_image_label.configure(image=smt_img, text="")
                self.left_image_label.image = smt_img

            sim_img = load_and_resize_image(sim_path)
            if sim_img:
                self.right_image_label.configure(image=sim_img, text="")
                self.right_image_label.image = sim_img

    def go_back(self):
        """Return to the previous frame."""
        self.destroy()
        self.on_back()
