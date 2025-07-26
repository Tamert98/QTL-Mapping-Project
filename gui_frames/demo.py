import os
from customtkinter import *
from PIL import Image

# === Utility: load and resize ===
def load_and_resize_image(path, size=(500, 380)):
    try:
        image = Image.open(path)
        image = image.resize(size, Image.Resampling.LANCZOS)
        return CTkImage(light_image=image, size=size)
    except Exception as e:
        print(f"Failed to load image {path}: {e}")
        return None

class SMT_SIM_CompareFrame(CTkFrame):
    def __init__(self, master, traits, on_back):
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
        screen_width = self.master.winfo_screenwidth()
        screen_height = self.master.winfo_screenheight()
        x = int((screen_width - width) / 2)
        y = int((screen_height - height) / 2) - 60  # Move window slightly upward
        self.master.geometry(f"{width}x{height}+{x}+{y}")

    def build_ui(self):
        CTkLabel(
            self,
            text="Comparison of SMT vs SIM outputs",
            font=("Arial", 24, "bold"),
            text_color="white"
        ).place(relx=0.5, rely=0.07, anchor="center")

        # Centered dropdown layout
        spacing = 20
        dropdown_width = 200
        label_width = 210
        total_width = label_width + dropdown_width + spacing + label_width + dropdown_width
        start_x = (1100 - total_width) // 2

        CTkLabel(self, text="Pick a trait from list:", font=("Arial", 16), text_color="white")\
            .place(x=start_x + 50, y=90)
        self.trait_dropdown = CTkOptionMenu(
            self,
            values=list(self.traits.keys()),
            command=self.update_trait,
            width=dropdown_width
        )
        self.trait_dropdown.place(x=start_x + label_width, y=90)

        CTkLabel(self, text="Pick a chromosome from list:", font=("Arial", 16), text_color="white")\
            .place(x=start_x + label_width + dropdown_width + spacing, y=90)
        self.chrom_dropdown = CTkOptionMenu(
            self,
            values=["chr01", "chr02", "chr03", "chr04", "chr05", "chr06", "chr07", "chr08", "chr09", "chr10"],
            command=self.update_chrom,
            width=dropdown_width
        )
        self.chrom_dropdown.place(x=start_x + label_width*2 + dropdown_width + spacing, y=90)

        # Image zones
        self.left_image_label = CTkLabel(self, text="")
        self.left_image_label.place(x=40, y=180)

        self.right_image_label = CTkLabel(self, text="")
        self.right_image_label.place(x=570, y=180)

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
        self.selected_trait = value
        self.try_update_images()

    def update_chrom(self, value):
        self.selected_chrom = value
        self.try_update_images()

    def try_update_images(self):
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
        self.destroy()
        self.on_back()
