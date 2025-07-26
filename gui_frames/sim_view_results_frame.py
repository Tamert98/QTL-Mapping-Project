from customtkinter import *
from PIL import Image
import os

class SIMViewResultsFrame(CTkFrame):
    def __init__(self, master, traits, styles, on_back):
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
        CTkLabel(self, text="LOD score curves for each chromosome and trait",
                 font=("Segoe UI", 22, "bold")).pack(pady=(20, 10))

        dropdown_frame = CTkFrame(self, fg_color="transparent")
        dropdown_frame.pack(pady=10)

        # === Trait Dropdown ===
        trait_frame = CTkFrame(dropdown_frame, fg_color="transparent")
        trait_frame.pack(side="left", padx=15)
        CTkLabel(trait_frame, text="Pick a trait from list:", font=("Segoe UI", 14)).pack()
        self.trait_dropdown = CTkOptionMenu(
            trait_frame,
            values=list(self.traits.keys()),
            width=300,
            fg_color="#ffffff",
            button_color="#ffffff",
            text_color="#000000",
            font=("Segoe UI", 16, "bold"),
            command=self.update_plot
        )
        self.trait_dropdown.pack()

        # === Chromosome Dropdown ===
        chrom_frame = CTkFrame(dropdown_frame, fg_color="transparent")
        chrom_frame.pack(side="left", padx=15)
        CTkLabel(chrom_frame, text="Pick a chromosome from list:", font=("Segoe UI", 14)).pack()
        self.chrom_dropdown = CTkOptionMenu(
            chrom_frame,
            values=[f"chr{str(i).zfill(2)}" for i in range(1, 27)],
            width=220,
            fg_color="#ffffff",
            button_color="#ffffff",
            text_color="#000000",
            font=("Segoe UI", 16, "bold"),
            command=self.update_plot
        )
        self.chrom_dropdown.pack()

        # === Image Display ===
        self.image_label = CTkLabel(self, text="")
        self.image_label.pack(expand=True, fill="both", pady=(10, 10), padx=20)

        # === Back Button ===
        CTkButton(self, text="Back", font=("Segoe UI", 16, "bold"), width=160, height=46,
                  command=self.on_back, **self.styles["white"]).pack(pady=(5, 15))

    def update_plot(self, *_):
        trait = self.trait_dropdown.get()
        chrom = self.chrom_dropdown.get()
        image_path = f"Results/SIM/Plots/{trait}/{chrom}_LOD_curve.png"
        self.display_image(image_path)

    def display_image(self, path):
        if os.path.exists(path):
            img = Image.open(path)
            resized = img.resize((880, 400))
            tk_img = CTkImage(light_image=resized, size=resized.size)
            self.image_label.configure(image=tk_img, text="")
            self.image_label.image = tk_img
        else:
            self.image_label.configure(image=None, text=f"Image not found: {os.path.basename(path)}")
