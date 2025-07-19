from customtkinter import *
from PIL import Image as PILImage
import os

class GeneticMapViewerFrame(CTkFrame):
    def __init__(self, master, vcf_data, styles, on_next):
        super().__init__(master, width=1100, height=650)
        self.master = master
        self.vcf_data = vcf_data
        self.styles = styles
        self.on_next = on_next

        self.filtered_selected = StringVar(value="unfiltered")
        self.chrom_selected = StringVar()
        self.image_label = None
        self.chromosomes = [c for c in vcf_data if c.startswith("chr")]

        self.pack(fill="both", expand=True)
        self.build_ui()

    def build_ui(self):
        CTkLabel(self, text="View Genetic Maps", font=("Segoe UI", 22, "bold")).pack(pady=(40, 20))

        selection_frame = CTkFrame(self, fg_color="transparent")
        selection_frame.pack(pady=10)

        CTkOptionMenu(selection_frame, values=["unfiltered", "filtered"],
                      variable=self.filtered_selected, command=self.update_chromosome_menu).grid(row=0, column=0, padx=10)

        self.chrom_menu = CTkOptionMenu(selection_frame, values=self.chromosomes,
                                        variable=self.chrom_selected, command=self.display_map)
        self.chrom_menu.grid(row=0, column=1, padx=10)

        if self.chromosomes:
            self.chrom_selected.set(self.chromosomes[0])
            self.display_map(self.chromosomes[0])

        self.image_label = CTkLabel(self, text="")
        self.image_label.pack(pady=(20, 10))

        # NEXT BUTTON (bottom right)
        CTkButton(self, text="Next", font=("Segoe UI", 16, "bold"),
                  width=180, height=48, command=self.on_next, **self.styles["red"]).place(
            relx=1.0, rely=1.0, x=-30, y=-30, anchor="se"
        )

    def update_chromosome_menu(self, *_):
        self.chrom_menu.configure(values=self.chromosomes)
        if self.chromosomes:
            self.chrom_selected.set(self.chromosomes[0])
            self.display_map(self.chromosomes[0])

    def display_map(self, chrom):
        if not chrom:
            return

        mode = self.filtered_selected.get()
        prefix = "filtered-" if mode == "filtered" else ""
        folder = "Genetic_Maps_Filtered" if mode == "filtered" else "Genetic_Maps"
        filename = f"{prefix}genetic-{chrom}.jpg"
        path = os.path.join("Results", folder, filename)

        if os.path.exists(path):
            image = PILImage.open(path).resize((900, 300))
            ctk_image = CTkImage(light_image=image, size=(900, 300))
            self.image_label.configure(image=ctk_image, text="")
            self.image_label.image = ctk_image
        else:
            self.image_label.configure(image=None, text="Map image not found.")
