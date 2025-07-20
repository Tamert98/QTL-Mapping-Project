from customtkinter import *
from PIL import Image as PILImage
import os

class GeneticMapViewerFrame(CTkFrame):
    def __init__(self, master, vcf_data, styles, on_back):
        super().__init__(master, width=1100, height=850)
        self.master = master
        self.vcf_data = vcf_data
        self.styles = styles
        self.on_back = on_back

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

        # Updated to include 'heatmap' option


        # Use the same red style logic as in SMTViewResultsFrame
        red_fg = self.styles["red"].get("fg_color", "#B22222")
        red_btn = self.styles["red"].get("button_color", red_fg)

        CTkOptionMenu(
            selection_frame,
            values=["unfiltered", "filtered", "heatmap"],
            variable=self.filtered_selected,
            command=self.update_chromosome_menu,
            fg_color=red_fg,
            button_color=red_fg
        ).grid(row=0, column=0, padx=10)

        self.chrom_menu = CTkOptionMenu(
            selection_frame,
            values=self.chromosomes,
            variable=self.chrom_selected,
            command=self.display_map,
            fg_color=red_fg,
            button_color=red_fg
        )
        self.chrom_menu.grid(row=0, column=1, padx=10)

        self.image_label = CTkLabel(self, text="")
        self.image_label.pack(pady=(20, 10))

        if self.chromosomes:
            self.chrom_selected.set(self.chromosomes[0])
            self.display_map(self.chromosomes[0])

        # Navigation button in upper left, but lower down
        CTkButton(self, text="Back", font=("Segoe UI", 14, "bold"),
                  width=140, height=40, command=self.on_back, **self.styles["white"]).place(x=20, y=50, anchor="nw")

    def update_chromosome_menu(self, *_):
        self.chrom_menu.configure(values=self.chromosomes)
        if self.chromosomes:
            self.chrom_selected.set(self.chromosomes[0])
            self.display_map(self.chromosomes[0])

    def display_map(self, chrom):
        if not chrom:
            return

        mode = self.filtered_selected.get()

        if mode == "heatmap":
            folder = "Results/HeatMaps"
            filename = f"heatmap-{chrom}.jpg"
        else:
            prefix = "filtered-" if mode == "filtered" else ""
            folder = "Results/Genetic_Maps_Filtered" if mode == "filtered" else "Results/Genetic_Maps"
            filename = f"{prefix}genetic-{chrom}.jpg"

        path = os.path.join(folder, filename)

        if os.path.exists(path):
            image = PILImage.open(path).resize((1000, 450))  # Enlarged image area
            ctk_image = CTkImage(light_image=image, size=(1000, 450))
            self.image_label.configure(image=ctk_image, text="")
            self.image_label.image = ctk_image
        else:
            self.image_label.configure(image=None, text="Map image not found.")
