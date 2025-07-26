from customtkinter import *
from PIL import Image as PILImage
import os

class GeneticMapViewerFrame(CTkFrame):
    def __init__(self, master, vcf_data, styles, on_back, mode):
        super().__init__(master, width=1100, height=850)
        self.master = master
        self.vcf_data = vcf_data or {}
        self.styles = styles
        self.on_back = on_back
        self.mode = mode

        self.chrom_selected = StringVar()
        self.image_label = None
        self.chromosomes = self.get_chromosomes_from_vcf(self.vcf_data)

        self.pack(fill="both", expand=True)
        self.build_ui()

    def get_chromosomes_from_vcf(self, vcf_data):
        return [chrom for chrom in vcf_data.keys() if chrom.startswith("chr")]

    def build_ui(self):
        self.title_label = CTkLabel(self, text="", font=("Segoe UI", 22, "bold"))
        self.title_label.pack(pady=(40, 20))

        # Check if VCF data is missing
        if not self.vcf_data:
            CTkLabel(self, text="No VCF data available.",
                     font=("Segoe UI", 16), text_color="red").pack(pady=20)
            CTkButton(self, text="Back", font=("Segoe UI", 14, "bold"),
                      width=140, height=40, command=self.on_back,
                      **self.styles["white"]).place(x=20, y=50, anchor="nw")
            return

        # Check if no chromosomes found
        if not self.chromosomes:
            CTkLabel(self, text="No valid chromosomes found in VCF data.",
                     font=("Segoe UI", 16), text_color="red").pack(pady=20)
            CTkButton(self, text="Back", font=("Segoe UI", 14, "bold"),
                      width=140, height=40, command=self.on_back,
                      **self.styles["white"]).place(x=20, y=50, anchor="nw")
            return

        selection_frame = CTkFrame(self, fg_color="transparent")
        selection_frame.pack(pady=10)

        CTkLabel(selection_frame, text="Pick Chromosome From List", font=("Segoe UI", 16, "bold")).grid(row=0, column=0, padx=10, pady=(0, 10))

        dropdown_style = {
            "fg_color": "#ffffff",
            "button_color": "#ffffff",
            "text_color": "#000000",
            "font": ("Segoe UI", 16, "bold")
        }

        self.chrom_menu = CTkOptionMenu(
            selection_frame,
            values=self.chromosomes,
            variable=self.chrom_selected,
            command=self.display_map,
            **dropdown_style
        )
        self.chrom_menu.grid(row=1, column=0, padx=10)

        self.image_label = CTkLabel(self, text="")
        self.image_label.pack(pady=(20, 10))

        if self.chromosomes:
            self.chrom_selected.set(self.chromosomes[0])
            self.display_map(self.chromosomes[0])

        CTkButton(self, text="Back", font=("Segoe UI", 14, "bold"),
                  width=140, height=40, command=self.on_back,
                  **self.styles["white"]).place(x=20, y=50, anchor="nw")

    def display_map(self, chrom):
        if not chrom:
            return

        filename = ""
        if self.mode == "heatmap":
            folder = "Results/HeatMaps"
            filename = f"heatmap-{chrom}.jpg"
            self.title_label.configure(text="Recombination Rate Heat Map")
        elif self.mode == "compare-filtered":
            folder = "Results/CompareOfDistances_filtered"
            filename = f"compare-filtered-{chrom}.jpg"
            self.title_label.configure(text="Filtered Genetic Map Comparison")
        elif self.mode == "compare":
            folder = "Results/CompareOfDistances_unfiltered"
            filename = f"compare-{chrom}.jpg"
            self.title_label.configure(text="Unfiltered Genetic Map Comparison")
        else:  # default to "genetic"
            folder = "Results/Genetic_Maps"
            filename = f"genetic-{chrom}.jpg"
            self.title_label.configure(text="Chromosome-wide Genetic Map")

        path = os.path.join(folder, filename)

        if os.path.exists(path):
            image = PILImage.open(path).resize((1000, 450))
            ctk_image = CTkImage(light_image=image, size=(1000, 450))
            self.image_label.configure(image=ctk_image, text="")
            self.image_label.image = ctk_image
        else:
            self.image_label.configure(image=None, text=f"Image not found: {filename}")
