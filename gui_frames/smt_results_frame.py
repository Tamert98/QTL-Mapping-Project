from customtkinter import *
from PIL import Image
import os
from single_marker_test import get_best_marker_info

class SMTConcatenatedResultsFrame(CTkFrame):
    def __init__(self, master, traits, all_smt_results, genetic_maps_unfiltered, styles, on_back):
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

        CTkButton(self, text="Back", font=("Segoe UI", 16, "bold"), width=160, height=46,
                  command=self.on_back, **self.styles["white"]).pack(pady=(5, 10))

    def show_overall_plot(self, trait):
        image_path = f"Results/SMT/Plots/smt_concatenated_{trait}.jpg"
        self.display_image(image_path)

        chr, cm, bp, logp = get_best_marker_info(trait, self.all_smt_results, self.genetic_maps_unfiltered)
        self.qtl_label.configure(
            text=f"Best QTL for '{trait}': Chromosome = {chr}, cM = {cm}, BP = {bp}, -log10(p) = {logp:.2f}"
        )

    def display_image(self, path):
        if os.path.exists(path):
            img = Image.open(path)
            resized = img.resize((880, 400))  # ✅ narrowed width for better fit
            tk_img = CTkImage(light_image=resized, size=resized.size)
            self.image_label.configure(image=tk_img, text="")
            self.image_label.image = tk_img
        else:
            self.image_label.configure(image=None, text="Image not found")


class SMTPerChromosomeResultsFrame(CTkFrame):
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
        CTkLabel(self, text="SMT p-value plots per chromosome", font=("Segoe UI", 22, "bold")).pack(pady=(20, 10))

        frame = CTkFrame(self, fg_color="transparent")
        frame.pack(pady=5)

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

        self.image_label = CTkLabel(self, text="")
        self.image_label.pack(expand=True, fill="both", pady=(10, 10))

        CTkButton(self, text="Back", font=("Segoe UI", 16, "bold"), width=160, height=46,
                  command=self.on_back, **self.styles["white"]).pack(pady=(10, 10))

    def update_plot(self, *_):
        trait = self.trait_dropdown.get()
        chrom = self.chrom_dropdown.get().replace("chr", "")
        image_path = f"Results/SMT/Pvalue_graphs/{trait}/smt-chr{chrom}.jpg"
        self.display_image(image_path)

    def display_image(self, path):
        if os.path.exists(path):
            img = Image.open(path)
            resized = img.resize((880, 400))  # Match width from above
            tk_img = CTkImage(light_image=resized, size=(880, 400))
            self.image_label.configure(image=tk_img, text="")
            self.image_label.image = tk_img
        else:
            self.image_label.configure(image=None, text=f"Image not found: {os.path.basename(path)}")
