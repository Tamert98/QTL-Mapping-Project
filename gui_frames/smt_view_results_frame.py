from customtkinter import *
from PIL import Image
import os
from single_marker_test import get_best_marker_info

class SMTViewResultsFrame(CTkFrame):
    def __init__(self, master, traits, all_smt_results, genetic_maps_unfiltered, styles, on_back):
        super().__init__(master, width=1100, height=650)
        self.master = master
        self.traits = traits
        self.all_smt_results = all_smt_results
        self.genetic_maps_unfiltered = genetic_maps_unfiltered
        self.styles = styles
        self.on_back = on_back

        self.trait_dropdown = None
        self.image_label = None
        self.qtl_label = None

        self.pack(fill="both", expand=True)
        self.build_ui()

    def build_ui(self):
        CTkLabel(self, text="View SMT Results", font=("Segoe UI", 22, "bold")).pack(pady=(20, 10))


        CTkLabel(self, text="Pick a trait from list, then click on Find Best QTL", font=("Segoe UI", 14, "bold")).pack(pady=(10, 0))
        dropdown_frame = CTkFrame(self, fg_color="transparent")
        dropdown_frame.pack(pady=10)


        # Style the dropdown in red (same as red buttons)
        self.trait_dropdown = CTkOptionMenu(
            dropdown_frame,
            values=list(self.traits.keys()),
            width=220,
            fg_color=self.styles["red"].get("fg_color", "#B22222"),
            button_color=self.styles["red"].get("fg_color", "#B22222"),
            button_hover_color=self.styles["red"].get("hover_color", "#8B1A1A")
        )
        self.trait_dropdown.pack(side="left", padx=10)

        CTkButton(dropdown_frame, text="Find Best QTL", font=("Segoe UI", 14, "bold"),
                  command=self.show_trait_plot, **self.styles["red"]).pack(side="left", padx=10)

        self.image_label = CTkLabel(self, text="")
        self.image_label.pack(expand=True, fill="both", pady=(10, 10))

        self.qtl_label = CTkLabel(self, text="", font=("Segoe UI", 14, "bold"))
        self.qtl_label.pack(pady=(5, 15))

        CTkButton(self, text="Back to Menu", font=("Segoe UI", 16, "bold"), width=160, height=46,
                  command=self.on_back, **self.styles["white"]).pack(pady=(10, 10))

    def show_trait_plot(self):
        trait = self.trait_dropdown.get()
        image_path = f"Results/SMT/Plots/smt_concatenated_{trait}.jpg"

        if os.path.exists(image_path):
            original_img = Image.open(image_path)
            screen_width = self.master.winfo_screenwidth()
            max_width = min(screen_width - 100, 1100)
            aspect_ratio = original_img.height / original_img.width
            new_height = int(max_width * aspect_ratio)
            resized_img = original_img.resize((max_width, new_height))

            img = CTkImage(light_image=resized_img, size=(max_width, new_height))
            self.image_label.configure(image=img, text="")
            self.image_label.image = img

        chr, cm, bp, logp = get_best_marker_info(trait, self.all_smt_results, self.genetic_maps_unfiltered)
        self.qtl_label.configure(
            text=f"Best QTL for {trait}:\nChromosome = {chr}, cM = {cm}, BP = {bp}, -log10(p) = {logp:.2f}"
        )
