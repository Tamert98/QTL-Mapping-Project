from customtkinter import *
from PIL import Image
import os
import threading
from single_marker_test import run_smt_for_all_traits, get_best_marker_info
from plot_utils import generate_concatenated_qtl_pdf

class SMTResultsFrame(CTkFrame):
    def __init__(self, master, vcf_data, traits, sample_names, genetic_maps_unfiltered, styles):
        super().__init__(master, width=1100, height=650)
        self.master = master
        self.vcf_data = vcf_data
        self.traits = traits
        self.sample_names = sample_names
        self.genetic_maps_unfiltered = genetic_maps_unfiltered
        self.styles = styles

        self.all_smt_results = None
        self.trait_dropdown = None
        self.image_label = None
        self.qtl_label = None

        self.pack(fill="both", expand=True)
        self.build_loading_ui()
        threading.Thread(target=self.run_smt_thread, daemon=True).start()

    def build_loading_ui(self):
        CTkLabel(self, text="Applying SMT", font=("Segoe UI", 22, "bold")).pack(pady=(30, 10))
        self.status_label = CTkLabel(self, text="The first algorithm is being applied to your data...", font=("Segoe UI", 16))
        self.status_label.pack(pady=10)

        self.plots_label = CTkLabel(self, text="✔ The plots have been saved in Results/SMT/Plots",
                                    font=("Segoe UI", 14), text_color=self.styles["red"]["fg_color"])
        self.reports_label = CTkLabel(self, text="✔ The reports have been saved in Results/SMT/Reports",
                                      font=("Segoe UI", 14), text_color=self.styles["red"]["fg_color"])

        self.plots_label.pack_forget()
        self.reports_label.pack_forget()

    def run_smt_thread(self):
        self.all_smt_results = run_smt_for_all_traits(self.vcf_data, self.traits, self.sample_names)
        generate_concatenated_qtl_pdf(self.all_smt_results)
        self.master.after(3000, self.on_smt_done)  # Delay by 3 seconds

    def on_smt_done(self):
        self.status_label.configure(text="✔ Finished applying SMT", text_color=self.styles["red"]["fg_color"])
        self.plots_label.pack()
        self.reports_label.pack(pady=(0, 15))
        self.after(500, self.build_results_ui)

    def build_results_ui(self):
        for widget in self.winfo_children():
            widget.destroy()

        CTkLabel(self, text="View Trait Results", font=("Segoe UI", 22, "bold")).pack(pady=(20, 10))

        dropdown_frame = CTkFrame(self, fg_color="transparent")
        dropdown_frame.pack(pady=10)

        self.trait_dropdown = CTkOptionMenu(dropdown_frame, values=list(self.traits.keys()), width=220)
        self.trait_dropdown.pack(side="left", padx=10)

        CTkButton(dropdown_frame, text="Find Best QTL", font=("Segoe UI", 14, "bold"),
                  command=self.show_trait_plot, **self.styles["red"]).pack(side="left", padx=10)

        self.image_label = CTkLabel(self, text="")
        self.image_label.pack(expand=True, fill="both", pady=(10, 10))

        self.qtl_label = CTkLabel(self, text="", font=("Segoe UI", 14, "bold"))
        self.qtl_label.pack(pady=(5, 15))

        CTkButton(self, text="Next", font=("Segoe UI", 16, "bold"), width=160, height=46,
                  **self.styles["red"]).pack(pady=(10, 10))

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
