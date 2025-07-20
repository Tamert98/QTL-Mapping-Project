from customtkinter import *
from PIL import Image
import os
import threading
import time
from vcf_data_handler import (
    parse_vcf_file,
    parse_trait_file,
    generate_genetic_map,
    generate_genetic_map_filtered,
    generate_pairwise_genetic_map
)
from plot_utils import (
    compute_global_physical_max,
    generate_genetic_map_images_and_pdf,
    generate_filtered_genetic_map_images_and_pdf,
    generate_heatmap_images_and_pdf,
    print_heatmap_pdf_path
)

class LoadingDataFrame(CTkFrame):
    def __init__(self, master, vcf_path, trait_path, on_done, on_back, styles):
        super().__init__(master, width=1100, height=650)
        self.master = master
        self.vcf_path = vcf_path
        self.trait_path = trait_path
        self.on_done = on_done
        self.on_back = on_back
        self.styles = styles

        self.vcf_data = None
        self.sample_names = None
        self.traits = None
        self.unfiltered_maps = None
        self.filtered_maps = None

        self.pack(fill="both", expand=True)
        self.data_ready = False
        self.next_btn = None
        self.maps_saved_label = None

        self.build_ui()
        threading.Thread(target=self.load_and_generate, daemon=True).start()

    def build_ui(self):
        # Center window
        width, height = 1100, 650
        self.master.geometry(f"{width}x{height}")
        self.master.update_idletasks()
        screen_width = self.master.winfo_screenwidth()
        screen_height = self.master.winfo_screenheight()
        x = (screen_width // 2) - (width // 2)
        y = (screen_height // 2) - (height // 2) - 100
        self.master.geometry(f"{width}x{height}+{x}+{y}")

        # Side image
        try:
            side_image = CTkImage(light_image=Image.open("sideimg.png"), size=(250, 650))
            CTkLabel(self, image=side_image, text="").place(x=0, y=0)
            self.side_image = side_image  # prevent garbage collection
        except:
            print("Side image not found or invalid format.")

        # Right content area (visually centered)
        content_frame = CTkFrame(self, fg_color="transparent")
        content_frame.place(relx=0.62, rely=0.5, anchor="center")

        CTkLabel(content_frame, text="Loading Data", font=("Segoe UI", 22, "bold")).pack(pady=(10, 10))

        intro_text = (
            "We are now loading and saving your VCF and Trait data.\n"
            "This process will prepare genetic maps and heatmaps used in the QTL analysis pipeline."
        )
        CTkLabel(content_frame, text=intro_text, font=("Segoe UI", 14), justify="left", wraplength=650).pack(pady=5)

        self.vcf_label = CTkLabel(content_frame, text="Loading VCF data...", font=("Segoe UI", 18, "bold"))
        self.vcf_label.pack(pady=(20, 5), anchor="w")

        self.trait_label = CTkLabel(content_frame, text="Waiting to load trait data...", font=("Segoe UI", 18, "bold"))
        self.trait_label.pack(pady=(10, 5), anchor="w")

        self.maps_label = CTkLabel(content_frame, text="Waiting to generate genetic maps and heatmaps...", font=("Segoe UI", 18, "bold"))
        self.maps_label.pack(pady=(20, 10), anchor="w")

        self.maps_saved_label = CTkLabel(content_frame, text="", font=("Segoe UI", 14, "bold"), text_color="#228B22")
        self.maps_saved_label.pack(pady=(10, 10), anchor="w")

        self.output_label = CTkLabel(content_frame, text="", font=("Segoe UI", 13), justify="left")
        self.output_label.pack(pady=(5, 20), anchor="w")

        nav_frame = CTkFrame(content_frame, fg_color="transparent")
        nav_frame.pack(pady=10)

        self.back_btn = CTkButton(nav_frame, text="Back", font=("Segoe UI", 14, "bold"),
                                  width=120, command=self.on_back, **self.styles["white"])
        self.back_btn.pack(side="left", padx=(0, 20))

        self.next_btn = CTkButton(nav_frame, text="Next", font=("Segoe UI", 14, "bold"),
                                  width=120, state="disabled", command=self.on_next_clicked, **self.styles["red"])
        self.next_btn.pack(side="left")

    def load_and_generate(self):
        self.vcf_data, self.sample_names = parse_vcf_file(self.vcf_path)
        self.vcf_label.configure(text="✔ VCF data has been loaded", text_color=self.styles["red"]["fg_color"])

        time.sleep(0.5)

        _, self.traits = parse_trait_file(self.trait_path)
        self.trait_label.configure(text="✔ Trait data has been loaded", text_color=self.styles["red"]["fg_color"])

        time.sleep(0.5)

        self.unfiltered_maps = generate_genetic_map(self.vcf_data, self.sample_names)
        self.filtered_maps = generate_genetic_map_filtered(self.vcf_data, self.sample_names)

        global_xmax = compute_global_physical_max(self.unfiltered_maps)
        generate_genetic_map_images_and_pdf(self.unfiltered_maps, global_xmax)
        generate_filtered_genetic_map_images_and_pdf(self.filtered_maps, global_xmax)

        # Optional heatmaps
        # pairwise_maps = {}
        # for chrom_id in self.filtered_maps:
        #     pairwise = generate_pairwise_genetic_map(self.vcf_data, self.sample_names, chrom_id)
        #     if pairwise:
        #         pairwise_maps[chrom_id] = pairwise
        # generate_heatmap_images_and_pdf(pairwise_maps)

        print_heatmap_pdf_path(output_dir="Results/Genetic_Maps")

        self.maps_label.configure(text="✔ Genetic maps and heatmaps have been generated", text_color=self.styles["red"]["fg_color"])
        self.maps_saved_label.configure(text="✔ Genetic maps and heatmaps have been saved into the Results folder.")

        self.output_label.configure(text=(
            "Genetic map PDF saved at: Results/Genetic_Maps/final_genetic_maps.pdf\n"
            "Filtered genetic map PDF saved at: Results/Genetic_Maps_Filtered/final_filtered_genetic_maps.pdf\n"
            "Heatmap PDF saved at: Results/HeatMaps/final_heatmaps.pdf"
        ))

        self.data_ready = True
        self.next_btn.configure(state="normal")

    def on_next_clicked(self):
        self.next_btn.configure(state="disabled")
        self.after(1000, self.finish_and_continue)

    def finish_and_continue(self):
        self.on_done(
            self.vcf_data,
            self.sample_names,
            self.traits,
            self.unfiltered_maps,
            self.filtered_maps
        )
