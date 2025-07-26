from customtkinter import *
from PIL import Image
import threading
import time

from single_marker_test import run_smt_for_all_traits
from single_interval_mapping import run_sim_on_selected_markers
from plot_utils import (
    generate_concatenated_qtl_pdf,
    plot_all_lod_curves,
    generate_smt_chromosome_images_all_traits,
    stitch_all_smt_chromosome_images,
)


class SMT_SIM_RunFrame(CTkFrame):
    def __init__(self, master, vcf_data, selected_markers, traits, sample_names, on_done, styles, genetic_maps_unfiltered):
        super().__init__(master, width=1100, height=650)
        self.master = master
        self.vcf_data = vcf_data
        self.selected_markers = selected_markers
        self.traits = traits
        self.sample_names = sample_names
        self.styles = styles
        self.on_done = on_done
        self.genetic_maps_unfiltered = genetic_maps_unfiltered
        self.smt_results = None
        self.sim_results = None

        self.pack(fill="both", expand=True)
        self.message_lines = []
        self.textbox = None
        self.build_ui()

        self.status_label.configure(
            text="The SMT algorithm is being applied...",
            text_color=self.styles["red"]["fg_color"]
        )
        self.update_idletasks()
        threading.Thread(target=self.run_smt_then_sim, daemon=True).start()

    def build_ui(self):
        try:
            image = CTkImage(light_image=Image.open("gui_frames/Step3.png"), size=(250, 650))
            CTkLabel(self, image=image, text="").place(x=0, y=0)
            self.step_image = image
        except Exception as e:
            print(f"Step3.png not found: {e}")

        content_frame = CTkFrame(self, fg_color="transparent")
        content_frame.place(relx=0.6, rely=0.5, anchor="center")

        CTkLabel(content_frame, text="Applying SMT and SIM", font=("Segoe UI", 22, "bold")).pack(pady=(30, 10))
        self.status_label = CTkLabel(content_frame, text="", font=("Segoe UI", 16))
        self.status_label.pack(pady=10)

        self.textbox = CTkTextbox(content_frame, width=800, height=240, font=("Segoe UI", 13), wrap="word")
        self.textbox.pack(pady=20)
        self.textbox.insert("end", "Progress will appear here...\n")
        self.textbox.configure(state="disabled")

        self.final_label = CTkLabel(content_frame, text="", font=("Segoe UI", 14), text_color="#228B22")
        self.final_label.pack()

    def append_message(self, msg):
        self.message_lines.append(msg)

        def update_textbox():
            self.textbox.configure(state="normal")
            self.textbox.insert("end", msg + "\n")
            self.textbox.see("end")
            self.textbox.configure(state="disabled")

        self.master.after(0, update_textbox)

    def run_smt_then_sim(self):
        self.append_message("=== Starting SMT ===")
        self.smt_results = run_smt_for_all_traits(
            self.vcf_data, self.traits, self.sample_names, message_callback=self.append_message
        )

        self.master.after(0, lambda: self.status_label.configure(
            text="Generating SMT plots...", text_color=self.styles["red"]["fg_color"]
        ))
        self.append_message("Generating chromosome-level SMT plots...")
        generate_smt_chromosome_images_all_traits(self.smt_results, send_message=self.append_message)

        self.append_message("Stitching SMT plots into final PDFs...")
        stitch_all_smt_chromosome_images(self.smt_results, send_message=self.append_message)

        self.append_message("Generating concatenated QTL map...")
        generate_concatenated_qtl_pdf(self.smt_results, message_callback=self.append_message)

        self.master.after(0, lambda: self.status_label.configure(
            text="✔ SMT Completed. Starting SIM...", text_color="#228B22"
        ))
        self.append_message("\n=== SMT Completed ===")
        self.append_message("=== Starting SIM ===")

        self.sim_results = run_sim_on_selected_markers(
            self.vcf_data,
            self.selected_markers,
            self.traits,
            self.sample_names,
            message_callback=self.append_message
        )

        self.append_message("Generating SIM LOD score plots...")
        plot_all_lod_curves(
            self.sim_results,
            alpha=0.05,
            output_dir="Results/SIM/Plots",
            send_message=self.append_message  # ✅ Corrected
        )
        self.master.after(0, lambda: self.status_label.configure(
            text="✔ SMT and SIM Completed", text_color="#228B22"
        ))
        self.append_message("=== All Analysis Complete ===")
        self.final_label.configure(text="✔ Results saved in Results/SMT and Results/SIM")

        time.sleep(1)
        self.master.after(0, lambda: self.on_done(
            sim_results=self.sim_results,
            smt_results=self.smt_results,
            vcf_data=self.vcf_data,
            traits=self.traits,
            sample_names=self.sample_names,
            selected_markers=self.selected_markers,
            genetic_maps_unfiltered=self.genetic_maps_unfiltered
        ))
