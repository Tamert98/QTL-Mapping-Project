import os
from customtkinter import *
from PIL import Image
from gui_frames.genetic_data_viewer_frame import GeneticMapViewerFrame
from gui_frames.smt_results_frame import SMTConcatenatedResultsFrame
from gui_frames.smt_results_frame import SMTPerChromosomeResultsFrame
from gui_frames.sim_view_results_frame import SIMViewResultsFrame
from gui_frames.smt_sim_compare_frame import SMT_SIM_CompareFrame

class ViewResultsFrame(CTkFrame):
    def __init__(
        self,
        master,
        sim_results,
        smt_results,
        vcf_data,
        traits,
        sample_names,
        selected_markers,
        genetic_maps_unfiltered,
        go_back_callback=None
    ):
        super().__init__(master, fg_color="#1e1e1e")

        self.sim_results = sim_results
        self.smt_results = smt_results
        self.vcf_data = vcf_data
        self.traits = traits
        self.sample_names = sample_names
        self.selected_markers = selected_markers
        self.genetic_maps_unfiltered = genetic_maps_unfiltered
        self.go_back_callback = go_back_callback

        self.pack(fill="both", expand=True)
        self.build_ui()

    def build_ui(self):
        width, height = 1100, 650
        self.master.geometry(f"{width}x{height}")
        self.master.update_idletasks()
        screen_width = self.master.winfo_screenwidth()
        screen_height = self.master.winfo_screenheight()
        x = (screen_width // 2) - (width // 2)
        y = (screen_height // 2) - (height // 2) - 100
        self.master.geometry(f"{width}x{height}+{x}+{y}")

        try:
            step4_image = CTkImage(light_image=Image.open("gui_frames/Step4.png"), size=(250, 650))
            CTkLabel(self, image=step4_image, text="").place(x=0, y=0)
            self.step4_image = step4_image
        except Exception as e:
            print(f"Step4.png not found: {e}")

        content = CTkFrame(self, fg_color="#1e1e1e")
        content.place(relx=0.62, rely=0.5, anchor="center")
        content.grid_columnconfigure(0, weight=1)

        CTkLabel(content, text="View Analysis Results", font=("Arial", 24, "bold"), text_color="white").pack(pady=(20, 30))

        # === Genotype Data ===
        geno_box = CTkFrame(content, fg_color="#333333")
        geno_box.pack(pady=10, fill="x")
        CTkLabel(geno_box, text="View Genotype Data Analysis", font=("Arial", 18, "bold"), text_color="white").pack(pady=(10, 5))
        self.geno_dropdown = CTkOptionMenu(
            geno_box,
            values=[
                "Chromosome-wide genetic maps",
                "Filtered Markers Distance comparison",
                "Unfiltered Markers Distance comparison",
                "HeatMaps"
            ],
            fg_color="#00ff00",
            button_color="#00ff00",
            text_color="white",
            dropdown_text_color="white",
            font=("Arial", 15, "bold"),
            dropdown_font=("Arial", 15),
            width=520,
            command=self.launch_genetic_map_viewer
        )
        self.geno_dropdown.pack(pady=(0, 10), padx=20)

        # === SMT ===
        smt_box = CTkFrame(content, fg_color="#333333")
        smt_box.pack(pady=10, fill="x")
        CTkLabel(smt_box, text="View SMT Results", font=("Arial", 18, "bold"), text_color="white").pack(pady=(10, 5))
        self.smt_dropdown = CTkOptionMenu(
            smt_box,
            values=[
                "Overall pvalues across the chromosomes and locate best QTL location",
                "p-value plots per chromosome"
            ],
            fg_color="#00ff00",
            button_color="#00ff00",
            text_color="white",
            dropdown_text_color="white",
            font=("Arial", 15, "bold"),
            dropdown_font=("Arial", 15),
            width=520,
            command=self.launch_smt_viewer
        )
        self.smt_dropdown.pack(pady=(0, 10), padx=20)

        # === SIM ===
        sim_box = CTkFrame(content, fg_color="#333333")
        sim_box.pack(pady=10, fill="x")
        CTkLabel(sim_box, text="View SIM Results", font=("Arial", 18, "bold"), text_color="white").pack(pady=(10, 5))
        CTkButton(
            sim_box,
            text="Browse LOD score curves for each chromosome and trait",
            font=("Arial", 15, "bold"),
            fg_color="#00ff00",
            hover_color="#00cc00",
            text_color="white",
            height=44,
            width=520,
            command=self.launch_sim_viewer
        ).pack(pady=(5, 10), padx=20)

        # === Comparison ===
        compare_box = CTkFrame(content, fg_color="#333333")
        compare_box.pack(pady=10, fill="x")
        CTkLabel(compare_box, text="Compare SMT vs SIM Methods", font=("Arial", 18, "bold"), text_color="white").pack(pady=(10, 5))
        CTkButton(
            compare_box,
            text="Visual and statistical comparison of SMT vs SIM outputs",
            font=("Arial", 15, "bold"),
            fg_color="#00ff00",
            hover_color="#00cc00",
            text_color="white",
            height=44,
            width=520,
            command=self.launch_smt_sim_compare  # ✅ FIXED
        ).pack(pady=(5, 10), padx=20)

        # === Back Button ===
        CTkButton(
            content,
            text="Back to Start Page",
            command=self.go_back_callback,
            font=("Arial", 14, "bold"),
            fg_color="white",
            text_color="black",
            hover_color="#dddddd",
            height=40,
            width=180
        ).pack(pady=(30, 20))

    def launch_genetic_map_viewer(self, selection):
        selection_map = {
            "Chromosome-wide genetic maps": "genetic",
            "Filtered Markers Distance comparison": "compare-filtered",
            "Unfiltered Markers Distance comparison": "compare",
            "HeatMaps": "heatmap"
        }
        mode = selection_map.get(selection)
        if not mode:
            return

        self.destroy()
        GeneticMapViewerFrame(
            master=self.master,
            vcf_data=self.vcf_data,
            styles={"white": {"fg_color": "white", "text_color": "black", "hover_color": "#dddddd"}},
            on_back=self.rebuild_self,
            mode=mode
        )

    def launch_smt_viewer(self, selection):
        self.destroy()
        if selection == "Overall pvalues across the chromosomes and locate best QTL location":
            SMTConcatenatedResultsFrame(
                master=self.master,
                traits=self.traits,
                all_smt_results=self.smt_results,
                genetic_maps_unfiltered=self.genetic_maps_unfiltered,
                styles={"white": {"fg_color": "white", "text_color": "black", "hover_color": "#dddddd"}},
                on_back=self.rebuild_self
            )
        elif selection == "p-value plots per chromosome":
            SMTPerChromosomeResultsFrame(
                master=self.master,
                traits=self.traits,
                styles={"white": {"fg_color": "white", "text_color": "black", "hover_color": "#dddddd"}},
                on_back=self.rebuild_self
            )

    def launch_sim_viewer(self):
        self.destroy()
        SIMViewResultsFrame(
            master=self.master,
            traits=self.traits,
            styles={"white": {"fg_color": "white", "text_color": "black", "hover_color": "#dddddd"}},
            on_back=self.rebuild_self
        )

    def launch_smt_sim_compare(self):  # ✅ FIXED
        self.destroy()
        SMT_SIM_CompareFrame(
            master=self.master,
            traits=self.traits,
            on_back=self.rebuild_self
        )

    def rebuild_self(self):
        for widget in self.master.winfo_children():
            widget.destroy()
        ViewResultsFrame(
            master=self.master,
            sim_results=self.sim_results,
            smt_results=self.smt_results,
            vcf_data=self.vcf_data,
            traits=self.traits,
            sample_names=self.sample_names,
            selected_markers=self.selected_markers,
            genetic_maps_unfiltered=self.genetic_maps_unfiltered,
            go_back_callback=self.go_back_callback
        )
