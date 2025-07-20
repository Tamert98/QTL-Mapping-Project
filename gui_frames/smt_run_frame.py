from customtkinter import *
from PIL import Image
import threading
from single_marker_test import run_smt_for_all_traits
from plot_utils import generate_concatenated_qtl_pdf

class SMTRunFrame(CTkFrame):
    def __init__(self, master, vcf_data, traits, sample_names, on_back, styles):
        super().__init__(master, width=1100, height=650)
        self.master = master
        self.vcf_data = vcf_data
        self.traits = traits
        self.sample_names = sample_names
        self.styles = styles
        self.on_back = on_back
        self.results = None

        self.pack(fill="both", expand=True)
        self.message_lines = []
        self.textbox = None
        self.build_ui()
        threading.Thread(target=self.run_smt_thread, daemon=True).start()

    def build_ui(self):
        CTkLabel(self, text="Applying SMT", font=("Segoe UI", 22, "bold")).pack(pady=(30, 10))
        self.status_label = CTkLabel(self, text="The SMT algorithm is being applied...", font=("Segoe UI", 16))
        self.status_label.pack(pady=10)

        # Textbox for progress messages
        self.textbox = CTkTextbox(self, width=800, height=220, font=("Segoe UI", 13), wrap="word")
        self.textbox.pack(pady=20)
        self.textbox.insert("end", "SMT progress will appear here...\n")
        self.textbox.configure(state="disabled")

        green_color = "#228B22"
        self.plots_label = CTkLabel(self, text="✔ Plots saved in Results/SMT/Plots",
                                    font=("Segoe UI", 14), text_color=green_color)
        self.reports_label = CTkLabel(self, text="✔ Reports saved in Results/SMT/Reports",
                                      font=("Segoe UI", 14), text_color=green_color)
        self.plots_label.pack_forget()
        self.reports_label.pack_forget()

        self.back_button = CTkButton(self, text="Back to Menu", font=("Segoe UI", 14, "bold"),
                                     command=self.finish_and_return,
                                     width=160, height=44, state="disabled", **self.styles["red"])
        self.back_button.pack(pady=30)

    def run_smt_thread(self):
        def gui_message_callback(msg):
            self.message_lines.append(msg)
            def update_textbox():
                self.textbox.configure(state="normal")
                self.textbox.insert("end", msg + "\n")
                self.textbox.see("end")
                self.textbox.configure(state="disabled")
            self.master.after(0, update_textbox)

        self.textbox.configure(state="normal")
        self.textbox.delete("1.0", "end")
        self.textbox.insert("end", "SMT progress will appear here...\n")
        self.textbox.configure(state="disabled")

        self.results = run_smt_for_all_traits(self.vcf_data, self.traits, self.sample_names, message_callback=gui_message_callback)
        # Update status label to indicate plot generation
        self.master.after(0, lambda: self.status_label.configure(text="Generating plots for all traits... Please wait.", text_color=self.styles["red"]["fg_color"]))
        generate_concatenated_qtl_pdf(self.results, message_callback=gui_message_callback)
        self.master.after(2000, self.on_smt_done)

    def on_smt_done(self):
        self.status_label.configure(text="✔ SMT Completed", text_color="#228B22")
        self.plots_label.pack()
        self.reports_label.pack()
        self.back_button.configure(state="normal")

    def finish_and_return(self):
        self.on_back(self.results)
