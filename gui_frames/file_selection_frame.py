from customtkinter import *
from tkinter import filedialog
import os

class FileSelectionFrame(CTkFrame):
    def __init__(self, master, on_back, on_next, styles):
        super().__init__(master)
        self.master = master
        self.on_back = on_back  # Function to call when back is pressed
        self.on_next = on_next  # Function to call when next is pressed
        self.styles = styles

        self.vcf_file_path = None
        self.trait_file_path = None

        self.build_ui()

    def build_ui(self):
        # Resize and center window when this frame is shown
        width, height = 700, 500
        self.master.geometry(f"{width}x{height}")
        self.master.update_idletasks()
        screen_width = self.master.winfo_screenwidth()
        screen_height = self.master.winfo_screenheight()
        x = (screen_width // 2) - (width // 2)
        y = (screen_height // 2) - (height // 2)
        self.master.geometry(f"{width}x{height}+{x}+{y}")
        self.pack(fill="both", expand=True)

        CTkLabel(self, text="Select the VCF and Trait files", font=("Segoe UI", 22, "bold")).pack(pady=20)

        # VCF selection
        CTkLabel(self, text="Please select the VCF file from file directory:", font=("Segoe UI", 14)).pack(pady=(10, 5))
        CTkButton(self, text="Choose File", font=("Segoe UI", 14, "bold"),
                  command=self.select_vcf_file, **self.styles["red"]).pack()
        self.vcf_label = CTkLabel(self, text="", font=("Segoe UI", 12))
        self.vcf_label.pack(pady=(5, 2))
        self.vcf_success_label = CTkLabel(self, text="", font=("Segoe UI", 12),
                                          text_color=self.styles["red"]["fg_color"])
        self.vcf_success_label.pack(pady=(0, 8))

        # Trait selection
        CTkLabel(self, text="Please select the trait file from file directory:", font=("Segoe UI", 14)).pack(pady=(10, 5))
        CTkButton(self, text="Choose File", font=("Segoe UI", 14, "bold"),
                  command=self.select_trait_file, **self.styles["red"]).pack()
        self.trait_label = CTkLabel(self, text="", font=("Segoe UI", 12))
        self.trait_label.pack(pady=(5, 2))
        self.trait_success_label = CTkLabel(self, text="", font=("Segoe UI", 12),
                                            text_color=self.styles["red"]["fg_color"])
        self.trait_success_label.pack(pady=(0, 8))

        # Message label
        self.message_label = CTkLabel(self, text="", font=("Segoe UI", 12))
        self.message_label.pack(pady=(10, 10))

        # Navigation buttons
        nav_frame = CTkFrame(self, fg_color="transparent")
        nav_frame.pack(pady=20, fill="x", padx=20)

        self.back_btn = CTkButton(nav_frame, text="Back", font=("Segoe UI", 14, "bold"),
                                  width=120, command=self.on_back, **self.styles["white"])
        self.next_btn = CTkButton(nav_frame, text="Next", font=("Segoe UI", 14, "bold"),
                                  width=120, command=self.try_next)
        self.update_next_button_state()

        self.back_btn.pack(side="left", padx=5)
        self.next_btn.pack(side="right", padx=5)

    def select_vcf_file(self):
        file_path = filedialog.askopenfilename(filetypes=[("VCF files", "*.vcf")])
        if file_path:
            if not file_path.endswith(".vcf"):
                self.message_label.configure(text="Error: Please select a valid .vcf file", text_color="red")
                self.vcf_file_path = None
                self.vcf_label.configure(text="")
                self.vcf_success_label.configure(text="")
            else:
                self.vcf_file_path = file_path
                self.vcf_label.configure(text=os.path.basename(file_path))
                self.message_label.configure(text="")
                self.vcf_success_label.configure(text="File uploaded successfully")
        else:
            self.vcf_success_label.configure(text="")
        self.update_next_button_state()

    def select_trait_file(self):
        file_path = filedialog.askopenfilename(filetypes=[("Text files", "*.txt")])
        if file_path:
            if not file_path.endswith(".txt"):
                self.message_label.configure(text="Error: Please select a valid .txt file", text_color="red")
                self.trait_file_path = None
                self.trait_label.configure(text="")
                self.trait_success_label.configure(text="")
            else:
                self.trait_file_path = file_path
                self.trait_label.configure(text=os.path.basename(file_path))
                self.message_label.configure(text="")
                self.trait_success_label.configure(text="File uploaded successfully")
        else:
            self.trait_success_label.configure(text="")
        self.update_next_button_state()

    def update_next_button_state(self):
        if self.vcf_file_path and self.trait_file_path:
            self.next_btn.configure(state="normal", **self.styles["red"])
        else:
            dark_red = {"fg_color": "#700000", "hover_color": "#550000", "text_color": "#ffffff"}
            self.next_btn.configure(state="disabled", **dark_red)

    def try_next(self):
        if self.vcf_file_path and self.trait_file_path:
            self.on_next(self.vcf_file_path, self.trait_file_path)
