# =============================================================================
# file_selection_frame.py
#
# CustomTkinter Frame for selecting input files (VCF + Trait) for QTL analysis.
# =============================================================================

from customtkinter import *
from tkinter import filedialog
from PIL import Image
import os


class FileSelectionFrame(CTkFrame):
    """
    GUI frame for selecting input files required to start the QTL Mapping process.
    
    Files Required:
    - VCF File (.vcf)
    - Trait File (.txt)

    Provides feedback and disables/enables the 'Next' button based on selection status.
    """

    def __init__(self, master, on_back, on_next, styles):
        """
        Initialize the frame.

        Parameters:
            master (Tk): Parent window.
            on_back (function): Callback for 'Back' button.
            on_next (function): Callback for 'Next' button.
            styles (dict): Dictionary with button style themes.
        """
        super().__init__(master, width=1100, height=650)
        self.master = master
        self.on_back = on_back
        self.on_next = on_next
        self.styles = styles

        self.vcf_file_path = None
        self.trait_file_path = None

        self.pack(fill="both", expand=True)
        self.build_ui()

    def build_ui(self):
        """
        Construct and layout the UI components.
        """
        # === Center the main window
        width, height = 1100, 650
        self.master.geometry(f"{width}x{height}")
        self.master.update_idletasks()
        screen_width = self.master.winfo_screenwidth()
        screen_height = self.master.winfo_screenheight()
        x = (screen_width // 2) - (width // 2)
        y = (screen_height // 2) - (height // 2) - 100
        self.master.geometry(f"{width}x{height}+{x}+{y}")

        # === Side image (optional)
        try:
            step1_image = CTkImage(light_image=Image.open("gui_frames/Step1.png"), size=(250, 650))
            CTkLabel(self, image=step1_image, text="").place(x=0, y=0)
            self.step1_image = step1_image  # retain reference
        except Exception as e:
            print(f"Step1.png not found: {e}")

        # === Main content area
        content_frame = CTkFrame(self, fg_color="transparent")
        content_frame.place(relx=0.62, rely=0.5, anchor="center")

        CTkLabel(
            content_frame, text="Select the VCF and Trait files", font=("Segoe UI", 22, "bold")
        ).pack(pady=20)

        # === VCF file selection
        CTkLabel(
            content_frame, text="Please select the VCF file from file directory:", font=("Segoe UI", 14)
        ).pack(pady=(10, 5))

        CTkButton(
            content_frame, text="Choose File", font=("Segoe UI", 14, "bold"),
            command=self.select_vcf_file, **self.styles["red"]
        ).pack()

        self.vcf_label = CTkLabel(content_frame, text="", font=("Segoe UI", 12))
        self.vcf_label.pack(pady=(5, 2))

        self.vcf_success_label = CTkLabel(
            content_frame, text="", font=("Segoe UI", 16, "bold"), text_color="#008000"
        )
        self.vcf_success_label.pack(pady=(0, 8))

        # === Trait file selection
        CTkLabel(
            content_frame, text="Please select the trait file from file directory:", font=("Segoe UI", 14)
        ).pack(pady=(10, 5))

        CTkButton(
            content_frame, text="Choose File", font=("Segoe UI", 14, "bold"),
            command=self.select_trait_file, **self.styles["red"]
        ).pack()

        self.trait_label = CTkLabel(content_frame, text="", font=("Segoe UI", 12))
        self.trait_label.pack(pady=(5, 2))

        self.trait_success_label = CTkLabel(
            content_frame, text="", font=("Segoe UI", 16, "bold"), text_color="#008000"
        )
        self.trait_success_label.pack(pady=(0, 8))

        # === Message area for errors
        self.message_label = CTkLabel(content_frame, text="", font=("Segoe UI", 12))
        self.message_label.pack(pady=(10, 10))

        # === Navigation buttons
        nav_frame = CTkFrame(content_frame, fg_color="transparent")
        nav_frame.pack(pady=20, fill="x")

        self.back_btn = CTkButton(
            nav_frame, text="Back", font=("Segoe UI", 14, "bold"),
            width=120, command=self.on_back, **self.styles["white"]
        )
        self.back_btn.pack(side="left", padx=(0, 20))

        self.next_btn = CTkButton(
            nav_frame, text="Run Full QTL Analysis", font=("Segoe UI", 16, "bold"),
            width=180, command=self.try_next
        )
        self.next_btn.pack(side="left")

        self.update_next_button_state()

    # === File Selectors ===

    def select_vcf_file(self):
        file_path = filedialog.askopenfilename(filetypes=[("VCF files", "*.vcf")])
        if file_path and file_path.endswith(".vcf"):
            self.vcf_file_path = file_path
            self.vcf_label.configure(text=f"File Selected: {os.path.basename(file_path)}")
            self.vcf_success_label.configure(text="File uploaded successfully")
            self.message_label.configure(text="")
        else:
            self.vcf_file_path = None
            self.vcf_label.configure(text="")
            self.vcf_success_label.configure(text="")
            if file_path:
                self.message_label.configure(text="Error: Please select a valid .vcf file", text_color="red")
        self.update_next_button_state()

    def select_trait_file(self):
        file_path = filedialog.askopenfilename(filetypes=[("Text files", "*.txt")])
        if file_path and file_path.endswith(".txt"):
            self.trait_file_path = file_path
            self.trait_label.configure(text=f"File Selected: {os.path.basename(file_path)}")
            self.trait_success_label.configure(text="File uploaded successfully")
            self.message_label.configure(text="")
        else:
            self.trait_file_path = None
            self.trait_label.configure(text="")
            self.trait_success_label.configure(text="")
            if file_path:
                self.message_label.configure(text="Error: Please select a valid .txt file", text_color="red")
        self.update_next_button_state()

    def update_next_button_state(self):
        """
        Enable 'Next' only if both files are selected.
        """
        if self.vcf_file_path and self.trait_file_path:
            self.next_btn.configure(
                state="normal", fg_color="#008000", hover_color="#006400", text_color="#ffffff"
            )
        else:
            self.next_btn.configure(
                state="disabled", fg_color="#cccccc", hover_color="#bbbbbb", text_color="#222222"
            )

    def try_next(self):
        """
        Trigger callback if both files are selected.
        """
        if self.vcf_file_path and self.trait_file_path:
            self.on_next(self.vcf_file_path, self.trait_file_path)
