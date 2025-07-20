from customtkinter import *
from PIL import Image

class MainMenuFrame(CTkFrame):
    def __init__(self, master, on_back, on_view_maps, on_apply_smt, on_view_smt_results, styles, smt_ready=False):
        super().__init__(master, width=1100, height=650)
        self.master = master
        self.on_back = on_back
        self.on_view_maps = on_view_maps
        self.on_apply_smt = on_apply_smt
        self.on_view_smt_results = on_view_smt_results
        self.styles = styles
        self.smt_ready = smt_ready

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
            side_image = CTkImage(light_image=Image.open("sideimg.png"), size=(250, 650))
            CTkLabel(self, image=side_image, text="").place(x=0, y=0)
            self.side_image = side_image
        except:
            print("Side image not found.")

        content = CTkFrame(self, fg_color="transparent")
        content.place(relx=0.62, rely=0.4, anchor="center")

        CTkLabel(content, text="Analysis Menu", font=("Segoe UI", 22, "bold")).pack(pady=(10, 30))

        CTkButton(content, text="View Maps", font=("Segoe UI", 16, "bold"),
                  width=240, height=45, command=self.on_view_maps, **self.styles["red"]).pack(pady=(0, 45))

        smt_frame = CTkFrame(content, fg_color="transparent")
        smt_frame.pack(pady=20)
        CTkButton(smt_frame, text="Apply SMT on the data", font=("Segoe UI", 14, "bold"),
                  width=200, height=40, command=self.on_apply_smt, **self.styles["white"]).pack(side="left", padx=(0, 20))

        CTkButton(smt_frame, text="View Results", font=("Segoe UI", 14, "bold"),
                  width=140, height=40,
                  state="normal" if self.smt_ready else "disabled",
                  command=self.on_view_smt_results if self.smt_ready else None,
                  **(self.styles["red"] if self.smt_ready else self.get_dark_red())).pack(side="left")

        sim_frame = CTkFrame(content, fg_color="transparent")
        sim_frame.pack(pady=20)
        CTkButton(sim_frame, text="Apply SIM on the data", font=("Segoe UI", 14, "bold"),
                  width=200, height=40, command=lambda: None, **self.styles["white"]).pack(side="left", padx=(0, 20))
        CTkButton(sim_frame, text="View Results", font=("Segoe UI", 14, "bold"),
                  width=140, height=40, state="disabled", **self.get_dark_red()).pack(side="left")

        CTkButton(content, text="Back", font=("Segoe UI", 14, "bold"),
                  width=140, height=40, command=self.on_back, **self.styles["white"]).pack(pady=(30, 10))

    def get_dark_red(self):
        return {
            "fg_color": "#700000",
            "hover_color": "#550000",
            "text_color": "#ffffff"
        }
