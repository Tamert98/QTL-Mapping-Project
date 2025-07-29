# =============================================================================
# format_window.py
#
# CustomTkinter Format Requirements Viewer
# Displays instructions for formatting VCF and trait files.
# =============================================================================

from customtkinter import *

def open_format_window(current_style, parent=None):
    """
    Opens a modal window showing the VCF and Trait file format requirements.

    Parameters:
        current_style (dict): Dictionary with color styles for buttons.
                              Example:
                              {
                                  "red": { ... },   # Red close button
                                  "white": { ... }  # Navigation arrows
                              }
        parent (CTk | None): Optional parent window for modal behavior.
    """
    # === Extract button styles
    red_button_style = current_style["red"]
    white_button_style = current_style["white"]

    # === Format instructions for each file type
    pages = [
        (
            "VCF File Requirements",
            (
                "The VCF file contains genetic markers and genotypes. It must follow this format:\n\n"
                "1. Chromosome and Contig Info:\n"
                "Each contig must be defined like this:\n"
                "  ##contig=<ID=chr01,length=100000>\n\n"
                "2. Header Line:\n"
                "A single line must start with:\n"
                "  #CHROM  POS  ID  REF  ALT  QUAL  ...  FORMAT  sample1  sample2 ...\n\n"
                "3. Sample Names:\n"
                "Sample names appear after FORMAT.\n"
                "We extract them by removing extensions:\n"
                "  \"male-1.bowtie\" → \"male-1\"\n\n"
                "4. FORMAT and Genotype Fields:\n"
                "Each sample has a string like:\n"
                "  GT:PL:DP:GQ\n\n"
                "Example sample entry:\n"
                "  0/1:35,0,120:20:99\n\n"
                "Common genotype codes:\n"
                "  0/0 → Homozygous reference\n"
                "  0/1 → Heterozygous\n"
                "  1/1 → Homozygous alternate\n"
                "  ././ → Missing\n\n"
                "5. Parsed Output:\n"
                "Markers are grouped by chromosome.\n"
                "Each marker includes:\n"
                "  - POS, REF, ALT, QUAL\n"
                "  - Per-sample genotype data"
            )
        ),
        (
            "Trait File Requirements",
            (
                "The trait file contains phenotype values for each sample.\n\n"
                "1. First Line – Sample Names:\n"
                "A single line like:\n"
                "  male-1.bowtie  male-2.bowtie  male-3.bowtie\n\n"
                "Sample names are trimmed:\n"
                "  \"male-1.bowtie\" → \"male-1\"\n\n"
                "2. Trait Lines – Values per Trait:\n"
                "Each following line has:\n"
                "  - Trait name\n"
                "  - Value per sample\n"
                "  - Use \"$\" for missing data\n\n"
                "Example:\n"
                "  Height   0.52   $   0.64\n\n"
                "This maps to:\n"
                "  \"Height\": {\n"
                "    \"male-1\": 0.52,\n"
                "    \"male-2\": None,\n"
                "    \"male-3\": 0.64\n"
                "  }\n\n"
                "3. Format Rules:\n"
                "- Number of values must match number of samples\n"
                "- Values must be numeric or \"$\"\n"
                "- Lines must be space or tab-separated"
            )
        )
    ]

    # === Create the main format window
    format_win = CTkToplevel()
    format_win.title("File Format Requirements")
    format_win.geometry("900x550")

    # Center on screen
    format_win.update_idletasks()
    screen_width = format_win.winfo_screenwidth()
    screen_height = format_win.winfo_screenheight()
    x = int((screen_width - 900) / 2)
    y = int((screen_height - 550) / 2)
    format_win.geometry(f"900x550+{x}+{y}")

    # Make modal if needed
    if parent:
        format_win.transient(parent)
        format_win.grab_set()
        format_win.focus()

    # === Page Frame (title + content)
    page_frame = CTkFrame(format_win)
    page_frame.pack(padx=20, pady=(10, 0), expand=True, fill="both")

    title_label = CTkLabel(
        page_frame, text="", font=("Segoe UI", 18, "bold"), anchor="center"
    )
    title_label.pack(pady=(10, 5))

    text_area = CTkTextbox(
        page_frame, wrap="word", font=("Segoe UI", 14), state="disabled"
    )
    text_area.pack(expand=True, fill="both", padx=10)

    # === Navigation Controls
    nav_frame = CTkFrame(format_win, fg_color="transparent")
    nav_frame.pack(pady=10, fill="x", padx=20)

    def show_page(index):
        """Update content of the viewer to the selected page."""
        title, content = pages[index]
        title_label.configure(text=title)
        text_area.configure(state="normal")
        text_area.delete("1.0", "end")
        text_area.insert("1.0", content)
        text_area.configure(state="disabled")

    # Buttons: ←, →, Close
    back_btn = CTkButton(
        nav_frame, text="←", width=50, font=("Segoe UI", 16, "bold"),
        command=lambda: show_page(0), **white_button_style
    )
    next_btn = CTkButton(
        nav_frame, text="→", width=50, font=("Segoe UI", 16, "bold"),
        command=lambda: show_page(1), **white_button_style
    )
    close_btn = CTkButton(
        nav_frame, text="Close", font=("Segoe UI", 14, "bold"), width=120,
        command=format_win.destroy, **red_button_style
    )

    back_btn.pack(side="left", padx=5)
    next_btn.pack(side="left", padx=5)
    close_btn.pack(side="right", padx=5)

    # Initialize first page
    show_page(0)
