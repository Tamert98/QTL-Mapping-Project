# styles.py

red_button_style = {
    "fg_color": "#e50914",
    "hover_color": "#b00610",
    "text_color": "#ffffff"
}

white_button_style = {
    "fg_color": "#ffffff",
    "hover_color": "#e6e6e6",
    "text_color": "#000000"
}

def get_styles(mode):
    if mode == "dark":
        return {"red": red_button_style, "white": white_button_style}
    else:
        return {
            "red": red_button_style,
            "white": white_button_style
        }
