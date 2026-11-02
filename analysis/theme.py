"""Display labels and colours shared by the static figures and the dashboard.

Colours come from a validated, colour-blind-safe categorical palette: responders
always take slot 1 (blue) and non-responders slot 2 (orange), in every chart and
in both light and dark mode, so a colour always means the same group.
"""

POPULATION_LABELS = {
    "b_cell": "B cells",
    "cd8_t_cell": "CD8 T cells",
    "cd4_t_cell": "CD4 T cells",
    "nk_cell": "NK cells",
    "monocyte": "Monocytes",
}
RESPONSE_LABELS = {"yes": "Responders", "no": "Non-responders"}
SEX_LABELS = {"M": "Male", "F": "Female"}

LIGHT = {
    "surface": "#fcfcfb",
    "ink": "#0b0b0b",
    "ink_secondary": "#52514e",
    "ink_muted": "#898781",
    "grid": "#e1e0d9",
    "axis": "#c3c2b7",
    "responders": "#2a78d6",
    "non_responders": "#eb6834",
    "highlight": "#6250d6",  # violet: single highlighted series, kept apart from the response colours
}
DARK = {
    "surface": "#1a1a19",
    "ink": "#f0efec",
    "ink_secondary": "#c3c2b7",
    "ink_muted": "#898781",
    "grid": "#2c2c2a",
    "axis": "#383835",
    "responders": "#3987e5",
    "non_responders": "#d95926",
    "highlight": "#9085e9",
}


def population_label(code: str) -> str:
    return POPULATION_LABELS.get(code, code)
