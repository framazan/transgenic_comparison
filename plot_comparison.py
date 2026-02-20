"""
Produces two figures:
  1. gffcompare_comparison.png — Precision vs Recall scatter plots
     at Exon level (top row) and Gene/Locus level (bottom row).
  2. busco_comparison.png — BUSCO completeness stacked bar chart.
"""

import os
import re
import json
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


INCLUDE_GENOMES = [
    'A_thaliana', 'O_sativa', 'Z_mays', 'P_trichocarpa'
]


BASE_DIR = os.path.dirname(os.path.abspath(__file__))
GFFCOMPARE_DIR = os.path.join(BASE_DIR, "gffcompare_results")
BUSCO_DIR = os.path.join(BASE_DIR, "busco_results")
GFF_OUTPUT = os.path.join(BASE_DIR, "gffcompare_comparison.png")
BUSCO_OUTPUT = os.path.join(BASE_DIR, "busco_comparison.png")

TOOL_STYLES = {
    "annevo": {"color": "#2ca02c", "marker": "o", "label": "Annevo"},
    "helixer": {"color": "#d62728", "marker": "^", "label": "Helixer"},
    "tiberius": {"color": "#1f77b4", "marker": "v", "label": "Tiberius"},
    "tiberius_softmasked": {
        "color": "#7f7f7f",
        "marker": ">",
        "label": "Tiberius (softmasked)",
    },
}
TOOLS_ORDER = ["annevo", "helixer", "tiberius", "tiberius_softmasked"]

GFF_LEVELS = [
    ("Exon level", "Exon level"),
    ("Locus level", "Gene level"),
]

def parse_name(name: str):
    # Split 'A_thaliana_annevo' → ('A_thaliana', 'annevo').
    # The tool suffix is matched longest‑first so that
    # 'tiberius_softmasked' is not mistaken for 'tiberius'.
    for tool in ["tiberius_softmasked", "tiberius", "helixer", "annevo"]:
        suffix = "_" + tool
        if name.endswith(suffix):
            return name[: -len(suffix)], tool
    return None, None


def parse_gffcompare_stats(filepath: str) -> dict:
    # Return {level: {sensitivity, precision}} from a .stats file.
    results = {}
    with open(filepath) as fh:
        for line in fh:
            line = line.strip()
            for level in [
                "Base level",
                "Exon level",
                "Intron level",
                "Intron chain level",
                "Transcript level",
                "Locus level",
            ]:
                if line.startswith(level):
                    nums = re.findall(r"[\d.]+", line)
                    if len(nums) >= 2:
                        results[level] = {
                            "sensitivity": float(nums[0]),
                            "precision": float(nums[1]),
                        }
    return results


def parse_busco_json(filepath: str) -> dict:
    # Return BUSCO completeness percentages from a summary JSON.
    with open(filepath) as fh:
        data = json.load(fh)
    r = data.get("results", {})
    return {
        "complete": r.get("Complete percentage", 0),
        "single_copy": r.get("Single copy percentage", 0),
        "multi_copy": r.get("Multi copy percentage", 0),
        "fragmented": r.get("Fragmented percentage", 0),
        "missing": r.get("Missing percentage", 0),
    }


def format_species(key: str) -> str:
    # 'A_thaliana' → 'A. thaliana'
    parts = key.split("_")
    return f"{parts[0]}. {parts[1]}" if len(parts) == 2 else key.replace("_", " ")




def collect_data():
    gff_data = {}
    busco_data = {}
    species_set = set()

    # gffcompare
    for path in sorted(glob.glob(os.path.join(GFFCOMPARE_DIR, "*.stats"))):
        name = os.path.splitext(os.path.basename(path))[0]
        sp, tool = parse_name(name)
        if sp and tool:
            species_set.add(sp)
            gff_data[(sp, tool)] = parse_gffcompare_stats(path)

    # BUSCO
    for entry in sorted(os.listdir(BUSCO_DIR)):
        entry_path = os.path.join(BUSCO_DIR, entry)
        if os.path.isdir(entry_path):
            jsons = glob.glob(os.path.join(entry_path, "*.json"))
            if jsons:
                sp, tool = parse_name(entry)
                if sp and tool:
                    species_set.add(sp)
                    busco_data[(sp, tool)] = parse_busco_json(jsons[0])

    # Filter species_set if INCLUDE_GENOMES is set
    if INCLUDE_GENOMES:
        filtered = [sp for sp in sorted(species_set) if sp in INCLUDE_GENOMES]
    else:
        filtered = sorted(species_set)
    return gff_data, busco_data, filtered




def plot_gffcompare(gff_data, species_list):
    # Scatter plots: Precision vs Recall at two gffcompare levels.
    n_sp = len(species_list)
    n_rows = len(GFF_LEVELS)
    col_w, row_h = 2.6, 2.8

    fig, axes = plt.subplots(
        n_rows,
        n_sp,
        figsize=(col_w * n_sp, row_h * n_rows + 1.0),
        squeeze=False,
    )

    for ri, (gff_key, row_label) in enumerate(GFF_LEVELS):
        for ci, sp in enumerate(species_list):
            ax = axes[ri, ci]

            for tool in TOOLS_ORDER:
                key = (sp, tool)
                if key in gff_data and gff_key in gff_data[key]:
                    d = gff_data[key][gff_key]
                    s = TOOL_STYLES[tool]
                    ax.scatter(
                        d["precision"],
                        d["sensitivity"],
                        color=s["color"],
                        marker=s["marker"],
                        s=120,
                        edgecolors="black",
                        linewidth=0.5,
                        zorder=5,
                    )

            ax.set_xlim(0, 100)
            ax.set_ylim(0, 100)
            ax.grid(True, alpha=0.2, linestyle="--")
            ax.tick_params(labelsize=7)

            # Column header (species) on top row only
            if ri == 0:
                ax.set_title(
                    format_species(sp), fontsize=9, fontstyle="italic"
                )

            # x‑label on bottom row only
            if ri == n_rows - 1:
                ax.set_xlabel("Precision [%]", fontsize=8)
            else:
                ax.set_xticklabels([])

            # y‑label on first column only
            if ci == 0:
                ax.set_ylabel("Recall [%]", fontsize=9)

        # Row label on the far left
        axes[ri, 0].annotate(
            row_label,
            xy=(-0.45, 0.5),
            xycoords="axes fraction",
            fontsize=11,
            fontweight="bold",
            rotation=90,
            va="center",
            ha="center",
        )

    handles = [
        plt.Line2D(
            [],
            [],
            color=TOOL_STYLES[t]["color"],
            marker=TOOL_STYLES[t]["marker"],
            linestyle="None",
            markersize=10,
            markeredgecolor="black",
            markeredgewidth=0.5,
            label=TOOL_STYLES[t]["label"],
        )
        for t in TOOLS_ORDER
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        ncol=len(TOOLS_ORDER),
        fontsize=9,
        frameon=True,
        fancybox=True,
        edgecolor="gray",
    )

    plt.tight_layout(rect=[0.05, 0.06, 1, 1])
    fig.savefig(GFF_OUTPUT, dpi=300, bbox_inches="tight")
    print(f"Saved gffcompare figure → {GFF_OUTPUT}")
    plt.close(fig)


def plot_busco(busco_data, species_list):
    # Stacked bar chart of BUSCO completeness per tool and genome.
    n_sp = len(species_list)
    col_w = 2.6

    fig, axes = plt.subplots(
        1, n_sp, figsize=(col_w * n_sp, 4), squeeze=False
    )

    for ci, sp in enumerate(species_list):
        ax = axes[0, ci]
        tools_present = [t for t in TOOLS_ORDER if (sp, t) in busco_data]

        if not tools_present:
            ax.set_visible(False)
            continue

        x_pos = np.arange(len(tools_present))

        for i, tool in enumerate(tools_present):
            bd = busco_data[(sp, tool)]
            s = TOOL_STYLES[tool]

            # Complete (solid)
            ax.bar(
                i,
                bd["complete"],
                color=s["color"],
                edgecolor="black",
                linewidth=0.5,
            )
            # Fragmented (hatched, stacked on top)
            ax.bar(
                i,
                bd["fragmented"],
                bottom=bd["complete"],
                color=s["color"],
                edgecolor="black",
                linewidth=0.5,
                alpha=0.45,
                hatch="//",
            )

        ax.set_title(format_species(sp), fontsize=9, fontstyle="italic")
        ax.set_ylim(0, 105)
        ax.set_xticks([])
        ax.grid(True, axis="y", alpha=0.2, linestyle="--")
        ax.tick_params(labelsize=7)

        if ci == 0:
            ax.set_ylabel("BUSCO [%]", fontsize=9)

    handles = [
        plt.Line2D(
            [],
            [],
            color=TOOL_STYLES[t]["color"],
            marker="s",
            linestyle="None",
            markersize=10,
            markeredgecolor="black",
            markeredgewidth=0.5,
            label=TOOL_STYLES[t]["label"],
        )
        for t in TOOLS_ORDER
    ]
    # Add completeness category entries
    handles.append(Patch(facecolor="white", edgecolor="black", label="Complete"))
    handles.append(
        Patch(
            facecolor="white",
            edgecolor="black",
            alpha=0.45,
            hatch="//",
            label="Fragmented",
        )
    )

    fig.legend(
        handles=handles,
        loc="lower center",
        ncol=3,
        fontsize=9,
        frameon=True,
        fancybox=True,
        edgecolor="gray",
    )

    plt.tight_layout(rect=[0.02, 0.10, 1, 1])
    fig.savefig(BUSCO_OUTPUT, dpi=300, bbox_inches="tight")
    print(f"Saved BUSCO figure   → {BUSCO_OUTPUT}")
    plt.close(fig)




def main():
    gff_data, busco_data, species_list = collect_data()

    print(
        f"Found {len(species_list)} genomes: "
        + ", ".join(format_species(s) for s in species_list)
    )
    tools_found = sorted(
        {t for _, t in list(gff_data.keys()) + list(busco_data.keys())}
    )
    print(f"Tools: {', '.join(TOOL_STYLES[t]['label'] for t in TOOLS_ORDER if t in tools_found)}")

    # ── Console summary ──────────────────────────────────────────────────
    print("\n── BUSCO Completeness ─────────────────────────────────────────")
    print(
        f"{'Species':<22} {'Tool':<26} {'Complete':>8} {'Frag':>6} {'Miss':>6}"
    )
    print("─" * 72)
    for sp in species_list:
        for tool in TOOLS_ORDER:
            if (sp, tool) in busco_data:
                bd = busco_data[(sp, tool)]
                print(
                    f"{format_species(sp):<22} "
                    f"{TOOL_STYLES[tool]['label']:<26} "
                    f"{bd['complete']:>7.1f}% "
                    f"{bd['fragmented']:>5.1f}% "
                    f"{bd['missing']:>5.1f}%"
                )

    print("\n── gffcompare Exon level ──────────────────────────────────────")
    print(
        f"{'Species':<22} {'Tool':<26} {'Sensitivity':>12} {'Precision':>10}"
    )
    print("─" * 72)
    for sp in species_list:
        for tool in TOOLS_ORDER:
            key = (sp, tool)
            if key in gff_data and "Exon level" in gff_data[key]:
                d = gff_data[key]["Exon level"]
                print(
                    f"{format_species(sp):<22} "
                    f"{TOOL_STYLES[tool]['label']:<26} "
                    f"{d['sensitivity']:>11.1f}% "
                    f"{d['precision']:>9.1f}%"
                )

    # ── Generate plots ───────────────────────────────────────────────────
    plot_gffcompare(gff_data, species_list)
    plot_busco(busco_data, species_list)
    print("\nDone!")


if __name__ == "__main__":
    main()
