#!/usr/bin/env python3

# GEP2 Genome Stats Report Generator across all Assembly-Entries
# 2026
# This script is part of the GEP2 pipeline

# Generates a markdown / html report aggregating relevant results [with most important metrics and thresholds based on standards by EBP etc. over all species/asm_id-Entrys of one GEP2-Pipeline-Run to heatmap-style Qualtyrating]
# [ metric 1-11 explained ]
#   (e. g. "- contig N50 (gfastats) (assembly metrics)")

# Status: Skeleton-Version
# Uses whenever possible the existing parsing functions from make_gep2_repors.py


from __future__ import annotations
from pathlib import Path

import yaml # type: ignore
import argparse
import os
import sys
from typing import Any, Optional

# Import existing parsers.
from make_gep2_report import (
    parse_gfastats,
    parse_compleasm,
    parse_compleasm_full,
    parse_busco_summary,
    parse_merqury_qv,
    parse_merqury_completeness,
    parse_inspector,
    detect_kmer_tool,
    get_species_genomic_data_from_goat,
)

__version__ = "0.1.0"


# -------------------------------------------------------------------------------
# TYPE ALIASES
# -------------------------------------------------------------------------------
# Ein paar Alias-Typen, um die Signaturen unten lesbarer zu halten. Besonders
# SuperTable ist noch bewusst vage (Any) gehalten, bis Section 4 entschieden hat,
# ob das ein pandas.DataFrame, eine Liste von Dicts, o.ae. wird.

AssemblyKey = tuple[str, str]          # (species, asm_id)
AssemblyData = dict[str, Any]          # Rohwerte EINER Assembly
ThresholdSpec = dict[str, Any]         # ein Eintrag aus THRESHOLDS
SuperTable = Any                       # TODO: konkretisieren, sobald Section 4 steht


# -------------------------------------------------------------------------------
# HELPER-FUNCTIONS
# -------------------------------------------------------------------------------

def _is_reads_only_entry(species, asm_id, samples_config) -> bool:
    """Check if this is a reads-only entry (no assembly, just reads for profiling)."""
    if asm_id is None or str(asm_id).lower() in ("none", "", "-"):
        return True
    
    try:
        asm_data = samples_config["sp_name"][species]["asm_id"].get(asm_id, {})
        asm_files = asm_data.get("asm_files", {})
        
        # Check if all asm_files are None/empty
        has_assembly = False
        for path_key, path_value in asm_files.items():
            if path_value and str(path_value).lower() not in ("none", "", "-"):
                has_assembly = True
                break
        
        return not has_assembly
        
    except (KeyError, TypeError, AttributeError):
        return True
    
def get_assembly_files(species, assembly, samples_config):
    """Get all assembly files for a given species/assembly"""
    try:
        asm_data = samples_config["sp_name"][species]["asm_id"][assembly]
        asm_files = asm_data.get("asm_files", {})
        return {k: v for k, v in asm_files.items() if v and v != "None"}
    except (KeyError, TypeError, AttributeError):
        return {}
    
def get_assembly_basename(filepath):
    """Extract assembly basename from filepath"""
    basename = os.path.basename(filepath)
    for ext in ['.fna.gz', '.fa.gz', '.fasta.gz', '.fna', '.fa', '.fasta', '.gz']:
        if basename.endswith(ext):
            basename = basename[:-len(ext)]
    return basename
    
# -------------------------------------------------------------------------------
# SECTION 1: THRESHOLDS / KONFIGURATION
# -------------------------------------------------------------------------------

THRESHOLDS: dict[str, ThresholdSpec] = {
    "contig_n50": {
        "source": "EBP Report v7 (2026)",
        # To-Do: exakte Grenzen aus 6.C.Q40-Notation ableiten
    },
    "scaffold_n50": {
        "source": "EBP Report v7 (2026)",
    },
    "gaps_per_gbp": {
        "source": "TODO - VGP-Zitatstelle noch verifizieren (siehe Recherche-Notiz)",
    },
    "structural_errors": {
        "source": None,  # kein Konsens-Schwellenwert in der Literatur gefunden
    },
    "false_duplications_pct": {
        "source": "EBP Report v7 (2026)",
        "threshold": 5.0,  # < 5%
    },
    "merqury_qv": {
        "source": "EBP Report v7 (2026) / VGP-Notation",
        "threshold": 40.0,  # Q40
    },
    "kmer_completeness_pct": {
        "source": "EBP Report v7 (2026)",
        "threshold": 90.0,  # > 90%
    },
    "busco_single_pct": {
        "source": "EBP Report v7 (2026)",
        "threshold": 90.0,  # > 90%
    },
    "busco_duplicated_pct": {
        "source": "Pipeline-interne Heuristik, konzeptionell gestuetzt durch "
                   "Simao et al. 2015 ('should be rare'), kein exakter "
                   "Literatur-Zahlenwert - siehe Limitations-Kapitel",
    },
    "chrom_seq_pct": {
        "source": "EBP Report v7 (2026)",
        "threshold": 90.0,  # > 90%
    },
    "l90_haploid_proxy": {
        "source": "GEP2-intern, Proxy fuer chrom_seq_pct (siehe Herleitung in Thesis)",
    },
}


# -------------------------------------------------------------------------------
# SECTION 2: DISCOVERY - alle Assembly-Eintraege eines Pipeline-Laufs finden
# -------------------------------------------------------------------------------

def load_control_panel(path: str):
    control_panel_path = Path(path)
    if not control_panel_path.exists():
        raise FileNotFoundError(
            f"control_panel.yaml not not found in: {control_panel_path}"
        )
        
    with control_panel_path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f) or {}

def load_samples_config(path: str):
    
    config_path = Path(path) / "data_config.yaml"
    if not config_path.exists():
        raise FileNotFoundError(
            f"data_config.yaml not found in path: {config_path}"
        )

    with config_path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f) or {}

def discover_assemblies(config: dict) -> list[AssemblyKey]:
    """
    [...]
    """

    assemblies: list[AssemblyKey] = []

    for species, species_data in config["sp_name"].items():
        if not species_data or "asm_id" not in species_data:
            continue
        
        for assembly, assembly_data in species_data["asm_id"].items():
            if not assembly_data:
                continue
            
            # Skip reads-only entries (no assembly)
            if _is_reads_only_entry(species, assembly, config):
                continue
            
            assemblies.append((species, assembly))
            
    return assemblies
    
    


# -------------------------------------------------------------------------------
# SECTION 3: DATENSAMMLUNG PRO ASSEMBLY
# -------------------------------------------------------------------------------
# Eine Funktion pro Assembly-Eintrag, die die bestehenden Parser 1:1 wiederverwendet.
# Gibt ein flaches Dict mit Rohwerten zurueck - KEIN Rating/Styling passiert hier.

def collect_assembly_data(
    species: str,
    asm_id: str,
    out_folder: str,
    config: dict[str, Any],
) -> AssemblyData:
    """
    Sammelt alle Rohwerte fuer eine species/asm_id-Kombination.

    To-Do je Metrik:
      - gfastats (#1,#2,#3,#11): direkt parse_gfastats() wiederverwenden
      - compleasm/BUSCO (#8,#9): PREFER_BUSCO-Flag beachten, genau wie
        get_report_busco_inputs()/get_report_compleasm_inputs() in
        Z_reporting.smk das tun
      - Merqury (#5,#6,#7): kmer_tool-Erkennung (MERQ vs MERQ.FK) exakt wie
        im bestehenden Skript uebernehmen
      - False Duplications (#5): NOCH NICHT IMPLEMENTIERT
        -> neue Parsing-/Berechnungslogik noetig (aus Merqury .hist/spectra-cn
        ableiten, oder BUSCO-Duplicated-Rate als Naeherung verwenden -
        siehe Recherche-Notizen)
      - Inspector (#4): parse_inspector() wiederverwenden, aber KEIN
        Sterne-Rating anwenden (kein Threshold vorhanden - Rohwert anzeigen)
      - Chromosomen-% (#10): NOCH NICHT IMPLEMENTIERT
        -> neue Berechnung: Summe der n groessten Scaffolds / Gesamtlaenge,
        n = haploide Chromosomenzahl aus GoaT. Unabhaengig vom Hi-C-Modul
        halten (RUN_HIC kann False sein), z.B. eigenen seqkit-Aufruf auf die
        Assembly-Datei direkt.
      - L90-Proxy (#11): parse_gfastats() + get_species_genomic_data_from_goat()
    """
    
    data: AssemblyData = {"species": species, "asm_id": asm_id}
    
    # ---------- 1. gfastats: -------------------------------------------------
    asm_files = get_assembly_files(species, asm_id, config)
    gfastats_list = []
    for asm_key, asm_path in sorted(asm_files.items()):
        asm_basename = get_assembly_basename(filepath=asm_path)
        stats_path = Path(out_folder) / species / asm_id / "gfastats" / f"{asm_basename}_stats.txt"
        if stats_path.exists():
            gfastats_list.append(parse_gfastats(filepath=str(stats_path)))
        else:
            gfastats_list.append({})
    
    data["gfastats"] = gfastats_list
    
    # -------------------------------------------------------------------------
    
    return data


# -------------------------------------------------------------------------------
# SECTION 4: AGGREGATION UEBER ALLE ASSEMBLIES
# -------------------------------------------------------------------------------

def build_super_table(all_assembly_data: list[AssemblyData]) -> SuperTable:
    """
    Fuehrt die Liste der Pro-Assembly-Dicts zu einer aggregierten
    Tabellenstruktur zusammen.

    TODO: Datenstruktur festlegen. Ein pandas DataFrame wuerde CSV/Excel-Export
    und spaeteres Heatmap-Styling (z.B. via pandas Styler oder openpyxl)
    deutlich erleichtern gegenueber selbstgebauten Listen/Dicts - diese
    Entscheidung lohnt sich VOR Section 5/6. Sobald entschieden, den
    SuperTable-Alias oben von Any auf den konkreten Typ aendern
    (z.B. "pandas.DataFrame").
    """
    raise NotImplementedError


# -------------------------------------------------------------------------------
# SECTION 5: RATING / GUT-MITTEL-SCHLECHT-KLASSIFIZIERUNG
# -------------------------------------------------------------------------------
# Bewusst getrennt von der Ausgabe (Section 6), damit die Klassifizierungslogik
# unabhaengig vom spaeteren Ausgabeformat getestet werden kann.

def rate_value(
    metric_name: str,
    value: Optional[float],
    **context: Any,
) -> str:
    """
    Bildet einen Rohwert auf eine Gut/Mittel/Schlecht- (oder ****/***-/etc.)
    Klassifizierung ab, basierend auf THRESHOLDS oben.

    TODO - SPAETERE AUFGABE: Das ist der Kern der noch offenen
    Heatmap-Farblogik. Nicht noetig, solange die Metrik-Recherche/Thresholds
    noch nicht 100% final sind - erstmal nur Stub.
    context kann z.B. haploid_number enthalten (fuer l90_haploid_proxy).

    Rueckgabewert bewusst als str (z.B. "****", "···-", oder "N/A") statt als
    Enum/Klasse - haelt die Schnittstelle zu render_* simpel; falls spaeter
    mehr Struktur noetig ist (z.B. fuer Farbcodes), liesse sich hier ohne
    Bruch auf ein eigenes Rating-Objekt umstellen.
    """
    raise NotImplementedError


# -------------------------------------------------------------------------------
# SECTION 6: RENDERING / AUSGABE
# -------------------------------------------------------------------------------
# Format noch offen (Markdown+HTML? reines HTML? Bild?) - als eigene,
# austauschbare Schicht gehalten, damit diese Entscheidung Section 1-4 nicht
# blockiert.

def render_markdown_table(super_table: SuperTable) -> str:
    """Reine Markdown-Tabelle ohne Farben - guter erster Meilenstein/Fallback."""
    raise NotImplementedError


def render_html_heatmap(super_table: SuperTable) -> str:
    """
    TODO - SPAETERE AUFGABE: HTML-Tabelle mit farbcodierten Zellen
    (Gruen/Gelb/Rot). Fuer die erste lauffaehige Version nicht noetig - Stub.
    """
    raise NotImplementedError


# -------------------------------------------------------------------------------
# SECTION 7: CLI ENTRY POINT
# -------------------------------------------------------------------------------

def main() -> None:
    # --------------- CLI - Abschnitt: --------------------
    parser = argparse.ArgumentParser(description="Generates a GEP2 Super Report for all assemblies in a pipeline run.")
    
    parser.add_argument(
        "--control-panel", required=True,
        help="Path to 'control_panel.yaml'"
    )
    parser.add_argument(
        "--output", required=True,
        help="Output file for the aggregated super report"
    )
    parser.add_argument(
        "--format", choices=["markdown", "html"], default="markdown",
        help="Output Format"
    )
    
    args: argparse.Namespace = parser.parse_args()
    # -----------------------------------------------------
    
    control_panel = load_control_panel(path=args.control_panel)
    
    OUT_FOLDER = Path(control_panel["OUT_FOLDER"]) / "GEP2_results" 
    
    samples_config = load_samples_config(path=OUT_FOLDER)


    assemblies: list[AssemblyKey] = discover_assemblies(samples_config)
    # e.g.:     assemblies = [('Caenorhabditis_afra', 'nxCaeAfra1.1'), ('Russula_nigricans', 'gfRusNigr1.hap')]
    
    print(assemblies)
    
    # Über alle ASM's iterieren: Parsing-Funktionen für alle Werte
    all_data: list[AssemblyData] = [
        collect_assembly_data(sp, asm_id, str(OUT_FOLDER), config=samples_config)
        for sp, asm_id in assemblies
    ]
    
    print(all_data)
    
    # # Liste des ASM-Daten in Tabelle übertragen:
    # super_table: SuperTable = build_super_table(all_data)

    # # Output:
    # output_text: str
    
    # # Output-Text generieren:
    # if args.format == "markdown":
    #     output_text = render_markdown_table(super_table)
    # else:
    #     output_text = render_html_heatmap(super_table)

    # # Output-Text in File schreiben:
    # with open(args.output, "w") as f:
    #     f.write(output_text)


if __name__ == "__main__":
    main()