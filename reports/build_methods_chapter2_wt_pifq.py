from __future__ import annotations

import html
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT.parent / "workflow" / "capitulo2_wt_pifq.R"
OUT = ROOT / "methods_paper_chapter2_wt_pifq.html"


def esc(value: object) -> str:
    return html.escape("" if value is None else str(value))


def code_block(language: str, text: str) -> str:
    return (
        '<div class="code-wrap">'
        f'<div class="code-label">{esc(language)}</div>'
        f"<pre><code>{esc(text.strip())}</code></pre>"
        "</div>"
    )


def img(filename: str, caption: str, cls: str = "") -> str:
    return (
        f'<figure class="{esc(cls)}">'
        f'<img src="assets/{esc(filename)}" alt="{esc(caption)}">'
        f"<figcaption>{esc(caption)}</figcaption>"
        "</figure>"
    )


def strip_comment_lines(text: str) -> str:
    lines = []
    for line in text.splitlines():
        if line.strip().startswith("#"):
            continue
        lines.append(line.rstrip())
    return "\n".join(lines).strip()


def section_code(script: str, start: str, end: str | None = None) -> str:
    start_idx = script.index(start)
    if end is None:
        chunk = script[start_idx:]
    else:
        end_idx = script.index(end, start_idx + len(start))
        chunk = script[start_idx:end_idx]
    return strip_comment_lines(chunk)


def build() -> str:
    script = SCRIPT.read_text()

    config = strip_comment_lines(script.split("# =============================================================================\n# LOAD")[0])
    load = section_code(script, "# LOAD", "# =============================================================================\n# SECTION 13")
    s13 = section_code(script, "# SECTION 13", "# SECTION 14")
    s14 = section_code(script, "# SECTION 14", "# SECTION 15")
    s15 = section_code(script, "# SECTION 15", "# SECTION 16")
    s16 = section_code(script, "# SECTION 16", "# SECTION 17")
    s17 = section_code(script, "# SECTION 17", "# SECTION 18")
    s18 = section_code(script, "# SECTION 18", "# SECTION 19")
    s19 = section_code(script, "# SECTION 19", "# SECTION 20")
    s20 = section_code(script, "# SECTION 20", "# SECTION 21")
    s21 = section_code(script, "# SECTION 21", "# SECTION 21b")
    s21b = section_code(script, "# SECTION 21b", None)

    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Chapter 2 WT pifq Methods Tutorial</title>
  <style>
    @page {{
      size: Letter;
      margin: 0.72in 0.78in;
      @bottom-center {{
        content: counter(page);
        color: #6b7280;
        font-size: 8.5pt;
      }}
    }}

    body {{
      color: #1f2933;
      font-family: Georgia, "Times New Roman", serif;
      font-size: 10.8pt;
      line-height: 1.43;
    }}

    .title-page {{
      margin-top: 0.55in;
    }}

    .title {{
      color: #123f5a;
      font-size: 23pt;
      font-weight: 700;
      line-height: 1.12;
      margin: 0 0 0.2em;
      text-align: center;
    }}

    .title em {{
      display: block;
      font-size: 20pt;
    }}

    .subtitle {{
      color: #2f6b49;
      font-size: 13.2pt;
      font-style: italic;
      margin: 1.25em 0 1.15em;
    }}

    .meta {{
      color: #374151;
      margin: 0.2em 0;
    }}

    .rule {{
      border-top: 2px solid #7ba7c0;
      margin: 1.05em 0 1.65em;
    }}

    h1, h2 {{
      color: #123f5a;
      font-family: Georgia, "Times New Roman", serif;
      font-weight: 700;
      line-height: 1.16;
      page-break-after: avoid;
    }}

    h1 {{
      border-bottom: 1.6px solid #a9c5d6;
      font-size: 17pt;
      margin: 1.2em 0 0.75em;
      padding-bottom: 0.22em;
    }}

    h2 {{
      border-bottom: 1px solid #b8ccd9;
      font-size: 13.3pt;
      margin: 1.05em 0 0.55em;
      padding-bottom: 0.15em;
    }}

    p {{
      margin: 0.42em 0 0.62em;
      text-align: justify;
    }}

    code {{
      color: #174a67;
      font-family: "SFMono-Regular", Menlo, Consolas, monospace;
      font-size: 9.4pt;
    }}

    .code-wrap {{
      border: 1px solid #b7b7b7;
      margin: 0.62em 0 1.05em;
      page-break-inside: auto;
      position: relative;
    }}

    .code-label {{
      background: white;
      color: #2f6b49;
      font-size: 8.2pt;
      font-style: italic;
      left: 1.1em;
      padding: 0 0.35em;
      position: absolute;
      top: -0.68em;
    }}

    pre {{
      color: #174a67;
      font-family: "SFMono-Regular", Menlo, Consolas, monospace;
      font-size: 9.2pt;
      line-height: 1.33;
      margin: 0;
      overflow-wrap: break-word;
      padding: 0.72em 0.85em;
      white-space: pre-wrap;
    }}

    .page-section {{
      break-before: page;
    }}

    .small-note {{
      color: #4b5563;
      font-size: 9.4pt;
    }}

    figure {{
      margin: 0.45em auto 0.75em;
      page-break-inside: avoid;
      text-align: center;
    }}

    figure img {{
      display: block;
      height: auto;
      margin: 0 auto;
      max-height: 3.25in;
      max-width: 50%;
    }}

    figure.wide img {{
      max-height: 3.5in;
      max-width: 64%;
    }}

    figure.tall img {{
      max-height: 4.6in;
      max-width: 64%;
    }}

    .two-col {{
      display: grid;
      gap: 0.35in;
      grid-template-columns: 1fr 1fr;
      page-break-inside: avoid;
    }}

    .two-col figure img {{
      max-height: 2.85in;
      max-width: 94%;
    }}

    figcaption {{
      color: #4b5563;
      font-size: 8.9pt;
      margin-top: 0.35em;
    }}
  </style>
</head>
<body>
  <section class="title-page">
    <h1 class="title">Computational Methods for Single-Cell RNA-Seq Analysis of<br><em>Arabidopsis thaliana</em></h1>
    <p class="subtitle">Chapter 2 - WT vs pifq: pseudobulk DE, GO enrichment, heatmaps and network sections</p>
    <p class="meta">SingleCell Pipeline</p>
    <p class="meta">2026-07-15</p>
    <div class="rule"></div>

    <h1>Purpose</h1>
    <p>
      This Chapter 2 draft follows the structure of
      <code>workflow/capitulo2_pseudobulk_de.R</code>. The WT/pifq adaptation
      changes only the dataset-specific pieces: output folder
      <code>resultados_wt</code>, input object
      <code>pbmc_harmony_annotated.rds</code>, direct annotation column
      <code>celltype</code>, and the <code>WT_vs_pifq</code> contrast.
    </p>
    <p>
      No grouped or curated labels are used here. This keeps Chapter 2 aligned
      with the direct bibliography annotation from Chapter 1.
    </p>
  </section>

  <section class="page-section">
    <h1>Configuration</h1>
    <p>
      Configuration mirrors the original Chapter 2 script, with the WT/pifq
      result directory and the same function sources.
    </p>
    {code_block("R", config)}
  </section>

  <section class="page-section">
    <h1>Load - Annotated Object From Part 1</h1>
    <p>
      The input is the Chapter 1 bibliography-annotated object. The active
      annotation column for Chapter 2 is <code>celltype</code>.
    </p>
    {code_block("R", load)}
  </section>

  <section class="page-section">
    <h1>Section 13 - Cell-Type Subsets</h1>
    <p>
      One Seurat subset is created per direct <code>celltype</code> label.
      Object names are sanitized only for list names and filenames.
    </p>
    {code_block("R", s13)}
  </section>

  <section class="page-section">
    <h1>Section 14 - Pseudo-Replicate Assignment</h1>
    <p>
      The same pseudo-replicate assignment structure is retained. The QC check
      uses <code>Mesophyll</code> as an existing label in the WT/pifq dataset.
    </p>
    {code_block("R", s14)}
  </section>

  <section class="page-section">
    <h1>Section 15 - Pseudobulk Tables and DESeq2</h1>
    <p>
      The contrast is <code>WT_vs_pifq</code>. Positive log2FC values correspond
      to higher expression in pifq relative to WT.
    </p>
    {code_block("R", s15)}
  </section>

  <section class="page-section">
    <h1>Section 16 - Volcano Plots</h1>
    {code_block("R", s16)}
    <div class="two-col">
      {img("cap2_volcano_mesophyll_WT_vs_pifq.png", "Representative volcano output: Mesophyll.")}
      {img("cap2_volcano_epidermis_hypocotyl_1_WT_vs_pifq.png", "Representative volcano output: Epidermis Hypocotyl.1.")}
    </div>
  </section>

  <section class="page-section">
    <h1>Section 17 - Differential Gene Tables</h1>
    {code_block("R", s17)}
  </section>

  <section class="page-section">
    <h1>Section 18 - GO Enrichment</h1>
    {code_block("R", s18)}
    {img("cap2_go_mesophyll_bp.png", "Representative GO biological-process bubble plot for Mesophyll.", "wide")}
  </section>

  <section class="page-section">
    <h1>Section 19 - Log2FC Heatmap</h1>
    {code_block("R", s19)}
    {img("cap2_heatmap_WT_vs_pifq.png", "log2FC heatmap produced by Section 19.", "tall")}
  </section>

  <section class="page-section">
    <h1>Section 20 - Coexpression Network</h1>
    <p>
      The hdWGCNA section is preserved from the original Chapter 2 structure.
      It uses the significant-gene table generated in Section 19.
    </p>
    {code_block("R", s20)}
  </section>

  <section class="page-section">
    <h1>Section 21 - Network Export and Visualization</h1>
    {code_block("R", s21)}
  </section>

  <section class="page-section">
    <h1>Section 21b - TF Coexpression Network</h1>
    {code_block("R", s21b)}
  </section>
</body>
</html>
"""


if __name__ == "__main__":
    OUT.write_text(build(), encoding="utf-8")
    print(OUT)
