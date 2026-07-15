from pathlib import Path
import re


ROOT = Path(__file__).resolve().parent
cap1 = (ROOT / "methods_paper_chapter1_wt_pifq.html").read_text()
cap2 = (ROOT / "methods_paper_chapter2_wt_pifq.html").read_text()

style1 = re.search(r"<style>(.*?)</style>", cap1, re.S).group(1)
style2 = re.search(r"<style>(.*?)</style>", cap2, re.S).group(1)
body1 = re.search(r"<body>(.*?)</body>", cap1, re.S).group(1).strip()
body2 = re.search(r"<body>(.*?)</body>", cap2, re.S).group(1).strip()


def replace_text(html: str, replacements: dict[str, str]) -> str:
    for old, new in replacements.items():
        html = html.replace(old, new)
    return html


def strip_cap2_front_matter(html: str) -> str:
    html = re.sub(r"^\s*<section class=\"title-page\">.*?</section>\s*", "", html, flags=re.S)
    match = re.search(r"(<section class=\"page-section\">\s*<h1>Section 13\b.*)$", html, re.S)
    if not match:
        raise RuntimeError("Could not find Section 13 in Chapter 2 HTML")
    return match.group(1).strip()


versions_match = re.search(
    r"\s*(<section class=\"page-section versions-section\">.*?</section>)\s*$",
    body1,
    re.S,
)
if not versions_match:
    raise RuntimeError("Could not find software versions section")

versions_section = versions_match.group(1)
body1 = body1[: versions_match.start()].rstrip()
body2 = strip_cap2_front_matter(body2)
body1 = re.sub(
    r"\s*<h2>Docker quick start</h2>.*?<div class=\"code-wrap\"><div class=\"code-label\">bash</div><pre><code>docker pull.*?</pre></div>\s*",
    "\n",
    body1,
    count=1,
    flags=re.S,
)

body1 = replace_text(
    body1,
    {
        "<h2>Working directory and repository</h2>":
            "<h1>Section 0.1 - Working directory and repository</h1>",
        "<h1>Raw data download</h1>":
            "<h1>Section 0.2 - Raw data download</h1>",
        "<h1>Reference preparation for Cell Ranger</h1>":
            "<h1>Section 0.3 - Reference preparation for Cell Ranger</h1>",
        "<h1>Cell Ranger count generation</h1>":
            "<h1>Section 0.4 - Cell Ranger count generation</h1>",
        "<h1>Chapter 1 in R</h1>":
            "<h1>Section 0.5 - Seurat workflow in R</h1>",
        "<h1>Step 1 - Initialization</h1>":
            "<h1>Section 0.6 - Initialization</h1>",
        "<h1>Step 2 - Sample manifest</h1>":
            "<h1>Section 0.7 - Sample manifest</h1>",
        "Chapter 1 - WT vs pifq tutorial: data download, Cell Ranger processing, and Seurat analysis":
            "WT vs pifq workflow: data download, Cell Ranger processing, Seurat analysis, pseudobulk DE, GO, heatmap and networks",
        "This document rewrites Chapter 1 as a practical tutorial for the":
            "This document presents a practical tutorial for the",
        "Cell Ranger outputs, Docker execution, and Chapter 1 results are":
            "Cell Ranger outputs, Docker execution, and workflow results are",
        "all placed under that same Git working directory. The document stops at\n      Chapter 1.":
            "all placed under that same Git working directory. The sections then continue without restarting the environment.",
        "all placed under that same Git working directory. The document stops at Chapter 1.":
            "all placed under that same Git working directory. The sections then continue without restarting the environment.",
        "The adapted Chapter 1 script is":
            "The adapted Seurat preprocessing script is",
        "The output folder\n      <code>resultados_wt/</code>":
            "The output folder\n      <code>resultados_wt/</code>",
        "Before launching Chapter 1, place":
            "Before launching the R workflow, place",
        "This keeps Chapter 1\n      inspectable":
            "This keeps the workflow\n      inspectable",
        "Chapter 1 uses the standard filtered matrices":
            "The workflow uses the standard filtered matrices",
        "The Chapter 1 run finished":
            "The Seurat preprocessing run finished",
        "adapting the chapter\n      to another pair or set":
            "adapting the workflow\n      to another pair or set",
    },
)

body2 = replace_text(
    body2,
    {
        "The hdWGCNA section is preserved from the original Chapter 2 structure.":
            "The hdWGCNA section is preserved from the original pseudobulk workflow structure.",
        "It uses the significant-gene table generated in Section 19.":
            "It uses the significant-gene table generated in Section 19.",
    },
)

versions_section = replace_text(
    versions_section,
    {
        "<h1>Software versions</h1>":
            "<h1>Supplementary software versions</h1>",
        "while this appendix keeps only the versions\n    needed to reproduce Chapter 1 without filling the PDF with long tables.":
            "while this supplementary section keeps only the versions needed to reproduce the complete workflow without filling the PDF with long tables.",
        "needed to reproduce Chapter 1":
            "needed to reproduce the complete workflow",
    },
)

html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>WT pifq Methods Paper - Continuous Workflow</title>
  <style>
{style1}

/* Pseudobulk sections */
{style2}
  </style>
</head>
<body>
{body1}
{body2}
{versions_section}
</body>
</html>
"""

(ROOT / "methods_paper_wt_pifq_chapter1_2.html").write_text(html)
print(ROOT / "methods_paper_wt_pifq_chapter1_2.html")
