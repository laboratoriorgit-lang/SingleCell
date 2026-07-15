from pathlib import Path
import re

ROOT = Path(__file__).resolve().parent
cap1 = (ROOT / 'methods_paper_chapter1_wt_pifq.html').read_text()
cap2 = (ROOT / 'methods_paper_chapter2_wt_pifq.html').read_text()
style1 = re.search(r'<style>(.*?)</style>', cap1, re.S).group(1)
style2 = re.search(r'<style>(.*?)</style>', cap2, re.S).group(1)
body1 = re.search(r'<body>(.*?)</body>', cap1, re.S).group(1).strip()
body2 = re.search(r'<body>(.*?)</body>', cap2, re.S).group(1).strip()
body2 = re.sub(r'^\s*<section class="title-page">.*?</section>\s*', '', body2, flags=re.S)
chapter2_intro = '''
  <section class="page-section step-section">
    <h1>Chapter 2 - Pseudobulk DE and GO enrichment</h1>
    <p>
      The following sections continue directly from Chapter 1 and preserve the
      structure of <code>workflow/capitulo2_pseudobulk_de.R</code>. The WT/pifq
      adaptation uses <code>pbmc_harmony_annotated.rds</code>, direct
      <code>celltype</code> labels, and the <code>WT_vs_pifq</code> contrast.
    </p>
  </section>
'''
html = f'''<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>WT pifq Methods Paper - Chapters 1 and 2</title>
  <style>
{style1}

/* Chapter 2 additions */
{style2}
  </style>
</head>
<body>
{body1}
{chapter2_intro}
{body2}
</body>
</html>
'''
(ROOT / 'methods_paper_wt_pifq_chapter1_2.html').write_text(html)
print(ROOT / 'methods_paper_wt_pifq_chapter1_2.html')
