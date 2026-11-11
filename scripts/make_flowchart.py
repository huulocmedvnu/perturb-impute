"""Figure 1: study design as an SVG schematic, exported to PNG (600 dpi) and PDF.
Each row: data (left), a pictogram of the comparison (middle) and a small plot of the readout (right).
Readout marks use numbers.json. Canvas units: 10 per mm, 180 mm wide, 8 pt is about 28 units."""
import json, math, os, random, subprocess
import sys
HERE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "source_data")
F = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "output"); os.makedirs(F, exist_ok=True)
N = json.load(open(os.path.join(HERE, "numbers.json")))
W = 1800; INK, GREY, LINE = "#111111", "#8C96A8", "#DCE1EA"
BLUE, RED, GREEN, PURPLE, ORANGE = "#636EFA", "#EF553B", "#00CC96", "#AB63FA", "#FFA15A"
FONT = "Arial, Helvetica, sans-serif"; svg = []; add = svg.append
def text(x, y, s, size=27, weight="normal", anchor="start", fill=INK):
    add(f"<text x='{x:.1f}' y='{y:.1f}' font-family='{FONT}' font-size='{size}' font-weight='{weight}' text-anchor='{anchor}' fill='{fill}'>{s}</text>")
def line(x1, y1, x2, y2, col=INK, w=3, dash=None, head=False):
    d = f" stroke-dasharray='{dash}'" if dash else ""; h = " marker-end='url(#ah)'" if head else ""
    add(f"<line x1='{x1:.1f}' y1='{y1:.1f}' x2='{x2:.1f}' y2='{y2:.1f}' stroke='{col}' stroke-width='{w}'{d}{h}/>")
def rect(x, y, w, h, fill, op=1.0, stroke=None):
    s = f" stroke='{stroke}' stroke-width='2'" if stroke else ""
    add(f"<rect x='{x:.1f}' y='{y:.1f}' width='{w:.1f}' height='{h:.1f}' fill='{fill}' fill-opacity='{op}'{s}/>")
def cell(x, y, col, r=12):
    add(f"<circle cx='{x:.1f}' cy='{y:.1f}' r='{r}' fill='{col}' fill-opacity='0.25' stroke='{col}' stroke-width='3'/><circle cx='{x:.1f}' cy='{y:.1f}' r='{r * 0.35:.1f}' fill='{col}'/>")
def axes(x0, y0, w, h, ylab=None):
    line(x0, y0, x0 + w, y0, INK, 2.5); line(x0, y0, x0, y0 - h, INK, 2.5)
    if ylab: add(f"<text transform='translate({x0 - 14},{y0 - h / 2}) rotate(-90)' font-family='{FONT}' font-size='24' text-anchor='middle' fill='{INK}'>{ylab}</text>")
def grid(x, y, nr, nc, cw, ch, cols_on, col, missing=False):
    """matrix of perturbations (rows) x genes (columns); cols_on highlighted, others grey or hatched"""
    for i in range(nr):
        for j in range(nc):
            on = j in cols_on
            if missing and not on: rect(x + j * cw, y + i * ch, cw - 2, ch - 2, "white", 1, stroke=GREY)
            else: rect(x + j * cw, y + i * ch, cw - 2, ch - 2, col if on else GREY, 0.85 if on else 0.3)
add("<svg xmlns='http://www.w3.org/2000/svg' width='180mm' height='@HMM@mm' viewBox='0 0 1800 @H@'>")
add(f"<defs><marker id='ah' viewBox='0 0 10 10' refX='8' refY='5' markerWidth='4' markerHeight='4' orient='auto'><path d='M0,0 L10,5 L0,10 Z' fill='{GREY}'/></marker></defs>")
add("<rect width='1800' height='@H@' fill='white'/>")
C1, C2, C3 = 40, 560, 1250
text(C1, 58, "Data", 34, "bold"); text(C2, 58, "Comparison", 34, "bold"); text(C3, 58, "Readout", 34, "bold")
RH, Y0, FIGX = 260, 90, 1760
def row_frame(k, col, title, lines_, fig):
    y = Y0 + k * RH
    rect(C1 - 20, y + 12, 8, RH - 34, col); text(C1, y + 44, title, 29, "bold", fill=col)
    for i, s in enumerate(lines_): text(C1, y + 84 + i * 34, s, 26)
    text(FIGX, y + 44, fig, 27, "bold", "end", col)
    if k: line(C1, y, W - 40, y, LINE, 2)
    return y
PAN = {1, 4, 7}
# ---- row 1: reference atlas, panel, imputation
k = N["K562"]
y = row_frame(0, BLUE, "Reference atlas", ["Replogle, K562 genome-wide", f"{N['noise']['n_perts_gwps']:,} perturbations", f"{N['noise']['n_genes_universe']:,} genes"], "Fig. 2, 3")
grid(C2, y + 70, 5, 9, 22, 22, PAN, BLUE)
text(C2 + 99, y + 210, "choose panel genes", 24, anchor="middle")
line(C2 + 215, y + 125, C2 + 285, y + 125, GREY, 4, head=True)
grid(C2 + 300, y + 70, 5, 9, 22, 22, PAN, BLUE, missing=True)
text(C2 + 399, y + 210, "new screen: panel measured,", 24, anchor="middle"); text(C2 + 399, y + 238, "other genes imputed", 24, anchor="middle")
bx, by, bh = C3 + 40, y + 205, 130
axes(bx, by, 330, bh, ylab="recovered")
for i, (lab, v, col) in enumerate([("greedy", k["ratio_greedy_680"], BLUE), ("L1000", k["ratio_l1000_680"], RED), ("random", k["ratio_random1_680"], GREY)]):
    rect(bx + 25 + i * 100, by - v * bh, 64, v * bh, col, 0.85); text(bx + 57 + i * 100, by + 30, lab, 24, anchor="middle")
    text(bx + 57 + i * 100, by - v * bh - 10, f"{v:.2f}", 24, anchor="middle")
# ---- row 2: four cell lines
y = row_frame(1, RED, "Four cell lines", ["Replogle: K562, RPE1", "Nadig: HepG2, Jurkat", "Six panel selectors"], "Fig. 2, 3")
cx = C2 + 70; text(cx, y + 70, "K562", 24, "bold", "middle")
for i in range(4): cell(cx - 30 + (i % 2) * 60, y + 115 + (i // 2) * 50, BLUE)
text(cx, y + 220, "reference", 24, anchor="middle")
line(cx + 75, y + 140, C2 + 240, y + 140, GREY, 4, head=True); text(C2 + 175, y + 125, "transfer", 22, anchor="middle")
for j, (lab, col) in enumerate([("RPE1", RED), ("HepG2", GREEN), ("Jurkat", PURPLE)]):
    x = C2 + 300 + j * 130
    text(x, y + 70, lab, 24, "bold", "middle", col)
    for i in range(4): cell(x - 25 + (i % 2) * 50, y + 115 + (i // 2) * 50, col)
bx, by, bh, bw = C3 + 40, y + 205, 130, 330
axes(bx, by, bw, bh, ylab="recovered")
cts = ["K562", "RPE1", "HepG2", "Jurkat"]; lk = {"K562": "680", "RPE1": "680", "HepG2": "622", "Jurkat": "622"}
for i, ct in enumerate(cts):
    x = bx + 40 + i * 85
    for v, col in ((N[ct]["ratio_greedy_500"], BLUE), (N[ct][f"ratio_l1000_{lk[ct]}"], RED)):
        yy = by - (v - 0.4) / 0.55 * bh; add(f"<path d='M{x},{yy - 10} L{x + 10},{yy} L{x},{yy + 10} L{x - 10},{yy} Z' fill='{col}'/>")
    text(x, by + 30, ct, 22, anchor="middle")
text(bx + bw + 16, by - 105, "greedy", 24, fill=BLUE); text(bx + bw + 16, by - 30, "L1000", 24, fill=RED)
# ---- row 3: DE presence, direction, size
m = N["hybrid"]["K562"]
y = row_frame(2, GREEN, "Differential expression", ["Same screens, cells split", "into two independent halves", "DE calls, complexes, pathways"], "Fig. 4, 5")
for i in range(3): cell(C2 + 30 + i * 36, y + 100, GREEN)
for i in range(3): cell(C2 + 30 + i * 36, y + 170, GREEN)
text(C2 + 175, y + 106, "half A: panel measured, rest imputed", 24); text(C2 + 175, y + 176, "half B: all genes measured", 24)
line(C2 + 10, y + 135, C2 + 520, y + 135, GREY, 2, "6 5")
text(C2 + 175, y + 225, "DE call, sign and size compared", 24)
bx, by, bh = C3 + 40, y + 205, 130
axes(bx, by, 330, bh, ylab="fraction")
for i, (lab, a, b) in enumerate([("DE found", m["direct_all"]["recall"], m["hybrid_all"]["recall"]), ("wrong sign", m["direct_all"]["opp"], m["hybrid_all"]["opp"])]):
    x = bx + 30 + i * 160
    rect(x, by - a / 0.4 * bh, 55, a / 0.4 * bh, GREY, 0.8); rect(x + 60, by - b / 0.4 * bh, 55, b / 0.4 * bh, GREEN, 0.85)
    text(x + 57, by + 30, lab, 24, anchor="middle")
text(bx + 346, by - 105, "measured", 24, fill=GREY); text(bx + 346, by - 75, "hybrid", 24, fill=GREEN)
# ---- row 4: equal sequencing budget
S = N["sim"]
y = row_frame(3, PURPLE, "Equal sequencing budget", ["Replogle, K562 essential genes", f"{S['n_cells_total']:,} cells, reads removed", "at random (2% to 100%)"], "Fig. 6")
random.seed(2)
hs = [random.random() for _ in range(18)]
for j, h in enumerate(hs): rect(C2 + j * 13, y + 190 - 18 - h * 40, 10, 18 + h * 40, GREY, 0.6)
text(C2 + 117, y + 225, "whole transcriptome", 24, anchor="middle")
for j in range(6): rect(C2 + 300 + j * 30, y + 190 - 40 - hs[j] * 70 - 30, 24, 40 + hs[j] * 70 + 30, PURPLE, 0.8)
text(C2 + 385, y + 225, "same reads on panel genes", 24, anchor="middle")
bx, by, bh, bw = C3 + 40, y + 205, 130, 330
axes(bx, by, bw, bh, ylab="specific r")
Q = ["0.02", "0.05", "0.1", "0.25", "1.0"]; px = lambda q: bx + 20 + (bw - 40) * (math.log10(float(q)) + 1.7) / 1.7
for key, col in (("direct", GREY), ("greedyfull_500", PURPLE)):
    vals = [S[q]["direct_rp"] if key == "direct" else S[q][key + "|targeted"]["rp"] for q in Q]
    add("<polyline points='" + " ".join(f"{px(q):.1f},{by - max(v, 0) / 0.1 * bh:.1f}" for q, v in zip(Q, vals)) + f"' fill='none' stroke='{col}' stroke-width='4'/>")
    for q, v in zip(Q, vals): add(f"<circle cx='{px(q):.1f}' cy='{by - max(v, 0) / 0.1 * bh:.1f}' r='6' fill='{col}'/>")
text(px("0.02"), by + 30, "2%", 22, anchor="middle"); text(px("1.0"), by + 30, "100%", 22, anchor="middle"); text(bx + bw / 2, by + 30, "reads", 22, anchor="middle")
text(bx + bw + 16, by - 110, "panel", 24, fill=PURPLE); text(bx + bw + 16, by - 40, "whole", 24, fill=GREY)
PINK = "#FF6692"
# ---- row 5: paired whole-transcriptome and L1000-enriched libraries of the same cells
CSV = {}
for l in open(os.path.join(HERE, "real", "cells.summary.tsv")):
    t = l.rstrip("\n").split("\t")
    if len(t) == 6 and t[0] != "m": CSV[(t[0], t[1])] = float(t[3])
y = row_frame(4, PINK, "Paired libraries", ["Same K562 cells sequenced twice:", "whole transcriptome and", "L1000-enriched library"], "Fig. 7")
random.seed(11); hp = [random.random() for _ in range(18)]
for j, h in enumerate(hp): rect(C2 + j * 13, y + 190 - 18 - h * 40, 10, 18 + h * 40, PINK if j % 6 == 2 else GREY, 0.8 if j % 6 == 2 else 0.6)
text(C2 + 117, y + 225, "whole transcriptome", 24, anchor="middle")
for j in range(6): rect(C2 + 300 + j * 30, y + 190 - 110 - hp[j] * 30, 24, 110 + hp[j] * 30, PINK, 0.8)
text(C2 + 385, y + 225, "enriched landmark genes", 24, anchor="middle")
bx, by, bh, bw = C3 + 40, y + 205, 130, 330
axes(bx, by, bw, bh, ylab="specific r")
MSC = ["25", "50", "100", "200", "all"]; pxm = lambda i: bx + 25 + (bw - 50) * i / 4
for key, col, dash in (("direct", GREY, None), ("l1000_enriched", PINK, None), ("l1000_extract", PINK, "8,6")):
    vals = [CSV[(m, key)] for m in MSC]
    d = f" stroke-dasharray='{dash}'" if dash else ""
    add("<polyline points='" + " ".join(f"{pxm(i):.1f},{by - v / 0.2 * bh:.1f}" for i, v in enumerate(vals)) + f"' fill='none' stroke='{col}' stroke-width='4'{d}/>")
    for i, v in enumerate(vals): add(f"<circle cx='{pxm(i):.1f}' cy='{by - v / 0.2 * bh:.1f}' r='6' fill='{col}'/>")
text(pxm(0), by + 30, "25", 22, anchor="middle"); text(pxm(4), by + 30, "all", 22, anchor="middle"); text(bx + bw / 2, by + 30, "input cells", 22, anchor="middle")
text(bx + bw + 16, by - 110, "whole", 24, fill=GREY); text(bx + bw + 16, by - 60, "enriched", 24, fill=PINK)
# ---- row 6: independent atlases (X-Atlas/Orion) and pooled reference
def cvrec(f, panel):
    for l in open(os.path.join(HERE, "xatlas", f)):
        t = l.split()
        if t and t[0] == panel: return float(t[2])
XA = {ct: (cvrec(f"{ct.lower()}_cv.txt", "greedy_500"), cvrec(f"{ct.lower()}_cv.txt", "random1_500")) for ct in ("HCT116", "HEK293T")}
y = row_frame(5, ORANGE, "Independent atlases", ["X-Atlas/Orion: HCT116, HEK293T", "7.9 million cells, different protocol", "Atlases pooled as one reference"], "Fig. 8")
for j, (lab, col) in enumerate([("K562", BLUE), ("HCT116", ORANGE), ("HEK293T", "#B6801F")]):
    x = C2 + j * 120
    text(x + 44, y + 70, lab, 22, "bold", "middle", col)
    grid(x, y + 85, 4, 4, 22, 22, {0, 2}, col)
for j in range(3): line(C2 + 44 + j * 120, y + 178, C2 + 44 + j * 120, y + 200, GREY, 2.5)
line(C2 + 44, y + 200, C2 + 284, y + 200, GREY, 2.5); line(C2 + 284, y + 200, C2 + 390, y + 200, GREY, 4, head=True)
rect(C2 + 400, y + 150, 120, 90, ORANGE, 0.25, stroke=ORANGE)
text(C2 + 460, y + 190, "pooled", 22, anchor="middle"); text(C2 + 460, y + 218, "reference", 22, anchor="middle")
bx, by, bh = C3 + 40, y + 205, 130
axes(bx, by, 330, bh, ylab="recovered")
for i, ct in enumerate(("HCT116", "HEK293T")):
    x = bx + 30 + i * 160
    for k2, (v, col) in enumerate(((XA[ct][0], ORANGE), (XA[ct][1], GREY))):
        rect(x + k2 * 60, by - v * bh, 55, v * bh, col, 0.85)
    text(x + 57, by + 30, ct, 22, anchor="middle")
text(bx + 346, by - 105, "greedy", 24, fill=ORANGE); text(bx + 346, by - 75, "random", 24, fill=GREY)
# ---- row 7: profiles transferred from other cell lines
VC = {}
for l in open(os.path.join(HERE, "vcm", "vcmsum.txt")):
    if l.startswith("== scope all"): break
    t = l.split()
    if len(t) > 5 and t[0] == "RPE1": VC[t[1]] = (float(t[3]), float(t[5]))
CYAN = "#19D3F3"
y = row_frame(6, CYAN, "Transfer between cell lines", ["Same knockdown measured in", "other cell lines, as a virtual", "cell model could learn it"], "Fig. 9")
for i, col in enumerate((BLUE, ORANGE, GREEN, PURPLE)): cell(C2 + 25 + (i % 2) * 45, y + 100 + (i // 2) * 50, col)
text(C2 + 47, y + 215, "other lines", 22, anchor="middle")
line(C2 + 100, y + 125, C2 + 170, y + 125, GREY, 4, head=True)
random.seed(5)
line(C2 + 180, y + 130, C2 + 305, y + 130, GREY, 1.5)
for j in range(9):
    v = (random.random() - 0.45) * 90; rect(C2 + 185 + j * 13, y + 130 - max(v, 0), 10, abs(v), ORANGE, 0.8)
text(C2 + 243, y + 215, "transferred", 22, anchor="middle")
text(C2 + 330, y + 140, "+", 40, "bold", "middle")
grid(C2 + 360, y + 85, 4, 7, 22, 22, {1, 4}, BLUE, missing=True)
text(C2 + 437, y + 215, "panel in target", 22, anchor="middle")
bx, by, bh = C3 + 40, y + 205, 130
axes(bx, by, 330, bh)
text(bx + 12, by - bh - 8, "RPE1, genes outside panel", 20)
for i, (lab, idx) in enumerate((("specific r", 0), ("wrong sign", 1))):
    x = bx + 20 + i * 165
    for k2, (m, col) in enumerate((("copy_pooled", ORANGE), ("panel", BLUE), ("panel+copy_pooled", GREEN))):
        v = VC[m][idx]; rect(x + k2 * 45, by - v / 0.22 * bh, 40, v / 0.22 * bh, col, 0.85)
    text(x + 65, by + 30, lab, 22, anchor="middle")
text(bx + 346, by - 115, "transfer", 24, fill=ORANGE); text(bx + 346, by - 85, "panel", 24, fill=BLUE); text(bx + 346, by - 55, "both", 24, fill=GREEN)
y_end = Y0 + 7 * RH
for x in (C2 - 30, C3 - 30): line(x, 30, x, y_end - 10, LINE, 3)
add("</svg>")
Hf = y_end + 10; out = "\n".join(svg).replace("@HMM@", f"{Hf / 10:.1f}").replace("@H@", str(Hf))
p = os.path.join(F, "Figure1.svg"); open(p, "w").write(out + "\n")
subprocess.run(["rsvg-convert", "-w", str(round(180 / 25.4 * 600)), "-b", "white", "-o", os.path.join(F, "Figure1.png"), p], check=True)
subprocess.run(["rsvg-convert", "-f", "pdf", "-o", os.path.join(F, "Figure1.pdf"), p], check=True); print("wrote", p)
