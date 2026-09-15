#!/usr/bin/env python3

import argparse
import gzip
import html
import json
from datetime import datetime
from pathlib import Path

VERSION = "0.2.0"
PIPELINE_NAME = "modular-human-genomics"
AUTHOR = "Hasan Efe Mihci"


def open_text(path):
    path = str(path)
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def parse_fai(path):
    contigs = {}
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                contigs[parts[0]] = int(parts[1])
    return contigs


def parse_format(fmt, sample):
    if not fmt or not sample:
        return {}
    keys = fmt.split(":")
    values = sample.split(":")
    return dict(zip(keys, values))


def variant_type(ref, alt):
    if len(ref) == 1 and len(alt) == 1:
        return "SNP"
    if len(ref) != len(alt):
        return "INDEL"
    return "OTHER"


def parse_vcf(path):
    variants = []
    sample_name = None

    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 10:
                    sample_name = fields[9]
                continue

            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue

            chrom = fields[0]
            pos = int(fields[1])
            vid = fields[2]
            ref = fields[3]
            alt_field = fields[4]
            qual = fields[5]
            filt = fields[6]

            fmt = fields[8] if len(fields) > 8 else ""
            sample = fields[9] if len(fields) > 9 else ""
            sample_data = parse_format(fmt, sample)

            for alt in alt_field.split(","):
                variants.append(
                    {
                        "chrom": chrom,
                        "pos": pos,
                        "id": vid,
                        "ref": ref,
                        "alt": alt,
                        "type": variant_type(ref, alt),
                        "qual": qual,
                        "filter": filt,
                        "gt": sample_data.get("GT", "."),
                        "dp": sample_data.get("DP", "."),
                        "gq": sample_data.get("GQ", "."),
                    }
                )

    return sample_name, variants


def synthetic_variants():
    return [
        {
            "chrom": "chr20",
            "pos": 10234881,
            "id": ".",
            "ref": "A",
            "alt": "G",
            "type": "SNP",
            "qual": "68.2",
            "filter": "PASS",
            "gt": "0/1",
            "dp": "31",
            "gq": "99",
        },
        {
            "chrom": "chr20",
            "pos": 18177234,
            "id": ".",
            "ref": "C",
            "alt": "T",
            "type": "SNP",
            "qual": "54.7",
            "filter": "PASS",
            "gt": "1/1",
            "dp": "27",
            "gq": "91",
        },
        {
            "chrom": "chr20",
            "pos": 29881123,
            "id": ".",
            "ref": "G",
            "alt": "GA",
            "type": "INDEL",
            "qual": "42.1",
            "filter": "PASS",
            "gt": "0/1",
            "dp": "24",
            "gq": "76",
        },
        {
            "chrom": "chr7",
            "pos": 55249071,
            "id": ".",
            "ref": "C",
            "alt": "T",
            "type": "SNP",
            "qual": "71.8",
            "filter": "PASS",
            "gt": "0/1",
            "dp": "36",
            "gq": "99",
        },
        {
            "chrom": "chr17",
            "pos": 43071077,
            "id": ".",
            "ref": "A",
            "alt": "G",
            "type": "SNP",
            "qual": "39.5",
            "filter": "PASS",
            "gt": "0/1",
            "dp": "21",
            "gq": "63",
        },
    ]


def ti_tv(variants):
    transitions = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
    ti = 0
    tv = 0

    for v in variants:
        if v["type"] != "SNP":
            continue

        pair = (v["ref"].upper(), v["alt"].upper())
        if pair in transitions:
            ti += 1
        else:
            tv += 1

    return ti, tv


def qc_status(path):
    if not path or not Path(path).exists():
        return "UNKNOWN"

    try:
        data = json.loads(Path(path).read_text(encoding="utf-8"))
    except Exception:
        return "UNKNOWN"

    for key in ("gate", "overall_status", "status", "qc_status"):
        if key in data:
            return str(data[key]).upper()

    return "UNKNOWN"


def build_html(
    sample_name,
    variants,
    contigs,
    demo=False,
    raw_qc_status="UNKNOWN",
    alignment_qc_status="UNKNOWN",
):
    canonical = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    canonical = [c for c in canonical if c in contigs]

    total = len(variants)
    snps = sum(v["type"] == "SNP" for v in variants)
    indels = sum(v["type"] == "INDEL" for v in variants)

    het = sum(v["gt"] in ("0/1", "1/0", "0|1", "1|0") for v in variants)
    hom_alt = sum(v["gt"] in ("1/1", "1|1") for v in variants)

    ti, tv = ti_tv(variants)
    ratio = round(ti / tv, 2) if tv else None

    metadata = {
        "sample": sample_name or "Unknown sample",
        "date": datetime.now().strftime("%Y-%m-%d %H:%M"),
        "reference": "GRCh38 no-alt analysis set",
        "pipeline_version": VERSION,
        "author": AUTHOR,
    }

    data_json = json.dumps(variants, separators=(",", ":"))
    contig_json = json.dumps({c: contigs[c] for c in canonical}, separators=(",", ":"))
    metadata_json = json.dumps(metadata)

    demo_label = "SYNTHETIC DEMO" if demo else "REAL VCF"

    raw_qc = html.escape(raw_qc_status)
    aln_qc = html.escape(alignment_qc_status)

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">

<title>Genomics Platform · Variant Browser</title>

<style>
:root {{
    --bg:#050914;
    --panel:#08111f;
    --panel2:#0c1728;
    --panel3:#0f1d31;
    --line:#18304f;
    --line2:#24496f;
    --text:#f2f7ff;
    --muted:#8ea5c3;
    --blue:#35a7ff;
    --blue2:#1c7de8;
    --green:#48d597;
    --orange:#ff934d;
    --yellow:#ffd166;
    --red:#ff5d73;
}}

* {{ box-sizing:border-box; }}

body {{
    margin:0;
    background:linear-gradient(180deg,#040914,#06101b 60%,#030811);
    color:var(--text);
    font-family:Inter,system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;
}}

.app {{
    display:grid;
    grid-template-columns:205px minmax(0,1fr) 330px;
    min-height:100vh;
}}

.sidebar {{
    border-right:1px solid var(--line);
    background:#07101d;
    padding:18px 10px;
    display:flex;
    flex-direction:column;
}}

.brand {{
    padding:5px 12px 22px;
}}

.brand-title {{
    font-size:19px;
    font-weight:800;
}}

.brand-sub {{
    color:#62b7ff;
    margin-top:4px;
    font-size:13px;
}}

.nav {{
    display:grid;
    gap:5px;
}}

.nav-item {{
    padding:11px 12px;
    border-radius:8px;
    color:#b8cae0;
    font-size:14px;
}}

.nav-item.active {{
    background:#0c4f98;
    color:white;
}}

.qc-box {{
    margin-top:auto;
    background:#08182a;
    border:1px solid var(--line2);
    border-radius:10px;
    padding:14px;
}}

.qc-title {{
    font-weight:700;
    margin-bottom:10px;
    color:#6fc2ff;
}}

.qc-row {{
    display:flex;
    justify-content:space-between;
    gap:8px;
    font-size:12px;
    margin:9px 0;
}}

.pass {{ color:var(--green); }}
.warn {{ color:var(--yellow); }}
.fail {{ color:var(--red); }}
.unknown {{ color:var(--muted); }}

.main {{
    padding:14px 16px 24px;
    min-width:0;
}}

.topbar {{
    display:grid;
    grid-template-columns:repeat(4,minmax(0,1fr));
    background:#07111f;
    border:1px solid var(--line);
    border-radius:10px;
    margin-bottom:12px;
}}

.meta {{
    padding:12px 16px;
    border-right:1px solid var(--line);
}}

.meta:last-child {{
    border-right:none;
}}

.meta-label {{
    color:#8fb0d3;
    font-size:11px;
}}

.meta-value {{
    font-size:13px;
    margin-top:4px;
    font-weight:600;
}}

.demo-banner {{
    margin-bottom:12px;
    padding:8px 12px;
    border:1px solid #735b19;
    color:#ffd166;
    background:#2a2209;
    border-radius:8px;
    font-size:12px;
}}

.kpis {{
    display:grid;
    grid-template-columns:repeat(8,minmax(0,1fr));
    gap:9px;
    margin-bottom:12px;
}}

.kpi {{
    min-height:86px;
    padding:12px;
    border:1px solid var(--line2);
    background:linear-gradient(180deg,#091629,#07111f);
    border-radius:9px;
}}

.kpi-label {{
    color:#b6c7dc;
    font-size:11px;
}}

.kpi-value {{
    font-size:23px;
    font-weight:800;
    margin-top:7px;
}}

.kpi-sub {{
    color:#7fa0c3;
    font-size:10px;
    margin-top:3px;
}}

.panel {{
    background:linear-gradient(180deg,#071321,#06101c);
    border:1px solid var(--line2);
    border-radius:10px;
    margin-bottom:12px;
    overflow:hidden;
}}

.panel-head {{
    padding:12px 14px;
    font-weight:800;
    border-bottom:1px solid var(--line);
}}

.genome {{
    display:flex;
    align-items:flex-end;
    gap:11px;
    padding:18px 14px 14px;
    overflow-x:auto;
}}

.chr {{
    cursor:pointer;
    text-align:center;
    min-width:28px;
}}

.ideogram {{
    width:17px;
    margin:auto;
    border:1px solid #a9b5c4;
    border-radius:7px;
    background:
      repeating-linear-gradient(
        to bottom,
        #edf1f5 0 8px,
        #777 8px 13px,
        #dfe4e9 13px 20px,
        #444 20px 25px
      );
    position:relative;
    box-shadow:inset 0 0 6px #000;
}}

.ideogram::after {{
    content:"";
    position:absolute;
    top:47%;
    left:-2px;
    right:-2px;
    height:7px;
    background:#071321;
    border-radius:50%;
}}

.chr.selected .ideogram {{
    outline:2px solid var(--blue);
    box-shadow:0 0 12px #1e88ff;
}}

.chr-label {{
    margin-top:8px;
    font-size:11px;
}}

.chromosome-view {{
    padding:14px;
}}

.big-ideogram {{
    position:relative;
    height:28px;
    border:1px solid #9ca8b7;
    border-radius:14px;
    overflow:hidden;
    background:
      repeating-linear-gradient(
        to right,
        #eef2f4 0 7%,
        #707780 7% 13%,
        #d7dce1 13% 20%,
        #4d5258 20% 27%
      );
}}

.centromere {{
    position:absolute;
    left:48%;
    top:-2px;
    width:32px;
    height:32px;
    transform:translateX(-50%) rotate(45deg);
    background:#be3f42;
}}

.track {{
    position:relative;
    height:68px;
    margin-top:18px;
    border-bottom:1px solid var(--line2);
}}

.track-line {{
    position:absolute;
    left:0;
    right:0;
    top:35px;
    height:3px;
    background:#2d83cb;
}}

.variant-marker {{
    position:absolute;
    top:23px;
    width:4px;
    height:25px;
    background:var(--yellow);
    cursor:pointer;
}}

.variant-marker.snp {{
    background:var(--green);
}}

.variant-marker.indel {{
    background:var(--orange);
}}

.axis {{
    display:flex;
    justify-content:space-between;
    color:var(--muted);
    font-size:10px;
    margin-top:6px;
}}

.toolbar {{
    display:flex;
    gap:10px;
    padding:10px 14px;
    border-bottom:1px solid var(--line);
}}

input,select,button {{
    background:#0b1a2c;
    color:var(--text);
    border:1px solid var(--line2);
    border-radius:7px;
    padding:8px 10px;
}}

button {{
    cursor:pointer;
}}

table {{
    width:100%;
    border-collapse:collapse;
    font-size:12px;
}}

th {{
    color:#9db3ce;
    text-align:left;
    padding:9px 10px;
    border-bottom:1px solid var(--line);
}}

td {{
    padding:9px 10px;
    border-bottom:1px solid #11243a;
}}

tbody tr {{
    cursor:pointer;
}}

tbody tr:hover {{
    background:#0d2b4a;
}}

tbody tr.selected {{
    background:#0b55a5;
}}

.badge {{
    padding:3px 7px;
    border-radius:999px;
    font-size:10px;
    font-weight:700;
}}

.badge-pass {{
    background:#144d36;
    color:#7ef1bb;
}}

.badge-snp {{
    background:#103b65;
    color:#6ec5ff;
}}

.badge-indel {{
    background:#5e3717;
    color:#ffb780;
}}

.empty {{
    padding:45px 20px;
    text-align:center;
    color:var(--muted);
}}

.charts {{
    display:grid;
    grid-template-columns:repeat(3,minmax(0,1fr));
    gap:10px;
}}

.chart {{
    background:#071321;
    border:1px solid var(--line2);
    border-radius:10px;
    padding:12px;
    min-height:165px;
}}

.chart-title {{
    font-size:12px;
    font-weight:700;
    margin-bottom:12px;
}}

.donut-wrap {{
    display:flex;
    align-items:center;
    gap:18px;
}}

.donut {{
    width:88px;
    height:88px;
    border-radius:50%;
    background:conic-gradient(
        var(--blue) 0 60%,
        var(--orange) 60% 90%,
        var(--yellow) 90% 100%
    );
    position:relative;
}}

.donut::after {{
    content:"";
    position:absolute;
    inset:18px;
    background:#071321;
    border-radius:50%;
}}

.legend {{
    color:#b3c6db;
    font-size:11px;
    line-height:1.9;
}}

.rightbar {{
    border-left:1px solid var(--line);
    background:#06101b;
    padding:14px 12px;
}}

.detail-panel {{
    position:sticky;
    top:14px;
    background:#071321;
    border:1px solid var(--line2);
    border-radius:10px;
    overflow:hidden;
}}

.detail-head {{
    padding:13px;
    border-bottom:1px solid var(--line);
    font-weight:800;
    color:#6fc2ff;
}}

.detail-body {{
    padding:14px;
}}

.detail-title {{
    font-size:20px;
    font-weight:800;
    margin-bottom:8px;
}}

.section {{
    margin-top:22px;
}}

.section-title {{
    color:#67bdff;
    font-size:13px;
    font-weight:800;
    margin-bottom:8px;
}}

.detail-row {{
    display:flex;
    justify-content:space-between;
    gap:12px;
    padding:7px 0;
    border-bottom:1px solid #132741;
    font-size:12px;
}}

.detail-row span:first-child {{
    color:#9cb1ca;
}}

.note {{
    padding:10px;
    background:#0b1b2d;
    color:#9eb7d2;
    font-size:11px;
    border-radius:7px;
    line-height:1.5;
}}

@media(max-width:1300px) {{
    .app {{
        grid-template-columns:180px minmax(0,1fr);
    }}

    .rightbar {{
        grid-column:1/-1;
        border-left:none;
        border-top:1px solid var(--line);
    }}

    .kpis {{
        grid-template-columns:repeat(4,1fr);
    }}
}}

</style>
</head>

<body>

<div class="app">

<aside class="sidebar">

    <div class="brand">
        <div class="brand-title">🧬 Genomics Platform</div>
        <div class="brand-sub">Variant Browser · v{VERSION}</div>
    </div>

    <div class="nav">
        <div class="nav-item">⌂ Overview</div>
        <div class="nav-item active">▥ Genome Browser</div>
        <div class="nav-item">☷ Variant Table</div>
        <div class="nav-item">⌁ Variant QC</div>
        <div class="nav-item">⌬ Alignment QC</div>
        <div class="nav-item">▣ Raw QC</div>
        <div class="nav-item">▤ Reports</div>
        <div class="nav-item">⚙ Settings</div>
    </div>

    <div class="qc-box">
        <div class="qc-title">Pipeline QC Summary</div>

        <div class="qc-row">
            <span>Raw QC</span>
            <strong class="{raw_qc.lower() if raw_qc in ['PASS','WARN','FAIL'] else 'unknown'}">{raw_qc}</strong>
        </div>

        <div class="qc-row">
            <span>QC Gate</span>
            <strong class="{raw_qc.lower() if raw_qc in ['PASS','WARN','FAIL'] else 'unknown'}">{raw_qc}</strong>
        </div>

        <div class="qc-row">
            <span>Alignment</span>
            <strong class="pass">COMPLETED</strong>
        </div>

        <div class="qc-row">
            <span>Alignment QC</span>
            <strong class="{aln_qc.lower() if aln_qc in ['PASS','WARN','FAIL'] else 'unknown'}">{aln_qc}</strong>
        </div>

        <div class="qc-row">
            <span>Variant Calling</span>
            <strong class="pass">COMPLETED</strong>
        </div>
    </div>

</aside>


<main class="main">

    <div class="topbar">
        <div class="meta">
            <div class="meta-label">Sample ID</div>
            <div class="meta-value" id="metaSample"></div>
        </div>

        <div class="meta">
            <div class="meta-label">Analysis Date</div>
            <div class="meta-value" id="metaDate"></div>
        </div>

        <div class="meta">
            <div class="meta-label">Reference Genome</div>
            <div class="meta-value" id="metaReference"></div>
        </div>

        <div class="meta">
            <div class="meta-label">Pipeline / Author</div>
            <div class="meta-value" id="metaPipeline"></div>
        </div>
    </div>

    {"<div class='demo-banner'>SYNTHETIC DEMO DATA — not biological results.</div>" if demo else ""}

    <div class="kpis">
        <div class="kpi">
            <div class="kpi-label">Total Variants</div>
            <div class="kpi-value">{total}</div>
            <div class="kpi-sub">across loaded VCF</div>
        </div>

        <div class="kpi">
            <div class="kpi-label">SNPs</div>
            <div class="kpi-value">{snps}</div>
        </div>

        <div class="kpi">
            <div class="kpi-label">INDELs</div>
            <div class="kpi-value">{indels}</div>
        </div>

        <div class="kpi">
            <div class="kpi-label">Transitions</div>
            <div class="kpi-value">{ti}</div>
        </div>

        <div class="kpi">
            <div class="kpi-label">Transversions</div>
            <div class="kpi-value">{tv}</div>
        </div>

        <div class="kpi">
            <div class="kpi-label">Ti/Tv</div>
            <div class="kpi-value">{ratio if ratio is not None else "-"}</div>
        </div>

        <div class="kpi">
            <div class="kpi-label">Heterozygous</div>
            <div class="kpi-value">{het}</div>
        </div>

        <div class="kpi">
            <div class="kpi-label">Homozygous ALT</div>
            <div class="kpi-value">{hom_alt}</div>
        </div>
    </div>

    <div class="panel">
        <div class="panel-head">
            Human Genome (GRCh38)
        </div>
        <div id="genome" class="genome"></div>
    </div>

    <div class="panel">
        <div class="panel-head" id="chromTitle">
            Chromosome
        </div>

        <div class="chromosome-view">
            <div class="big-ideogram">
                <div class="centromere"></div>
            </div>

            <div id="variantTrack" class="track">
                <div class="track-line"></div>
            </div>

            <div class="axis">
                <span>0 Mb</span>
                <span id="chromEnd"></span>
            </div>
        </div>
    </div>

    <div class="panel">

        <div class="toolbar">
            <select id="typeFilter">
                <option value="ALL">All variants</option>
                <option value="SNP">SNP</option>
                <option value="INDEL">INDEL</option>
                <option value="OTHER">Other</option>
            </select>

            <input
                id="searchBox"
                type="text"
                placeholder="Search position / REF / ALT">

            <button id="downloadCsv">
                Download CSV
            </button>
        </div>

        <div id="variantTable"></div>

    </div>

    <div class="charts">

        <div class="chart">
            <div class="chart-title">
                Variant Type Distribution
            </div>

            <div class="donut-wrap">
                <div class="donut"></div>

                <div class="legend">
                    SNP · {snps}<br>
                    INDEL · {indels}<br>
                    Other · {total-snps-indels}
                </div>
            </div>
        </div>

        <div class="chart">
            <div class="chart-title">
                Genotype Distribution
            </div>

            <div class="donut-wrap">
                <div class="donut"></div>

                <div class="legend">
                    Heterozygous · {het}<br>
                    Homozygous ALT · {hom_alt}<br>
                    Other · {max(total-het-hom_alt,0)}
                </div>
            </div>
        </div>

        <div class="chart">
            <div class="chart-title">
                Quality Distribution
            </div>

            <div class="note">
                QUAL histogram will be upgraded in Variant QC v0.3.
                Current browser exposes raw QUAL values per variant.
            </div>
        </div>

    </div>

</main>


<aside class="rightbar">

    <div class="detail-panel">

        <div class="detail-head">
            Variant Details
        </div>

        <div id="details" class="detail-body">
            <div class="note">
                Select a variant from the chromosome track or variant table.
            </div>
        </div>

    </div>

</aside>

</div>


<script>

const variants = {data_json};
const contigs = {contig_json};
const metadata = {metadata_json};

const chromosomes = Object.keys(contigs);

let selectedChrom =
    chromosomes.includes("chr20")
        ? "chr20"
        : chromosomes[0];

let selectedVariant = null;


document.getElementById("metaSample").textContent =
    metadata.sample;

document.getElementById("metaDate").textContent =
    metadata.date;

document.getElementById("metaReference").textContent =
    metadata.reference;

document.getElementById("metaPipeline").textContent =
    `${{metadata.pipeline_version}} · ${{metadata.author}}`;


function chrVariants(chrom) {{
    return variants.filter(v => v.chrom === chrom);
}}


function renderGenome() {{

    const root =
        document.getElementById("genome");

    root.innerHTML = "";

    const maxLen =
        Math.max(...Object.values(contigs));

    chromosomes.forEach(chrom => {{

        const box =
            document.createElement("div");

        box.className =
            "chr" +
            (chrom === selectedChrom
                ? " selected"
                : "");

        const id =
            document.createElement("div");

        id.className =
            "ideogram";

        id.style.height =
            `${{45 + (contigs[chrom] / maxLen) * 55}}px`;

        const label =
            document.createElement("div");

        label.className =
            "chr-label";

        label.textContent =
            chrom.replace("chr","");

        box.appendChild(id);
        box.appendChild(label);

        box.onclick = () => {{
            selectedChrom = chrom;
            selectedVariant = null;
            renderGenome();
            renderChromosome();
            renderTable();
            clearDetails();
        }};

        root.appendChild(box);
    }});
}}


function renderChromosome() {{

    const length = contigs[selectedChrom];
    const vars = chrVariants(selectedChrom);

    document.getElementById("chromTitle").textContent =
        `${{selectedChrom.replace("chr","Chromosome ")}} · ${{length.toLocaleString()}} bp`;

    document.getElementById("chromEnd").textContent =
        `${{(length / 1e6).toFixed(1)}} Mb`;

    const track =
        document.getElementById("variantTrack");

    track.innerHTML =
        '<div class="track-line"></div>';

    vars.forEach(v => {{

        const marker =
            document.createElement("div");

        marker.className =
            "variant-marker " +
            v.type.toLowerCase();

        marker.style.left =
            `${{(v.pos / length) * 100}}%`;

        marker.title =
            `${{v.chrom}}:${{v.pos}} ${{v.ref}}>${{v.alt}}`;

        marker.onclick = () =>
            selectVariant(v);

        track.appendChild(marker);
    }});
}}


function filteredVariants() {{

    const type =
        document.getElementById("typeFilter").value;

    const q =
        document.getElementById("searchBox")
        .value
        .toLowerCase();

    return chrVariants(selectedChrom)
        .filter(v => {{

            const typeOk =
                type === "ALL" ||
                v.type === type;

            const haystack =
                `${{v.pos}} ${{v.ref}} ${{v.alt}}`
                .toLowerCase();

            const searchOk =
                !q || haystack.includes(q);

            return typeOk && searchOk;
        }});
}}


function renderTable() {{

    const rows = filteredVariants();

    const root =
        document.getElementById("variantTable");

    if (!rows.length) {{

        root.innerHTML = `
            <div class="empty">
                No variants found in ${{selectedChrom}}
            </div>
        `;

        return;
    }}

    let out = `
        <table>
        <thead>
        <tr>
            <th>Position</th>
            <th>REF</th>
            <th>ALT</th>
            <th>Type</th>
            <th>Genotype</th>
            <th>DP</th>
            <th>GQ</th>
            <th>QUAL</th>
            <th>Filter</th>
        </tr>
        </thead>
        <tbody>
    `;

    rows.forEach((v,i) => {{

        const selected =
            selectedVariant === v
                ? "selected"
                : "";

        out += `
            <tr class="${{selected}}"
                data-row="${{i}}">
                <td>${{v.pos.toLocaleString()}}</td>
                <td>${{v.ref}}</td>
                <td>${{v.alt}}</td>
                <td>
                    <span class="badge badge-${{v.type.toLowerCase()}}">
                        ${{v.type}}
                    </span>
                </td>
                <td>${{v.gt}}</td>
                <td>${{v.dp}}</td>
                <td>${{v.gq}}</td>
                <td>${{v.qual}}</td>
                <td>
                    <span class="badge badge-pass">
                        ${{v.filter}}
                    </span>
                </td>
            </tr>
        `;
    }});

    out += "</tbody></table>";

    root.innerHTML = out;

    root.querySelectorAll("tbody tr")
        .forEach((tr,i) => {{
            tr.onclick = () =>
                selectVariant(rows[i]);
        }});
}}


function clearDetails() {{
    document.getElementById("details").innerHTML = `
        <div class="note">
            Select a variant from the chromosome track or table.
        </div>
    `;
}}


function selectVariant(v) {{

    selectedVariant = v;

    renderTable();

    document.getElementById("details").innerHTML = `

        <div class="detail-title">
            ${{v.chrom}}:${{v.pos.toLocaleString()}}
        </div>

        <div style="font-size:18px;margin-bottom:10px">
            ${{v.ref}} → ${{v.alt}}
        </div>

        <div>
            <span class="badge badge-${{v.type.toLowerCase()}}">
                ${{v.type}}
            </span>

            <span class="badge badge-pass">
                ${{v.filter}}
            </span>
        </div>


        <div class="section">

            <div class="section-title">
                Technical Information
            </div>

            <div class="detail-row">
                <span>Chromosome</span>
                <strong>${{v.chrom}}</strong>
            </div>

            <div class="detail-row">
                <span>Position</span>
                <strong>${{v.pos.toLocaleString()}}</strong>
            </div>

            <div class="detail-row">
                <span>Reference</span>
                <strong>${{v.ref}}</strong>
            </div>

            <div class="detail-row">
                <span>Alternative</span>
                <strong>${{v.alt}}</strong>
            </div>

            <div class="detail-row">
                <span>Variant Type</span>
                <strong>${{v.type}}</strong>
            </div>

            <div class="detail-row">
                <span>Genotype</span>
                <strong>${{v.gt}}</strong>
            </div>

            <div class="detail-row">
                <span>Depth (DP)</span>
                <strong>${{v.dp}}</strong>
            </div>

            <div class="detail-row">
                <span>Genotype Quality</span>
                <strong>${{v.gq}}</strong>
            </div>

            <div class="detail-row">
                <span>Variant Quality</span>
                <strong>${{v.qual}}</strong>
            </div>

            <div class="detail-row">
                <span>Filter</span>
                <strong>${{v.filter}}</strong>
            </div>

        </div>


        <div class="section">

            <div class="section-title">
                Gene & Consequence
            </div>

            <div class="note">
                Annotation engine not connected yet.
                Gene, transcript, consequence and protein change
                will appear here after VEP/SnpEff integration.
            </div>

        </div>


        <div class="section">

            <div class="section-title">
                Population Frequency
            </div>

            <div class="note">
                gnomAD population evidence not connected yet.
            </div>

        </div>


        <div class="section">

            <div class="section-title">
                Clinical Evidence
            </div>

            <div class="note">
                ClinVar, ClinGen, phenotype, inheritance and
                evidence classification layers are intentionally
                not connected yet.
            </div>

        </div>


        <div class="section">

            <div class="section-title">
                Interpretation
            </div>

            <div class="note">
                Technical variant call available.
                Clinical interpretation remains unresolved until
                annotation and evidence layers are connected.
            </div>

        </div>
    `;
}}


function downloadCsv() {{

    const rows = filteredVariants();

    let csv =
        "chrom,pos,ref,alt,type,gt,dp,gq,qual,filter\\n";

    rows.forEach(v => {{
        csv += [
            v.chrom,
            v.pos,
            v.ref,
            v.alt,
            v.type,
            v.gt,
            v.dp,
            v.gq,
            v.qual,
            v.filter
        ].join(",") + "\\n";
    }});

    const blob =
        new Blob([csv], {{
            type:"text/csv"
        }});

    const url =
        URL.createObjectURL(blob);

    const a =
        document.createElement("a");

    a.href = url;

    a.download =
        `${{selectedChrom}}_variants.csv`;

    a.click();

    URL.revokeObjectURL(url);
}}


document.getElementById("typeFilter")
    .addEventListener(
        "change",
        renderTable
    );

document.getElementById("searchBox")
    .addEventListener(
        "input",
        renderTable
    );

document.getElementById("downloadCsv")
    .addEventListener(
        "click",
        downloadCsv
    );


renderGenome();
renderChromosome();
renderTable();

</script>

</body>
</html>
"""


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--vcf")
    parser.add_argument("--fai", required=True)
    parser.add_argument("--output", required=True)

    parser.add_argument("--raw-qc-json")
    parser.add_argument("--alignment-qc-json")

    parser.add_argument(
        "--demo",
        action="store_true"
    )

    args = parser.parse_args()

    contigs = parse_fai(args.fai)

    if args.demo:
        sample_name = "HG002_BROWSER_DEMO"
        variants = synthetic_variants()
    else:
        if not args.vcf:
            parser.error(
                "--vcf is required unless --demo is used"
            )

        sample_name, variants = parse_vcf(args.vcf)

    raw_status = qc_status(args.raw_qc_json)
    alignment_status = qc_status(args.alignment_qc_json)

    page = build_html(
        sample_name=sample_name,
        variants=variants,
        contigs=contigs,
        demo=args.demo,
        raw_qc_status=raw_status,
        alignment_qc_status=alignment_status,
    )

    output = Path(args.output)

    output.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    output.write_text(
        page,
        encoding="utf-8"
    )

    print(
        f"Variant browser written to: {output}"
    )

    print(
        f"Variants loaded: {len(variants)}"
    )


if __name__ == "__main__":
    main()
