#!/usr/bin/env python3

import argparse
import html
import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path

import yaml


QC_ENGINE_VERSION = "0.2.0"
REPORT_SCHEMA_VERSION = "0.2.0"
AUTHOR = "Efe Mıhcı"


# ============================================================
# Generic helpers
# ============================================================

def read_text(path):
    return Path(path).read_text(errors="replace")


def safe_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def fmt(value, digits=2):
    if value is None:
        return "N/A"

    if isinstance(value, int):
        return f"{value:,}"

    try:
        value = float(value)
    except Exception:
        return str(value)

    if math.isnan(value):
        return "N/A"

    return f"{value:.{digits}f}"


def pct(value):
    return "N/A" if value is None else f"{value:.2f}%"


# ============================================================
# samtools flagstat
# ============================================================

def parse_flagstat_count(text, label):
    pattern = rf"^(\d+)\s+\+\s+(\d+)\s+{re.escape(label)}(?:\s|\(|$)"

    match = re.search(
        pattern,
        text,
        flags=re.MULTILINE,
    )

    if not match:
        return None

    return int(match.group(1)) + int(match.group(2))


def parse_flagstat(path):
    text = read_text(path)

    metrics = {}

    labels = {
        "total_reads": "in total",
        "primary_reads": "primary",
        "secondary_reads": "secondary",
        "supplementary_reads": "supplementary",
        "duplicates": "duplicates",
        "primary_duplicates": "primary duplicates",
        "mapped_reads": "mapped",
        "primary_mapped_reads": "primary mapped",
        "paired_in_sequencing": "paired in sequencing",
        "properly_paired_reads": "properly paired",
        "singletons": "singletons",
    }

    for key, label in labels.items():
        value = parse_flagstat_count(text, label)

        if value is not None:
            metrics[key] = value

    primary = metrics.get("primary_reads")
    primary_mapped = metrics.get("primary_mapped_reads")

    if (
        primary is not None
        and primary_mapped is not None
        and primary > 0
    ):
        metrics["mapping_rate_pct"] = (
            primary_mapped / primary * 100.0
        )

        metrics["primary_unmapped_reads"] = (
            primary - primary_mapped
        )

    elif metrics.get("total_reads") and metrics.get("mapped_reads") is not None:
        metrics["mapping_rate_pct"] = (
            metrics["mapped_reads"]
            / metrics["total_reads"]
            * 100.0
        )

    paired = metrics.get("paired_in_sequencing")
    properly = metrics.get("properly_paired_reads")

    if paired and properly is not None:
        metrics["properly_paired_pct"] = (
            properly / paired * 100.0
        )

    duplicate_numerator = metrics.get(
        "primary_duplicates",
        metrics.get("duplicates"),
    )

    duplicate_denominator = metrics.get(
        "primary_reads",
        metrics.get("total_reads"),
    )

    if (
        duplicate_numerator is not None
        and duplicate_denominator
    ):
        metrics["duplicate_rate_pct"] = (
            duplicate_numerator
            / duplicate_denominator
            * 100.0
        )

    return metrics


# ============================================================
# samtools stats
# ============================================================

def parse_samtools_stats(path):
    stats = {}

    for line in read_text(path).splitlines():

        if not line.startswith("SN\t"):
            continue

        fields = line.split("\t")

        if len(fields) < 3:
            continue

        key = fields[1].rstrip(":").strip()
        value = fields[2].strip()

        number = safe_float(value)

        stats[key] = number if number is not None else value

    return stats


# ============================================================
# idxstats
# ============================================================

def parse_idxstats(path):
    rows = []

    for line in read_text(path).splitlines():

        fields = line.split("\t")

        if len(fields) != 4:
            continue

        contig, length, mapped, unmapped = fields

        rows.append(
            {
                "contig": contig,
                "length": int(length),
                "mapped": int(mapped),
                "unmapped": int(unmapped),
            }
        )

    return rows


# ============================================================
# mosdepth summary
# ============================================================

def parse_mosdepth_summary(path):
    autosomes = {f"chr{i}" for i in range(1, 23)}

    rows = []

    weighted_autosomal_bases = 0.0
    autosomal_length = 0

    for line in read_text(path).splitlines():

        fields = line.split("\t")

        if not fields:
            continue

        if fields[0] == "chrom":
            continue

        if len(fields) < 4:
            continue

        chrom = fields[0]

        try:
            length = int(fields[1])
            mean = float(fields[3])
        except ValueError:
            continue

        rows.append(
            {
                "contig": chrom,
                "length": length,
                "mean_coverage": mean,
            }
        )

        if chrom in autosomes:
            weighted_autosomal_bases += (
                length * mean
            )

            autosomal_length += length

    mean_autosomal = None

    if autosomal_length:
        mean_autosomal = (
            weighted_autosomal_bases
            / autosomal_length
        )

    return rows, mean_autosomal


# ============================================================
# mosdepth cumulative global distribution
# ============================================================

def parse_global_dist(path):
    """
    mosdepth *.global.dist.txt is cumulative:
        chrom   depth   fraction_of_bases_at_or_above_depth

    We return percentage of total genome >= 1x, 10x, 20x, 30x.
    """

    total_distribution = {}

    for line in read_text(path).splitlines():

        fields = line.split("\t")

        if len(fields) != 3:
            continue

        chrom, depth, fraction = fields

        if chrom != "total":
            continue

        try:
            depth = int(depth)
            fraction = float(fraction)
        except ValueError:
            continue

        total_distribution[depth] = fraction

    thresholds = [1, 10, 20, 30]

    result = {}

    if not total_distribution:
        return {
            threshold: None
            for threshold in thresholds
        }

    observed_depths = sorted(
        total_distribution
    )

    max_depth = max(
        observed_depths
    )

    for threshold in thresholds:

        if threshold in total_distribution:
            fraction = total_distribution[
                threshold
            ]

        elif threshold > max_depth:
            fraction = 0.0

        else:
            higher_depths = [
                depth
                for depth in observed_depths
                if depth > threshold
            ]

            if higher_depths:
                fraction = total_distribution[
                    min(higher_depths)
                ]
            else:
                fraction = 0.0

        result[threshold] = (
            fraction * 100.0
        )

    return result


# ============================================================
# QC logic
# ============================================================

def evaluate_lower(metric, value, thresholds):
    rule = thresholds.get(metric)

    if not rule:
        return "PASS", None

    if value is None:
        return "FAIL", "Metric unavailable"

    fail = rule.get("fail_below")
    warn = rule.get("warn_below")

    if fail is not None and value < float(fail):
        return (
            "FAIL",
            f"{value:.2f} < fail threshold {fail}",
        )

    if warn is not None and value < float(warn):
        return (
            "WARN",
            f"{value:.2f} < warning threshold {warn}",
        )

    return "PASS", None


def evaluate_upper(metric, value, thresholds):
    rule = thresholds.get(metric)

    if not rule:
        return "PASS", None

    if value is None:
        return "FAIL", "Metric unavailable"

    fail = rule.get("fail_above")
    warn = rule.get("warn_above")

    if fail is not None and value > float(fail):
        return (
            "FAIL",
            f"{value:.2f} > fail threshold {fail}",
        )

    if warn is not None and value > float(warn):
        return (
            "WARN",
            f"{value:.2f} > warning threshold {warn}",
        )

    return "PASS", None


def combine_gate(statuses):
    statuses = [
        value
        for value in statuses
        if value not in {
            "N/A",
            "NOT_APPLICABLE",
        }
    ]

    if "FAIL" in statuses:
        return "FAIL"

    if "WARN" in statuses:
        return "WARN"

    return "PASS"


# ============================================================
# HTML helpers
# ============================================================

def status_class(status):
    return {
        "PASS": "pass",
        "WARN": "warn",
        "FAIL": "fail",
        "N/A": "na",
        "NOT_APPLICABLE": "na",
    }.get(status, "na")


def metric_card(
    title,
    value,
    suffix,
    status,
    subtitle,
):
    return f"""
    <div class="metric-card">
      <div class="metric-top">
        <span class="metric-label">{html.escape(title)}</span>
        <span class="mini-status {status_class(status)}">
          {html.escape(status)}
        </span>
      </div>

      <div class="metric-value">
        {html.escape(str(value))}
        <span>{html.escape(suffix)}</span>
      </div>

      <div class="metric-subtitle">
        {html.escape(subtitle)}
      </div>
    </div>
    """


def donut_svg(
    pass_count,
    warn_count,
    fail_count,
):
    total = (
        pass_count
        + warn_count
        + fail_count
    )

    if total == 0:
        pass_pct = warn_pct = fail_pct = 0
    else:
        pass_pct = pass_count / total * 100
        warn_pct = warn_count / total * 100
        fail_pct = fail_count / total * 100

    first = pass_pct
    second = pass_pct + warn_pct

    return f"""
    <div
      class="donut"
      style="
        background:
        conic-gradient(
          var(--pass) 0% {first:.2f}%,
          var(--warn) {first:.2f}% {second:.2f}%,
          var(--fail) {second:.2f}% 100%
        );
      "
    >
      <div class="donut-inner">
        <div class="donut-number">{total}</div>
        <div class="donut-label">QC checks</div>
      </div>
    </div>
    """


def stacked_bar(segments):
    total = sum(
        max(value, 0)
        for _, value, _ in segments
    )

    pieces = []

    legend = []

    for label, value, css_class in segments:

        percentage = (
            value / total * 100.0
            if total
            else 0
        )

        pieces.append(
            f"""
            <div
              class="stack-segment {css_class}"
              style="width:{percentage:.4f}%"
              title="{html.escape(label)}: {value:,} ({percentage:.2f}%)">
            </div>
            """
        )

        legend.append(
            f"""
            <div class="legend-item">
              <span class="legend-dot {css_class}"></span>
              <span>{html.escape(label)}</span>
              <strong>{value:,}</strong>
              <small>{percentage:.2f}%</small>
            </div>
            """
        )

    return (
        '<div class="stack-bar">'
        + "".join(pieces)
        + "</div>"
        + '<div class="legend">'
        + "".join(legend)
        + "</div>"
    )


def chromosome_chart(coverage_rows):
    wanted = [
        *(f"chr{i}" for i in range(1, 23)),
        "chrX",
        "chrY",
        "chrM",
    ]

    coverage_map = {
        row["contig"]: row["mean_coverage"]
        for row in coverage_rows
    }

    values = [
        coverage_map.get(chrom, 0.0)
        for chrom in wanted
    ]

    maximum = max(values) if values else 0.0

    if maximum <= 0:
        maximum = 1.0

    bars = []

    for chrom in wanted:

        value = coverage_map.get(
            chrom,
            0.0,
        )

        width = (
            value / maximum * 100
            if maximum
            else 0
        )

        visual_width = (
            max(width, 0.8)
            if value > 0
            else 0
        )

        bars.append(
            f"""
            <div class="chrom-row">
              <div class="chrom-name">
                {html.escape(chrom)}
              </div>

              <div class="chrom-track">
                <div
                  class="chrom-bar"
                  style="width:{visual_width:.3f}%"
                  title="{html.escape(chrom)}: {value:.4f}×">
                </div>
              </div>

              <div class="chrom-value">
                {value:.3f}×
              </div>
            </div>
            """
        )

    return "".join(bars)


def threshold_chart(metrics):
    thresholds = [
        (
            "≥1×",
            metrics.get(
                "genome_1x_fraction_pct"
            ),
        ),
        (
            "≥10×",
            metrics.get(
                "genome_10x_fraction_pct"
            ),
        ),
        (
            "≥20×",
            metrics.get(
                "genome_20x_fraction_pct"
            ),
        ),
        (
            "≥30×",
            metrics.get(
                "genome_30x_fraction_pct"
            ),
        ),
    ]

    rows = []

    for label, value in thresholds:

        if value is None:
            width = 0
            value_text = "N/A"
        else:
            width = max(
                0,
                min(
                    100,
                    value,
                ),
            )

            value_text = (
                f"{value:.2f}%"
            )

        rows.append(
            f"""
            <div class="threshold-row">
              <div class="threshold-label">
                {label}
              </div>

              <div class="threshold-track">
                <div
                  class="threshold-bar"
                  style="width:{width:.3f}%"
                  title="{label}: {value_text}">
                </div>
              </div>

              <div class="threshold-value">
                {value_text}
              </div>
            </div>
            """
        )

    return "".join(rows)


# ============================================================
# HTML report
# ============================================================

def render_html(report):
    metrics = report["metrics"]
    events = report["events"]

    technical_gate = report[
        "technical_alignment_gate"
    ]

    coverage_gate = report[
        "coverage_gate"
    ]

    overall_gate = report[
        "gate"
    ]

    profile = report[
        "profile"
    ]

    coverage_applicable = report[
        "coverage_applicable"
    ]

    event_rows = ""

    pass_count = 0
    warn_count = 0
    fail_count = 0

    for event in events:

        status = event["status"]

        if status == "PASS":
            pass_count += 1

        elif status == "WARN":
            warn_count += 1

        elif status == "FAIL":
            fail_count += 1

        reason = (
            event.get("reason")
            or "Within configured threshold"
        )

        value = event.get("value")

        value_text = (
            fmt(value)
            if value is not None
            else "N/A"
        )

        event_rows += f"""
        <tr>
          <td>
            {html.escape(event["label"])}
          </td>

          <td>
            {html.escape(value_text)}
          </td>

          <td>
            <span class="table-status {status_class(status)}">
              {html.escape(status)}
            </span>
          </td>

          <td>
            {html.escape(reason)}
          </td>
        </tr>
        """

    coverage_rows = report[
        "coverage_by_contig"
    ]

    contig_rows = ""

    for row in coverage_rows:

        contig_rows += f"""
        <tr>
          <td>{html.escape(row["contig"])}</td>
          <td>{row["length"]:,}</td>
          <td>{row["mean_coverage"]:.4f}×</td>
        </tr>
        """

    mapping_status = next(
        (
            event["status"]
            for event in events
            if event["metric"]
            == "mapping_rate_pct"
        ),
        "N/A",
    )

    pairing_status = next(
        (
            event["status"]
            for event in events
            if event["metric"]
            == "properly_paired_pct"
        ),
        "N/A",
    )

    duplicate_status = next(
        (
            event["status"]
            for event in events
            if event["metric"]
            == "duplicate_rate_pct"
        ),
        "N/A",
    )

    coverage_status = (
        next(
            (
                event["status"]
                for event in events
                if event["metric"]
                == "mean_autosomal_coverage"
            ),
            "N/A",
        )
        if coverage_applicable
        else "N/A"
    )

    cards = "".join(
        [
            metric_card(
                "Primary mapping",
                fmt(
                    metrics.get(
                        "mapping_rate_pct"
                    )
                ),
                "%",
                mapping_status,
                "Primary reads mapped to reference",
            ),

            metric_card(
                "Properly paired",
                fmt(
                    metrics.get(
                        "properly_paired_pct"
                    )
                ),
                "%",
                pairing_status,
                "Expected paired-end orientation",
            ),

            metric_card(
                "Duplicate rate",
                fmt(
                    metrics.get(
                        "duplicate_rate_pct"
                    )
                ),
                "%",
                duplicate_status,
                "Marked primary duplicate reads",
            ),

            metric_card(
                "Mean autosomal",
                fmt(
                    metrics.get(
                        "mean_autosomal_coverage"
                    )
                ),
                "×",
                coverage_status,
                (
                    "Coverage gate active"
                    if coverage_applicable
                    else "Not evaluated in smoke mode"
                ),
            ),

            metric_card(
                "Genome ≥20×",
                fmt(
                    metrics.get(
                        "genome_20x_fraction_pct"
                    )
                ),
                "%",
                (
                    next(
                        (
                            event["status"]
                            for event in events
                            if event["metric"]
                            == "genome_20x_fraction_pct"
                        ),
                        "N/A",
                    )
                    if coverage_applicable
                    else "N/A"
                ),
                (
                    "Genome fraction above threshold"
                    if coverage_applicable
                    else "Not evaluated in smoke mode"
                ),
            ),

            metric_card(
                "Primary reads",
                fmt(
                    metrics.get(
                        "primary_reads"
                    )
                ),
                "",
                "PASS",
                "Reads evaluated for core alignment QC",
            ),
        ]
    )

    primary = metrics.get(
        "primary_reads",
        0,
    )

    secondary = metrics.get(
        "secondary_reads",
        0,
    )

    supplementary = metrics.get(
        "supplementary_reads",
        0,
    )

    mapped = metrics.get(
        "primary_mapped_reads",
        0,
    )

    unmapped = metrics.get(
        "primary_unmapped_reads",
        0,
    )

    duplicate_reads = metrics.get(
        "primary_duplicates",
        metrics.get(
            "duplicates",
            0,
        ),
    )

    unique_reads = max(
        primary - duplicate_reads,
        0,
    )

    smoke_notice = ""

    if not coverage_applicable:
        smoke_notice = """
        <div class="notice">
          <div class="notice-icon">i</div>
          <div>
            <strong>Smoke-test coverage mode</strong>
            <p>
              Genome-wide coverage sufficiency is intentionally not
              evaluated for this downsampled test dataset.
              Low coverage is expected and does not represent an
              alignment failure.
            </p>
          </div>
        </div>
        """

    technical_message = {
        "PASS": "Core alignment metrics are within the configured development thresholds.",
        "WARN": "One or more alignment metrics require review.",
        "FAIL": "One or more core alignment metrics failed configured thresholds.",
    }.get(
        technical_gate,
        "",
    )

    coverage_message = (
        {
            "PASS": "Genome-wide coverage is within configured development thresholds.",
            "WARN": "Coverage requires review.",
            "FAIL": "Coverage does not meet the configured development profile.",
        }.get(
            coverage_gate,
            ""
        )
        if coverage_applicable
        else "Coverage sufficiency is not applicable to this smoke-test profile."
    )

    generated = html.escape(
        report["generated_at"]
    )

    return f"""<!doctype html>
<html lang="en">

<head>

<meta charset="utf-8">

<meta name="viewport"
      content="width=device-width,initial-scale=1">

<meta name="author"
      content="{html.escape(AUTHOR)}">

<title>
{html.escape(report["sample_id"])} — Alignment QC Dashboard
</title>

<style>

:root {{
    --bg:#07101f;
    --bg2:#0b1425;
    --panel:#101b2e;
    --panel2:#15233a;
    --border:#243552;
    --text:#edf3fb;
    --muted:#9aabc3;

    --pass:#41d17d;
    --warn:#f6c85f;
    --fail:#ff6b6b;
    --na:#8090a8;

    --accent:#5fa8ff;
    --accent2:#886dff;

    --primary:#5fa8ff;
    --secondary:#9f7aea;
    --supplementary:#44c7c7;

    --mapped:#41d17d;
    --unmapped:#ff6b6b;

    --unique:#5fa8ff;
    --duplicate:#f6c85f;
}}

* {{
    box-sizing:border-box;
}}

html {{
    scroll-behavior:smooth;
}}

body {{
    margin:0;
    background:
      radial-gradient(
        circle at top right,
        rgba(95,168,255,.09),
        transparent 30%
      ),
      radial-gradient(
        circle at top left,
        rgba(136,109,255,.06),
        transparent 25%
      ),
      var(--bg);

    color:var(--text);

    font-family:
      Inter,
      -apple-system,
      BlinkMacSystemFont,
      "Segoe UI",
      Arial,
      sans-serif;
}}

.container {{
    width:min(1380px,94vw);
    margin:auto;
    padding:36px 0 64px;
}}

.hero {{
    display:flex;
    justify-content:space-between;
    align-items:flex-start;
    gap:28px;
    margin-bottom:24px;
}}

.eyebrow {{
    color:var(--accent);
    text-transform:uppercase;
    letter-spacing:.14em;
    font-size:12px;
    font-weight:800;
    margin-bottom:8px;
}}

h1 {{
    margin:0;
    font-size:34px;
    letter-spacing:-.025em;
}}

.subtitle {{
    margin-top:9px;
    color:var(--muted);
    line-height:1.6;
}}

.hero-status {{
    min-width:180px;
    border-radius:18px;
    padding:18px 20px;
    border:1px solid var(--border);
    background:var(--panel);
}}

.hero-status span {{
    display:block;
    color:var(--muted);
    font-size:12px;
    margin-bottom:7px;
}}

.hero-status strong {{
    font-size:22px;
}}

.hero-status.pass strong {{
    color:var(--pass);
}}

.hero-status.warn strong {{
    color:var(--warn);
}}

.hero-status.fail strong {{
    color:var(--fail);
}}

.pipeline {{
    display:flex;
    align-items:center;
    gap:10px;
    overflow-x:auto;
    padding:15px 18px;
    margin-bottom:24px;

    background:var(--panel);
    border:1px solid var(--border);
    border-radius:16px;
}}

.pipe-stage {{
    white-space:nowrap;
    color:var(--muted);
    font-size:13px;
}}

.pipe-stage.done {{
    color:var(--pass);
    font-weight:700;
}}

.pipe-arrow {{
    color:#52627a;
}}

.gate-grid {{
    display:grid;
    grid-template-columns:1fr 1fr;
    gap:15px;
    margin-bottom:24px;
}}

.gate-card {{
    padding:20px;
    border:1px solid var(--border);
    border-radius:17px;
    background:var(--panel);
}}

.gate-card h3 {{
    margin:0 0 8px;
    color:var(--muted);
    font-size:13px;
    text-transform:uppercase;
    letter-spacing:.07em;
}}

.gate-big {{
    font-size:27px;
    font-weight:800;
}}

.gate-big.pass {{
    color:var(--pass);
}}

.gate-big.warn {{
    color:var(--warn);
}}

.gate-big.fail {{
    color:var(--fail);
}}

.gate-big.na {{
    color:var(--na);
}}

.gate-card p {{
    margin:9px 0 0;
    color:var(--muted);
    line-height:1.55;
    font-size:13px;
}}

.metrics {{
    display:grid;
    grid-template-columns:
      repeat(auto-fit,minmax(190px,1fr));

    gap:14px;
    margin-bottom:24px;
}}

.metric-card {{
    border:1px solid var(--border);
    border-radius:17px;
    background:
      linear-gradient(
        145deg,
        rgba(255,255,255,.018),
        transparent
      ),
      var(--panel);

    padding:18px;
    min-height:152px;
}}

.metric-top {{
    display:flex;
    justify-content:space-between;
    gap:8px;
    align-items:center;
}}

.metric-label {{
    color:var(--muted);
    font-size:13px;
}}

.metric-value {{
    font-size:30px;
    font-weight:800;
    margin:19px 0 7px;
    letter-spacing:-.025em;
}}

.metric-value span {{
    font-size:16px;
    color:var(--muted);
}}

.metric-subtitle {{
    color:var(--muted);
    font-size:12px;
    line-height:1.45;
}}

.mini-status,
.table-status {{
    font-size:10px;
    font-weight:800;
    letter-spacing:.05em;
    padding:4px 8px;
    border-radius:100px;
}}

.pass {{
    color:var(--pass);
}}

.warn {{
    color:var(--warn);
}}

.fail {{
    color:var(--fail);
}}

.na {{
    color:var(--na);
}}

.mini-status.pass,
.table-status.pass {{
    background:rgba(65,209,125,.12);
}}

.mini-status.warn,
.table-status.warn {{
    background:rgba(246,200,95,.12);
}}

.mini-status.fail,
.table-status.fail {{
    background:rgba(255,107,107,.12);
}}

.mini-status.na,
.table-status.na {{
    background:rgba(128,144,168,.12);
}}

.notice {{
    display:flex;
    gap:14px;
    border:1px solid rgba(95,168,255,.32);
    background:rgba(95,168,255,.07);
    border-radius:16px;
    padding:16px 18px;
    margin-bottom:24px;
}}

.notice-icon {{
    flex:0 0 30px;
    width:30px;
    height:30px;
    border-radius:50%;
    display:flex;
    align-items:center;
    justify-content:center;
    background:rgba(95,168,255,.18);
    color:var(--accent);
    font-weight:900;
}}

.notice p {{
    margin:5px 0 0;
    color:var(--muted);
    line-height:1.55;
}}

.section-grid {{
    display:grid;
    grid-template-columns:1.6fr 1fr;
    gap:16px;
    margin-bottom:16px;
}}

.section {{
    border:1px solid var(--border);
    border-radius:18px;
    background:var(--panel);
    padding:21px;
}}

.section h2 {{
    margin:0;
    font-size:18px;
}}

.section-heading {{
    display:flex;
    justify-content:space-between;
    align-items:center;
    gap:10px;
    margin-bottom:20px;
}}

.section-description {{
    color:var(--muted);
    font-size:12px;
}}

.chrom-row {{
    display:grid;
    grid-template-columns:58px 1fr 70px;
    gap:10px;
    align-items:center;
    margin:7px 0;
}}

.chrom-name {{
    color:var(--muted);
    font-size:11px;
}}

.chrom-track,
.threshold-track {{
    height:8px;
    background:#1c2a42;
    border-radius:100px;
    overflow:hidden;
}}

.chrom-bar {{
    height:100%;
    border-radius:100px;
    background:
      linear-gradient(
        90deg,
        var(--accent),
        var(--accent2)
      );

    min-width:0;
}}

.chrom-value {{
    text-align:right;
    font-size:11px;
    color:var(--muted);
    font-variant-numeric:tabular-nums;
}}

.threshold-row {{
    display:grid;
    grid-template-columns:50px 1fr 72px;
    gap:10px;
    align-items:center;
    margin:15px 0;
}}

.threshold-label {{
    font-size:12px;
    color:var(--muted);
}}

.threshold-bar {{
    height:100%;
    background:
      linear-gradient(
        90deg,
        var(--accent),
        var(--pass)
      );

    border-radius:100px;
}}

.threshold-value {{
    text-align:right;
    font-size:12px;
    font-variant-numeric:tabular-nums;
}}

.donut-wrapper {{
    display:flex;
    justify-content:center;
    padding:12px 0 18px;
}}

.donut {{
    width:190px;
    height:190px;
    border-radius:50%;
    display:flex;
    align-items:center;
    justify-content:center;
}}

.donut-inner {{
    width:132px;
    height:132px;
    background:var(--panel);
    border-radius:50%;
    display:flex;
    flex-direction:column;
    align-items:center;
    justify-content:center;
}}

.donut-number {{
    font-size:32px;
    font-weight:800;
}}

.donut-label {{
    color:var(--muted);
    font-size:12px;
}}

.stack-bar {{
    display:flex;
    width:100%;
    height:22px;
    overflow:hidden;
    border-radius:100px;
    background:#1c2a42;
    margin:16px 0;
}}

.stack-segment {{
    height:100%;
}}

.stack-segment.primary {{
    background:var(--primary);
}}

.stack-segment.secondary {{
    background:var(--secondary);
}}

.stack-segment.supplementary {{
    background:var(--supplementary);
}}

.stack-segment.mapped {{
    background:var(--mapped);
}}

.stack-segment.unmapped {{
    background:var(--unmapped);
}}

.stack-segment.unique {{
    background:var(--unique);
}}

.stack-segment.duplicate {{
    background:var(--duplicate);
}}

.legend {{
    display:grid;
    gap:9px;
}}

.legend-item {{
    display:grid;
    grid-template-columns:12px 1fr auto auto;
    gap:8px;
    align-items:center;
    font-size:12px;
    color:var(--muted);
}}

.legend-item strong {{
    color:var(--text);
}}

.legend-item small {{
    min-width:58px;
    text-align:right;
}}

.legend-dot {{
    width:9px;
    height:9px;
    border-radius:50%;
}}

.legend-dot.primary {{
    background:var(--primary);
}}

.legend-dot.secondary {{
    background:var(--secondary);
}}

.legend-dot.supplementary {{
    background:var(--supplementary);
}}

.legend-dot.mapped {{
    background:var(--mapped);
}}

.legend-dot.unmapped {{
    background:var(--unmapped);
}}

.legend-dot.unique {{
    background:var(--unique);
}}

.legend-dot.duplicate {{
    background:var(--duplicate);
}}

.triple-grid {{
    display:grid;
    grid-template-columns:repeat(3,1fr);
    gap:16px;
    margin-bottom:16px;
}}

details {{
    border:1px solid var(--border);
    border-radius:17px;
    background:var(--panel);
    margin-top:15px;
    overflow:hidden;
}}

summary {{
    cursor:pointer;
    padding:17px 20px;
    font-weight:700;
    user-select:none;
}}

.details-content {{
    border-top:1px solid var(--border);
    padding:18px 20px;
}}

table {{
    width:100%;
    border-collapse:collapse;
    font-size:12px;
}}

th,
td {{
    padding:10px 11px;
    border-bottom:1px solid var(--border);
    text-align:left;
}}

th {{
    color:var(--muted);
    font-size:11px;
    text-transform:uppercase;
    letter-spacing:.04em;
}}

td {{
    color:#dce5f2;
}}

.provenance-grid {{
    display:grid;
    grid-template-columns:
      repeat(auto-fit,minmax(210px,1fr));
    gap:10px;
}}

.provenance-item {{
    padding:12px 14px;
    border-radius:12px;
    background:var(--panel2);
}}

.provenance-item span {{
    display:block;
    color:var(--muted);
    font-size:10px;
    margin-bottom:5px;
}}

.provenance-item strong {{
    font-size:12px;
    word-break:break-all;
}}

.footer {{
    color:var(--muted);
    font-size:11px;
    line-height:1.6;
    margin-top:26px;
}}

@media(max-width:900px) {{

    .section-grid,
    .gate-grid,
    .triple-grid {{
        grid-template-columns:1fr;
    }}

    .hero {{
        flex-direction:column;
    }}

    .hero-status {{
        width:100%;
    }}
}}

</style>

</head>


<body>

<div class="container">

<section class="hero">

<div>

<div class="eyebrow">
Genomics platform · Alignment QC
</div>

<h1>
{html.escape(report["sample_id"])}
</h1>

<div class="subtitle">
{html.escape(profile.get("name","Unknown profile"))}
<br>
QC Engine v{QC_ENGINE_VERSION}
</div>

</div>


<div class="hero-status {status_class(overall_gate)}">

<span>
Overall QC status
</span>

<strong>
{html.escape(overall_gate)}
</strong>

</div>

</section>


<div class="pipeline">

<span class="pipe-stage done">
✓ Raw QC
</span>

<span class="pipe-arrow">→</span>

<span class="pipe-stage done">
✓ Alignment
</span>

<span class="pipe-arrow">→</span>

<span class="pipe-stage done">
✓ Alignment QC
</span>

<span class="pipe-arrow">→</span>

<span class="pipe-stage">
○ Variant calling
</span>

<span class="pipe-arrow">→</span>

<span class="pipe-stage">
○ Variant QC
</span>

<span class="pipe-arrow">→</span>

<span class="pipe-stage">
○ GIAB benchmark
</span>

</div>


<div class="gate-grid">

<div class="gate-card">

<h3>
Technical alignment quality
</h3>

<div class="gate-big {status_class(technical_gate)}">
{html.escape(technical_gate)}
</div>

<p>
{html.escape(technical_message)}
</p>

</div>


<div class="gate-card">

<h3>
Coverage sufficiency
</h3>

<div class="gate-big {status_class(coverage_gate)}">
{html.escape(coverage_gate)}
</div>

<p>
{html.escape(coverage_message)}
</p>

</div>

</div>


{smoke_notice}


<div class="metrics">

{cards}

</div>


<div class="section-grid">

<section class="section">

<div class="section-heading">

<div>
<h2>Chromosome coverage</h2>
<div class="section-description">
Mean coverage across primary GRCh38 chromosomes
</div>
</div>

</div>

{chromosome_chart(coverage_rows)}

</section>


<section class="section">

<div class="section-heading">

<div>
<h2>QC decision profile</h2>
<div class="section-description">
Technical checks evaluated by this profile
</div>
</div>

</div>

<div class="donut-wrapper">
{donut_svg(pass_count,warn_count,fail_count)}
</div>

<div class="legend">

<div class="legend-item">
<span class="legend-dot mapped"></span>
<span>PASS</span>
<strong>{pass_count}</strong>
<small>checks</small>
</div>

<div class="legend-item">
<span class="legend-dot duplicate"></span>
<span>WARN</span>
<strong>{warn_count}</strong>
<small>checks</small>
</div>

<div class="legend-item">
<span class="legend-dot unmapped"></span>
<span>FAIL</span>
<strong>{fail_count}</strong>
<small>checks</small>
</div>

</div>

</section>

</div>


<div class="triple-grid">

<section class="section">

<h2>Read composition</h2>

{stacked_bar([
    ("Primary",primary,"primary"),
    ("Secondary",secondary,"secondary"),
    ("Supplementary",supplementary,"supplementary"),
])}

</section>


<section class="section">

<h2>Primary mapping</h2>

{stacked_bar([
    ("Mapped",mapped,"mapped"),
    ("Unmapped",unmapped,"unmapped"),
])}

</section>


<section class="section">

<h2>Duplicate composition</h2>

{stacked_bar([
    ("Unique",unique_reads,"unique"),
    ("Duplicates",duplicate_reads,"duplicate"),
])}

</section>

</div>


<div class="section-grid">

<section class="section">

<div class="section-heading">

<div>
<h2>Coverage thresholds</h2>

<div class="section-description">
Fraction of reference at or above selected depths
</div>

</div>

</div>

{threshold_chart(metrics)}

</section>


<section class="section">

<h2>Interpretation</h2>

<p class="subtitle">
Technical alignment quality and genome-wide coverage
are evaluated independently.
</p>

<p class="subtitle">
This distinction prevents a deliberately downsampled
smoke-test dataset from being interpreted as a failed
alignment solely because genome-wide coverage is low.
</p>

<p class="subtitle">
These thresholds are intended for research,
development and pipeline validation only.
They are not clinically validated diagnostic
acceptance criteria.
</p>

</section>

</div>


<details open>

<summary>
QC decisions
</summary>

<div class="details-content">

<table>

<thead>
<tr>
<th>Metric</th>
<th>Observed</th>
<th>Status</th>
<th>Interpretation</th>
</tr>
</thead>

<tbody>
{event_rows}
</tbody>

</table>

</div>

</details>


<details>

<summary>
Contig-level coverage
</summary>

<div class="details-content">

<table>

<thead>
<tr>
<th>Contig</th>
<th>Length</th>
<th>Mean coverage</th>
</tr>
</thead>

<tbody>
{contig_rows}
</tbody>

</table>

</div>

</details>


<details>

<summary>
Provenance
</summary>

<div class="details-content">

<div class="provenance-grid">

<div class="provenance-item">
<span>Sample</span>
<strong>{html.escape(report["sample_id"])}</strong>
</div>

<div class="provenance-item">
<span>QC engine</span>
<strong>{QC_ENGINE_VERSION}</strong>
</div>

<div class="provenance-item">
<span>Report schema</span>
<strong>{REPORT_SCHEMA_VERSION}</strong>
</div>

<div class="provenance-item">
<span>QC profile</span>
<strong>{html.escape(report["profile_id"])}</strong>
</div>

<div class="provenance-item">
<span>Generated</span>
<strong>{generated}</strong>
</div>

<div class="provenance-item">
<span>Author</span>
<strong>{html.escape(AUTHOR)}</strong>
</div>

</div>

</div>

</details>


<div class="footer">

Generated by the Modular Human Genomics Pipeline.
Research, benchmarking and pipeline-development use only.
A PASS result does not establish clinical validity,
diagnostic accuracy or clinical interpretation.

</div>

</div>

</body>

</html>
"""


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="Alignment QC engine and HTML dashboard."
    )

    parser.add_argument(
        "--sample-id",
        required=True,
    )

    parser.add_argument(
        "--flagstat",
        required=True,
    )

    parser.add_argument(
        "--stats",
        required=True,
    )

    parser.add_argument(
        "--idxstats",
        required=True,
    )

    parser.add_argument(
        "--mosdepth-summary",
        required=True,
    )

    parser.add_argument(
        "--mosdepth-global-dist",
        required=True,
    )

    parser.add_argument(
        "--thresholds",
        required=True,
    )

    parser.add_argument(
        "--profile",
        required=True,
    )

    parser.add_argument(
        "--output-json",
        required=True,
    )

    parser.add_argument(
        "--output-html",
        required=True,
    )

    parser.add_argument(
        "--log",
        required=True,
    )

    args = parser.parse_args()

    with open(args.thresholds) as handle:
        config = yaml.safe_load(handle)

    profiles = config.get(
        "profiles",
        {},
    )

    if args.profile not in profiles:
        raise ValueError(
            f"Unknown alignment QC profile: {args.profile}. "
            f"Available profiles: {', '.join(profiles)}"
        )

    profile = profiles[
        args.profile
    ]

    thresholds = profile.get(
        "thresholds",
        {},
    )

    coverage_applicable = bool(
        profile.get(
            "coverage_applicable",
            True,
        )
    )

    flagstat = parse_flagstat(
        args.flagstat
    )

    samtools_stats = parse_samtools_stats(
        args.stats
    )

    idxstats = parse_idxstats(
        args.idxstats
    )

    (
        coverage_rows,
        mean_autosomal,
    ) = parse_mosdepth_summary(
        args.mosdepth_summary
    )

    depth_fractions = parse_global_dist(
        args.mosdepth_global_dist
    )

    metrics = dict(
        flagstat
    )

    metrics[
        "mean_autosomal_coverage"
    ] = mean_autosomal

    metrics[
        "genome_1x_fraction_pct"
    ] = depth_fractions.get(1)

    metrics[
        "genome_10x_fraction_pct"
    ] = depth_fractions.get(10)

    metrics[
        "genome_20x_fraction_pct"
    ] = depth_fractions.get(20)

    metrics[
        "genome_30x_fraction_pct"
    ] = depth_fractions.get(30)

    events = []

    technical_checks = [
        (
            "mapping_rate_pct",
            "Primary mapping rate",
            evaluate_lower,
        ),
        (
            "properly_paired_pct",
            "Properly paired reads",
            evaluate_lower,
        ),
        (
            "duplicate_rate_pct",
            "Duplicate rate",
            evaluate_upper,
        ),
    ]

    technical_statuses = []

    for (
        metric,
        label,
        evaluator,
    ) in technical_checks:

        value = metrics.get(
            metric
        )

        status, reason = evaluator(
            metric,
            value,
            thresholds,
        )

        technical_statuses.append(
            status
        )

        events.append(
            {
                "metric": metric,
                "label": label,
                "value": value,
                "status": status,
                "reason": reason,
                "category": "technical",
            }
        )

    technical_gate = combine_gate(
        technical_statuses
    )

    coverage_statuses = []

    coverage_checks = [
        (
            "mean_autosomal_coverage",
            "Mean autosomal coverage",
            evaluate_lower,
        ),
        (
            "genome_20x_fraction_pct",
            "Genome ≥20×",
            evaluate_lower,
        ),
    ]

    for (
        metric,
        label,
        evaluator,
    ) in coverage_checks:

        value = metrics.get(
            metric
        )

        if coverage_applicable:

            status, reason = evaluator(
                metric,
                value,
                thresholds,
            )

            coverage_statuses.append(
                status
            )

        else:

            status = "N/A"

            reason = (
                "Coverage sufficiency is intentionally "
                "not assessed in smoke-test mode"
            )

        events.append(
            {
                "metric": metric,
                "label": label,
                "value": value,
                "status": status,
                "reason": reason,
                "category": "coverage",
            }
        )

    coverage_gate = (
        combine_gate(
            coverage_statuses
        )
        if coverage_applicable
        else "NOT_APPLICABLE"
    )

    overall_gate = (
        combine_gate(
            [
                technical_gate,
                coverage_gate,
            ]
        )
        if coverage_applicable
        else technical_gate
    )

    report = {
        "schema_version": REPORT_SCHEMA_VERSION,
        "qc_engine_version": QC_ENGINE_VERSION,
        "author": AUTHOR,

        "generated_at": datetime.now(
            timezone.utc
        ).isoformat(),

        "sample_id": args.sample_id,

        "profile_id": args.profile,
        "profile": profile,

        "coverage_applicable": coverage_applicable,

        "gate": overall_gate,

        "technical_alignment_gate": technical_gate,

        "coverage_gate": coverage_gate,

        "metrics": metrics,

        "events": events,

        "samtools_stats": samtools_stats,

        "idxstats": idxstats,

        "coverage_by_contig": coverage_rows,
    }

    Path(
        args.output_json
    ).write_text(
        json.dumps(
            report,
            indent=2,
        )
    )

    Path(
        args.output_html
    ).write_text(
        render_html(
            report
        )
    )

    Path(
        args.log
    ).write_text(
        "\n".join(
            [
                f"sample_id={args.sample_id}",
                f"profile={args.profile}",
                f"overall_gate={overall_gate}",
                f"technical_alignment_gate={technical_gate}",
                f"coverage_gate={coverage_gate}",
                f"qc_engine_version={QC_ENGINE_VERSION}",
                "",
            ]
        )
    )

    print(
        f"[ALIGNMENT_QC] "
        f"{args.sample_id}: "
        f"overall={overall_gate}, "
        f"technical={technical_gate}, "
        f"coverage={coverage_gate}"
    )


if __name__ == "__main__":
    main()
