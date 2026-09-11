#!/usr/bin/env python3

import argparse
import html
import json
import logging
import math
import sys
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import yaml


LOGGER = logging.getLogger("raw_qc")

QC_ENGINE_VERSION = "0.3.0"
REPORT_SCHEMA_VERSION = "0.3.0"
AUTHOR = "Efe Mıhcı"


# ---------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------

def setup_logging(log_path: Path) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)

    LOGGER.setLevel(logging.INFO)

    formatter = logging.Formatter(
        "%(asctime)s | %(levelname)s | %(message)s"
    )

    file_handler = logging.FileHandler(log_path)
    file_handler.setFormatter(formatter)

    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(formatter)

    LOGGER.handlers.clear()
    LOGGER.addHandler(file_handler)
    LOGGER.addHandler(stream_handler)


# ---------------------------------------------------------------------
# FastQC archive reading
# ---------------------------------------------------------------------

def read_fastqc_zip(zip_path: Path):
    LOGGER.info("Reading FastQC archive: %s", zip_path)

    if not zip_path.exists():
        raise FileNotFoundError(
            f"FastQC ZIP not found: {zip_path}"
        )

    if not zipfile.is_zipfile(zip_path):
        raise ValueError(
            f"Invalid or corrupted FastQC ZIP: {zip_path}"
        )

    with zipfile.ZipFile(zip_path, "r") as zf:
        names = zf.namelist()

        data_files = [
            name
            for name in names
            if name.endswith("/fastqc_data.txt")
        ]

        summary_files = [
            name
            for name in names
            if name.endswith("/summary.txt")
        ]

        if len(data_files) != 1:
            raise RuntimeError(
                f"Expected exactly one fastqc_data.txt in "
                f"{zip_path}; found {len(data_files)}"
            )

        if len(summary_files) != 1:
            raise RuntimeError(
                f"Expected exactly one summary.txt in "
                f"{zip_path}; found {len(summary_files)}"
            )

        data_text = zf.read(
            data_files[0]
        ).decode(
            "utf-8",
            errors="replace"
        )

        summary_text = zf.read(
            summary_files[0]
        ).decode(
            "utf-8",
            errors="replace"
        )

    return data_text, summary_text


# ---------------------------------------------------------------------
# FastQC parsing
# ---------------------------------------------------------------------

def parse_basic_statistics(data_text: str):
    stats = {}
    inside_module = False

    for line in data_text.splitlines():

        if line.startswith(">>Basic Statistics"):
            inside_module = True
            continue

        if inside_module and line.startswith(">>END_MODULE"):
            break

        if (
            inside_module
            and line
            and not line.startswith("#")
        ):
            parts = line.split("\t")

            if len(parts) >= 2:
                stats[
                    parts[0].strip()
                ] = parts[1].strip()

    return stats


def parse_fastqc_modules(summary_text: str):
    modules = {}

    for line in summary_text.splitlines():
        parts = line.split("\t")

        if len(parts) >= 2:
            status = parts[0].strip()
            module = parts[1].strip()

            modules[module] = status

    return modules


def parse_module_tables(data_text: str):
    """
    Extract tabular FastQC module data.

    Returns:
        {
            "Per base sequence quality": {
                "status": "pass",
                "header": [...],
                "rows": [[...], ...]
            },
            ...
        }
    """

    modules = {}

    current_name = None
    current_status = None
    current_header = None
    current_rows = []

    for raw_line in data_text.splitlines():
        line = raw_line.rstrip("\n")

        if line.startswith(">>") and not line.startswith(
            ">>END_MODULE"
        ):
            parts = line[2:].split("\t")

            current_name = parts[0].strip()

            current_status = (
                parts[1].strip()
                if len(parts) > 1
                else "unknown"
            )

            current_header = None
            current_rows = []

            continue

        if line.startswith(">>END_MODULE"):

            if current_name:
                modules[current_name] = {
                    "status": current_status,
                    "header": current_header or [],
                    "rows": current_rows,
                }

            current_name = None
            current_status = None
            current_header = None
            current_rows = []

            continue

        if current_name is None:
            continue

        if line.startswith("#"):
            header = line[1:].split("\t")

            if len(header) > 1:
                current_header = [
                    value.strip()
                    for value in header
                ]

            continue

        if not line.strip():
            continue

        if current_header is not None:
            current_rows.append(
                [
                    value.strip()
                    for value in line.split("\t")
                ]
            )

    return modules


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------

def safe_int(value):
    try:
        return int(
            str(value).replace(",", "")
        )
    except Exception:
        return None


def safe_float(value):
    try:
        return float(value)
    except Exception:
        return None


def detect_read(filename: str):
    name = filename.upper()

    r1_patterns = [
        "_R1",
        ".R1",
        "_1.FASTQ",
        "_1.FQ",
    ]

    r2_patterns = [
        "_R2",
        ".R2",
        "_2.FASTQ",
        "_2.FQ",
    ]

    if any(
        pattern in name
        for pattern in r1_patterns
    ):
        return "R1"

    if any(
        pattern in name
        for pattern in r2_patterns
    ):
        return "R2"

    return "UNKNOWN"


def parse_base_position(value):
    """
    FastQC can report:
      1
      10-14
      145-149

    Return midpoint for plotting.
    """

    value = str(value).strip()

    if "-" in value:
        try:
            start, end = value.split("-", 1)

            return (
                float(start)
                + float(end)
            ) / 2
        except Exception:
            return None

    return safe_float(value)


def parse_archive(zip_path: Path):
    data_text, summary_text = read_fastqc_zip(
        zip_path
    )

    basic = parse_basic_statistics(data_text)

    modules = parse_fastqc_modules(
        summary_text
    )

    module_tables = parse_module_tables(
        data_text
    )

    filename = basic.get(
        "Filename",
        zip_path.name
    )

    return {
        "archive": str(zip_path),
        "filename": filename,
        "read": detect_read(filename),
        "total_sequences": safe_int(
            basic.get(
                "Total Sequences",
                ""
            )
        ),
        "sequence_length": basic.get(
            "Sequence length"
        ),
        "gc_percent": safe_float(
            basic.get(
                "%GC",
                ""
            )
        ),
        "encoding": basic.get(
            "Encoding"
        ),
        "file_type": basic.get(
            "File type"
        ),
        "modules": modules,
        "module_tables": module_tables,
    }


# ---------------------------------------------------------------------
# QC gate
# ---------------------------------------------------------------------

def severity_rank(status: str) -> int:
    return {
        "PASS": 0,
        "WARN": 1,
        "FAIL": 2,
    }.get(
        status,
        1
    )


def evaluate_qc(reads, thresholds):
    events = []
    gate = "PASS"

    def add_event(
        level,
        code,
        message
    ):
        nonlocal gate

        events.append({
            "level": level,
            "code": code,
            "message": message,
        })

        if (
            severity_rank(level)
            > severity_rank(gate)
        ):
            gate = level

    fastqc_config = thresholds.get(
        "fastqc_modules",
        {}
    )

    critical_modules = fastqc_config.get(
        "critical",
        []
    )

    warning_modules = fastqc_config.get(
        "warning_only",
        []
    )

    for read in reads:
        read_name = read["read"]

        for module in critical_modules:
            status = read[
                "modules"
            ].get(module)

            if status == "FAIL":
                add_event(
                    "FAIL",
                    "FASTQC_CRITICAL_MODULE_FAIL",
                    f"{read_name}: "
                    f"{module} failed",
                )

            elif status == "WARN":
                add_event(
                    "WARN",
                    "FASTQC_CRITICAL_MODULE_WARN",
                    f"{read_name}: "
                    f"{module} warning",
                )

        for module in warning_modules:
            status = read[
                "modules"
            ].get(module)

            if status in {
                "WARN",
                "FAIL"
            }:
                add_event(
                    "WARN",
                    "FASTQC_NONCRITICAL_MODULE",
                    f"{read_name}: "
                    f"{module} = {status}",
                )

    r1 = next(
        (
            x
            for x in reads
            if x["read"] == "R1"
        ),
        None,
    )

    r2 = next(
        (
            x
            for x in reads
            if x["read"] == "R2"
        ),
        None,
    )

    if not r1 or not r2:
        add_event(
            "FAIL",
            "PAIRED_READ_MISSING",
            "Could not identify both R1 and R2",
        )

    elif (
        r1["total_sequences"] is not None
        and r2["total_sequences"] is not None
    ):

        maximum = max(
            r1["total_sequences"],
            r2["total_sequences"],
        )

        difference = abs(
            r1["total_sequences"]
            - r2["total_sequences"]
        )

        difference_pct = (
            difference
            / maximum
            * 100
            if maximum
            else 0
        )

        pair_config = thresholds.get(
            "paired_end",
            {}
        ).get(
            "max_read_count_difference_pct",
            {}
        )

        warn_threshold = pair_config.get(
            "warn",
            0.5
        )

        fail_threshold = pair_config.get(
            "fail",
            2.0
        )

        if difference_pct >= fail_threshold:
            add_event(
                "FAIL",
                "READ_COUNT_MISMATCH",
                f"R1/R2 read count "
                f"difference "
                f"{difference_pct:.3f}%",
            )

        elif difference_pct >= warn_threshold:
            add_event(
                "WARN",
                "READ_COUNT_MISMATCH",
                f"R1/R2 read count "
                f"difference "
                f"{difference_pct:.3f}%",
            )

    minimum_sequences = thresholds.get(
        "basic_statistics",
        {}
    ).get(
        "minimum_total_sequences",
        {}
    )

    warn_min = minimum_sequences.get(
        "warn"
    )

    fail_min = minimum_sequences.get(
        "fail"
    )

    for read in reads:
        total = read.get(
            "total_sequences"
        )

        if total is None:
            add_event(
                "WARN",
                "TOTAL_SEQUENCES_UNKNOWN",
                f"{read['read']}: "
                f"total sequence count "
                f"could not be parsed",
            )
            continue

        if (
            fail_min is not None
            and total < fail_min
        ):
            add_event(
                "FAIL",
                "LOW_READ_COUNT",
                f"{read['read']}: only "
                f"{total:,} reads detected",
            )

        elif (
            warn_min is not None
            and total < warn_min
        ):
            add_event(
                "WARN",
                "LOW_READ_COUNT",
                f"{read['read']}: only "
                f"{total:,} reads detected",
            )

    if not events:
        events.append({
            "level": "PASS",
            "code": "RAW_QC_OK",
            "message": (
                "No raw sequencing QC "
                "problems detected"
            ),
        })

    return gate, events


# ---------------------------------------------------------------------
# Summary helpers
# ---------------------------------------------------------------------

def module_status_counts(reads):
    counts = {
        "PASS": 0,
        "WARN": 0,
        "FAIL": 0,
    }

    for read in reads:
        for status in read[
            "modules"
        ].values():

            if status in counts:
                counts[status] += 1

    return counts


def event_status_counts(events):
    counts = {
        "PASS": 0,
        "WARN": 0,
        "FAIL": 0,
    }

    for event in events:
        level = event["level"]

        if level in counts:
            counts[level] += 1

    return counts


def pretty_count(value):
    if value is None:
        return "NA"

    if value >= 1_000_000_000:
        return (
            f"{value / 1_000_000_000:.2f}B"
        )

    if value >= 1_000_000:
        return (
            f"{value / 1_000_000:.1f}M"
        )

    if value >= 1_000:
        return (
            f"{value / 1_000:.1f}K"
        )

    return str(value)


# ---------------------------------------------------------------------
# Extract representative chart data
# ---------------------------------------------------------------------

def extract_per_base_quality(read):
    table = read[
        "module_tables"
    ].get(
        "Per base sequence quality"
    )

    if not table:
        return []

    header = table.get(
        "header",
        []
    )

    try:
        base_index = header.index("Base")
        mean_index = header.index("Mean")
    except ValueError:
        return []

    points = []

    for row in table.get(
        "rows",
        []
    ):
        if (
            len(row)
            <= max(
                base_index,
                mean_index
            )
        ):
            continue

        x = parse_base_position(
            row[base_index]
        )

        y = safe_float(
            row[mean_index]
        )

        if x is None or y is None:
            continue

        points.append(
            (x, y)
        )

    return points


def extract_gc_distribution(read):
    table = read[
        "module_tables"
    ].get(
        "Per sequence GC content"
    )

    if not table:
        return []

    header = table.get(
        "header",
        []
    )

    try:
        gc_index = header.index(
            "GC Content"
        )
        count_index = header.index(
            "Count"
        )
    except ValueError:
        return []

    points = []

    for row in table.get(
        "rows",
        []
    ):

        if (
            len(row)
            <= max(
                gc_index,
                count_index
            )
        ):
            continue

        x = safe_float(
            row[gc_index]
        )

        y = safe_float(
            row[count_index]
        )

        if x is None or y is None:
            continue

        points.append(
            (x, y)
        )

    return points


def extract_adapter_max(read):
    table = read[
        "module_tables"
    ].get(
        "Adapter Content"
    )

    if not table:
        return []

    header = table.get(
        "header",
        []
    )

    if len(header) < 2:
        return []

    points = []

    for row in table.get(
        "rows",
        []
    ):

        if len(row) < 2:
            continue

        position = parse_base_position(
            row[0]
        )

        values = [
            safe_float(value)
            for value in row[1:]
        ]

        values = [
            value
            for value in values
            if value is not None
        ]

        if (
            position is None
            or not values
        ):
            continue

        points.append(
            (
                position,
                max(values)
            )
        )

    return points


# ---------------------------------------------------------------------
# SVG chart engine
# ---------------------------------------------------------------------

def svg_line_chart(
    series,
    title,
    x_label,
    y_label,
    y_min=None,
    y_max=None,
):
    """
    series:
    [
        {
            "name": "R1",
            "points": [(x, y), ...],
            "class": "line-r1"
        },
        ...
    ]
    """

    valid_series = [
        s
        for s in series
        if s.get("points")
    ]

    if not valid_series:
        return """
        <div class="empty-chart">
            Data unavailable for this chart.
        </div>
        """

    all_points = [
        point
        for item in valid_series
        for point in item["points"]
    ]

    xs = [
        point[0]
        for point in all_points
    ]

    ys = [
        point[1]
        for point in all_points
    ]

    x_min = min(xs)
    x_max = max(xs)

    calculated_y_min = min(ys)
    calculated_y_max = max(ys)

    if y_min is None:
        y_min = calculated_y_min

    if y_max is None:
        y_max = calculated_y_max

    if x_max == x_min:
        x_max = x_min + 1

    if y_max == y_min:
        y_max = y_min + 1

    width = 900
    height = 360

    left = 76
    right = 28
    top = 28
    bottom = 62

    plot_width = (
        width
        - left
        - right
    )

    plot_height = (
        height
        - top
        - bottom
    )

    def sx(x):
        return (
            left
            + (
                (x - x_min)
                / (x_max - x_min)
            )
            * plot_width
        )

    def sy(y):
        return (
            top
            + (
                1
                - (
                    (y - y_min)
                    / (y_max - y_min)
                )
            )
            * plot_height
        )

    grid = []

    for i in range(5):
        frac = i / 4

        y_value = (
            y_max
            - frac
            * (
                y_max
                - y_min
            )
        )

        y_pos = (
            top
            + frac
            * plot_height
        )

        grid.append(
            f"""
            <line
                x1="{left}"
                y1="{y_pos:.2f}"
                x2="{left + plot_width}"
                y2="{y_pos:.2f}"
                class="grid-line"
            />
            <text
                x="{left - 12}"
                y="{y_pos + 4:.2f}"
                text-anchor="end"
                class="axis-text"
            >
                {y_value:.1f}
            </text>
            """
        )

    x_ticks = []

    for i in range(6):
        frac = i / 5

        x_value = (
            x_min
            + frac
            * (
                x_max
                - x_min
            )
        )

        x_pos = (
            left
            + frac
            * plot_width
        )

        x_ticks.append(
            f"""
            <text
                x="{x_pos:.2f}"
                y="{height - 31}"
                text-anchor="middle"
                class="axis-text"
            >
                {x_value:.0f}
            </text>
            """
        )

    lines = []

    legend = []

    for index, item in enumerate(
        valid_series
    ):
        points = item[
            "points"
        ]

        css_class = item.get(
            "class",
            f"line-s{index}"
        )

        coords = " ".join(
            f"{sx(x):.2f},{sy(y):.2f}"
            for x, y in points
        )

        circles = []

        for x, y in points:
            circles.append(
                f"""
                <circle
                    cx="{sx(x):.2f}"
                    cy="{sy(y):.2f}"
                    r="7"
                    class="hover-point {css_class}-point"
                >
                    <title>
                        {html.escape(item['name'])}
                        — {x_label}: {x:g},
                        {y_label}: {y:.3f}
                    </title>
                </circle>
                """
            )

        lines.append(
            f"""
            <polyline
                points="{coords}"
                class="chart-line {css_class}"
            />
            {''.join(circles)}
            """
        )

        legend.append(
            f"""
            <span class="legend-item">
                <span
                    class="legend-line {css_class}"
                ></span>
                {html.escape(item['name'])}
            </span>
            """
        )

    return f"""
    <div class="chart-header">
        <div>
            <h3>{html.escape(title)}</h3>
        </div>
        <div class="legend">
            {''.join(legend)}
        </div>
    </div>

    <div class="chart-scroll">
        <svg
            viewBox="0 0 {width} {height}"
            class="chart-svg"
            role="img"
            aria-label="{html.escape(title)}"
        >

            {''.join(grid)}

            <line
                x1="{left}"
                y1="{top}"
                x2="{left}"
                y2="{top + plot_height}"
                class="axis-line"
            />

            <line
                x1="{left}"
                y1="{top + plot_height}"
                x2="{left + plot_width}"
                y2="{top + plot_height}"
                class="axis-line"
            />

            {''.join(x_ticks)}

            {''.join(lines)}

            <text
                x="{left + plot_width / 2}"
                y="{height - 5}"
                text-anchor="middle"
                class="axis-label"
            >
                {html.escape(x_label)}
            </text>

            <text
                transform="
                    translate(
                        18
                        {top + plot_height / 2}
                    )
                    rotate(-90)
                "
                text-anchor="middle"
                class="axis-label"
            >
                {html.escape(y_label)}
            </text>

        </svg>
    </div>
    """


def status_donut(counts):
    total = sum(
        counts.values()
    )

    if total == 0:
        return """
        <div class="empty-chart">
            No module data.
        </div>
        """

    pass_pct = (
        counts["PASS"]
        / total
        * 100
    )

    warn_pct = (
        counts["WARN"]
        / total
        * 100
    )

    fail_pct = (
        counts["FAIL"]
        / total
        * 100
    )

    return f"""
    <div class="donut-layout">

        <div
            class="donut"
            style="
                --pass:{pass_pct:.2f};
                --warn:{warn_pct:.2f};
                --fail:{fail_pct:.2f};
            "
        >
            <div class="donut-center">
                <strong>{total}</strong>
                <span>checks</span>
            </div>
        </div>

        <div class="donut-legend">

            <div>
                <span class="dot dot-pass"></span>
                PASS
                <strong>
                    {counts["PASS"]}
                </strong>
            </div>

            <div>
                <span class="dot dot-warn"></span>
                WARN
                <strong>
                    {counts["WARN"]}
                </strong>
            </div>

            <div>
                <span class="dot dot-fail"></span>
                FAIL
                <strong>
                    {counts["FAIL"]}
                </strong>
            </div>

        </div>

    </div>
    """


# ---------------------------------------------------------------------
# HTML report
# ---------------------------------------------------------------------

def build_html(
    sample_id,
    gate,
    reads,
    events,
    output_path,
):
    module_counts = module_status_counts(
        reads
    )

    r1 = next(
        (
            read
            for read in reads
            if read["read"] == "R1"
        ),
        None,
    )

    r2 = next(
        (
            read
            for read in reads
            if read["read"] == "R2"
        ),
        None,
    )

    representative = (
        r1
        or r2
        or (
            reads[0]
            if reads
            else None
        )
    )

    read_count = (
        pretty_count(
            representative[
                "total_sequences"
            ]
        )
        if representative
        else "NA"
    )

    read_length = (
        str(
            representative[
                "sequence_length"
            ]
        )
        if representative
        else "NA"
    )

    gc_percent = (
        f"{representative['gc_percent']:.1f}%"
        if (
            representative
            and representative[
                "gc_percent"
            ] is not None
        )
        else "NA"
    )

    total_pairs = None

    if (
        r1
        and r2
        and r1["total_sequences"] is not None
        and r2["total_sequences"] is not None
    ):
        total_pairs = min(
            r1["total_sequences"],
            r2["total_sequences"],
        )

    pair_count = pretty_count(
        total_pairs
    )

    # --------------------------------------------------------------
    # Representative charts
    # --------------------------------------------------------------

    quality_series = []

    if r1:
        quality_series.append({
            "name": "R1",
            "points": extract_per_base_quality(
                r1
            ),
            "class": "line-r1",
        })

    if r2:
        quality_series.append({
            "name": "R2",
            "points": extract_per_base_quality(
                r2
            ),
            "class": "line-r2",
        })

    quality_chart = svg_line_chart(
        quality_series,
        "Per-base mean sequence quality",
        "Read position (bp)",
        "Mean Phred quality",
        y_min=0,
        y_max=45,
    )

    gc_series = []

    if r1:
        gc_series.append({
            "name": "R1",
            "points": extract_gc_distribution(
                r1
            ),
            "class": "line-r1",
        })

    if r2:
        gc_series.append({
            "name": "R2",
            "points": extract_gc_distribution(
                r2
            ),
            "class": "line-r2",
        })

    gc_chart = svg_line_chart(
        gc_series,
        "Per-sequence GC distribution",
        "GC content (%)",
        "Sequence count",
    )

    adapter_series = []

    if r1:
        adapter_series.append({
            "name": "R1",
            "points": extract_adapter_max(
                r1
            ),
            "class": "line-r1",
        })

    if r2:
        adapter_series.append({
            "name": "R2",
            "points": extract_adapter_max(
                r2
            ),
            "class": "line-r2",
        })

    adapter_chart = svg_line_chart(
        adapter_series,
        "Maximum adapter content",
        "Read position (bp)",
        "Adapter content (%)",
        y_min=0,
    )

    # --------------------------------------------------------------
    # Read table
    # --------------------------------------------------------------

    read_rows = []

    for read in reads:

        total_sequences = (
            f"{read['total_sequences']:,}"
            if read[
                "total_sequences"
            ] is not None
            else "NA"
        )

        gc_value = (
            f"{read['gc_percent']:.1f}%"
            if read[
                "gc_percent"
            ] is not None
            else "NA"
        )

        read_rows.append(
            f"""
            <tr>
                <td>
                    <span class="read-pill">
                        {
                            html.escape(
                                str(
                                    read["read"]
                                )
                            )
                        }
                    </span>
                </td>

                <td class="filename">
                    {
                        html.escape(
                            str(
                                read["filename"]
                            )
                        )
                    }
                </td>

                <td>
                    {total_sequences}
                </td>

                <td>
                    {
                        html.escape(
                            str(
                                read[
                                    "sequence_length"
                                ]
                            )
                        )
                    }
                </td>

                <td>
                    {gc_value}
                </td>

                <td>
                    {
                        html.escape(
                            str(
                                read[
                                    "encoding"
                                ]
                            )
                        )
                    }
                </td>
            </tr>
            """
        )

    # --------------------------------------------------------------
    # QC events
    # --------------------------------------------------------------

    event_cards = []

    for event in events:

        icon = {
            "PASS": "✓",
            "WARN": "!",
            "FAIL": "×",
        }.get(
            event["level"],
            "•"
        )

        event_cards.append(
            f"""
            <div class="
                event-card
                event-{event['level'].lower()}
            ">
                <div class="event-icon">
                    {icon}
                </div>

                <div>
                    <div class="event-top">
                        <span
                            class="
                                badge
                                {event['level'].lower()}
                            "
                        >
                            {
                                html.escape(
                                    event[
                                        "level"
                                    ]
                                )
                            }
                        </span>

                        <code>
                            {
                                html.escape(
                                    event[
                                        "code"
                                    ]
                                )
                            }
                        </code>
                    </div>

                    <div class="event-message">
                        {
                            html.escape(
                                event[
                                    "message"
                                ]
                            )
                        }
                    </div>
                </div>
            </div>
            """
        )

    # --------------------------------------------------------------
    # Module table
    # --------------------------------------------------------------

    module_names = sorted(
        set(
            module
            for read in reads
            for module in read[
                "modules"
            ]
        )
    )

    module_rows = []

    for module in module_names:
        cells = []

        for read in reads:
            status = read[
                "modules"
            ].get(
                module,
                "NA"
            )

            cells.append(
                f"""
                <td>
                    <span
                        class="
                            badge
                            {status.lower()}
                        "
                    >
                        {
                            html.escape(
                                status
                            )
                        }
                    </span>
                </td>
                """
            )

        module_rows.append(
            f"""
            <tr>
                <td>
                    {
                        html.escape(
                            module
                        )
                    }
                </td>
                {''.join(cells)}
            </tr>
            """
        )

    generated = datetime.now(
        timezone.utc
    ).isoformat()

    warning_text = ""

    if gate == "WARN":
        warning_text = """
        <div class="interpretation interpretation-warn">
            <div class="interpretation-symbol">!</div>
            <div>
                <strong>
                    Proceed with caution — not a hard failure
                </strong>
                <p>
                    One or more non-critical raw sequencing
                    QC signals require review. A WARN gate
                    does not automatically prevent alignment
                    or variant calling. Alignment, coverage,
                    duplication, contamination and
                    variant-level QC remain required before
                    judging overall sample suitability.
                </p>
            </div>
        </div>
        """

    elif gate == "FAIL":
        warning_text = """
        <div class="interpretation interpretation-fail">
            <div class="interpretation-symbol">×</div>
            <div>
                <strong>
                    Raw sequencing QC failed
                </strong>
                <p>
                    At least one configured critical QC
                    condition failed. Review the QC events
                    before downstream processing.
                </p>
            </div>
        </div>
        """

    else:
        warning_text = """
        <div class="interpretation interpretation-pass">
            <div class="interpretation-symbol">✓</div>
            <div>
                <strong>
                    Raw sequencing QC passed
                </strong>
                <p>
                    No configured raw sequencing QC issue
                    requiring escalation was detected.
                    Alignment-level QC is still required.
                </p>
            </div>
        </div>
        """

    html_text = f"""
<!DOCTYPE html>
<html lang="en">

<head>

<meta charset="UTF-8">

<meta
    name="viewport"
    content="width=device-width, initial-scale=1.0"
>

<meta
    name="author"
    content="{html.escape(AUTHOR)}"
>

<title>
    {html.escape(sample_id)}
    · Raw Sequencing QC
</title>

<style>

:root {{
    --bg: #07101f;
    --bg2: #0a1325;

    --panel: rgba(
        16,
        27,
        48,
        0.92
    );

    --panel2: rgba(
        21,
        35,
        61,
        0.92
    );

    --border: rgba(
        148,
        163,
        184,
        0.14
    );

    --text: #eef4ff;
    --muted: #91a3bd;

    --pass: #27d59b;
    --warn: #ffb648;
    --fail: #ff6978;

    --r1: #68a7ff;
    --r2: #b78cff;

    --accent: #79aaff;
}}


* {{
    box-sizing: border-box;
}}


html {{
    scroll-behavior: smooth;
}}


body {{
    margin: 0;

    color: var(--text);

    font-family:
        Inter,
        ui-sans-serif,
        system-ui,
        -apple-system,
        BlinkMacSystemFont,
        "Segoe UI",
        sans-serif;

    background:
        radial-gradient(
            circle at 80% -10%,
            rgba(
                82,
                113,
                255,
                0.20
            ),
            transparent 36%
        ),
        radial-gradient(
            circle at -10% 28%,
            rgba(
                83,
                205,
                178,
                0.08
            ),
            transparent 30%
        ),
        linear-gradient(
            180deg,
            var(--bg),
            var(--bg2)
        );

    min-height: 100vh;
}}


.page {{
    max-width: 1380px;

    margin:
        0
        auto;

    padding:
        42px
        28px
        70px;
}}


/* --------------------------------------------------------------
   Header
-------------------------------------------------------------- */

.hero {{
    position: relative;

    overflow: hidden;

    border:
        1px
        solid
        var(--border);

    border-radius: 24px;

    padding:
        30px
        32px;

    margin-bottom: 20px;

    background:
        linear-gradient(
            115deg,
            rgba(
                26,
                44,
                77,
                0.86
            ),
            rgba(
                10,
                22,
                42,
                0.90
            )
        );

    box-shadow:
        0
        22px
        70px
        rgba(
            0,
            0,
            0,
            0.20
        );
}}


.hero:after {{
    content: "";

    position: absolute;

    width: 330px;
    height: 330px;

    border-radius: 50%;

    top: -190px;
    right: -70px;

    background:
        rgba(
            104,
            167,
            255,
            0.12
        );

    filter: blur(2px);
}}


.hero-top {{
    display: flex;

    align-items: center;

    justify-content:
        space-between;

    gap: 18px;
}}


.product {{
    font-size: 12px;

    letter-spacing:
        0.15em;

    text-transform:
        uppercase;

    color:
        var(--accent);

    font-weight: 800;
}}


.hero h1 {{
    font-size:
        clamp(
            31px,
            5vw,
            48px
        );

    letter-spacing:
        -0.035em;

    margin:
        9px
        0
        6px;
}}


.hero p {{
    color:
        var(--muted);

    margin: 0;

    font-size: 15px;
}}


.sample-name {{
    color:
        var(--text);

    font-weight: 800;
}}


/* --------------------------------------------------------------
   Cards
-------------------------------------------------------------- */

.metric-grid {{
    display: grid;

    grid-template-columns:
        repeat(
            5,
            minmax(
                0,
                1fr
            )
        );

    gap: 14px;

    margin-bottom: 20px;
}}


.metric {{
    padding:
        18px
        19px;

    border:
        1px
        solid
        var(--border);

    border-radius: 18px;

    background:
        linear-gradient(
            180deg,
            rgba(
                255,
                255,
                255,
                0.025
            ),
            rgba(
                255,
                255,
                255,
                0
            )
        ),
        var(--panel);

    min-height: 113px;
}}


.metric-label {{
    color:
        var(--muted);

    font-size: 12px;

    font-weight: 700;

    text-transform:
        uppercase;

    letter-spacing:
        0.055em;
}}


.metric-value {{
    font-size: 27px;

    font-weight: 800;

    letter-spacing:
        -0.025em;

    margin-top: 12px;
}}


.metric-sub {{
    color:
        var(--muted);

    margin-top: 3px;

    font-size: 11px;
}}


/* --------------------------------------------------------------
   Panels
-------------------------------------------------------------- */

.panel {{
    border:
        1px
        solid
        var(--border);

    background:
        var(--panel);

    border-radius: 20px;

    padding:
        22px;

    margin-bottom: 20px;

    box-shadow:
        0
        17px
        50px
        rgba(
            0,
            0,
            0,
            0.13
        );
}}


.panel-title {{
    display: flex;

    align-items: center;

    justify-content:
        space-between;

    gap: 14px;

    margin-bottom: 14px;
}}


.panel h2 {{
    font-size: 19px;

    margin: 0;

    letter-spacing:
        -0.015em;
}}


.panel-sub {{
    color:
        var(--muted);

    font-size: 12px;
}}


.chart-grid {{
    display: grid;

    grid-template-columns:
        2fr
        1fr;

    gap: 20px;

    margin-bottom: 20px;
}}


.chart-grid .panel {{
    margin-bottom: 0;
}}


/* --------------------------------------------------------------
   Badges
-------------------------------------------------------------- */

.badge {{
    display:
        inline-flex;

    align-items:
        center;

    justify-content:
        center;

    min-width: 63px;

    padding:
        5px
        10px;

    border-radius:
        999px;

    font-size:
        11px;

    font-weight:
        900;

    letter-spacing:
        0.035em;
}}


.pass {{
    color:
        #7bf0c5;

    background:
        rgba(
            39,
            213,
            155,
            0.12
        );

    border:
        1px
        solid
        rgba(
            39,
            213,
            155,
            0.28
        );
}}


.warn {{
    color:
        #ffc86e;

    background:
        rgba(
            255,
            182,
            72,
            0.12
        );

    border:
        1px
        solid
        rgba(
            255,
            182,
            72,
            0.29
        );
}}


.fail {{
    color:
        #ff9ba7;

    background:
        rgba(
            255,
            105,
            120,
            0.12
        );

    border:
        1px
        solid
        rgba(
            255,
            105,
            120,
            0.29
        );
}}


.na {{
    color:
        #bac6d7;

    background:
        rgba(
            148,
            163,
            184,
            0.10
        );

    border:
        1px
        solid
        rgba(
            148,
            163,
            184,
            0.20
        );
}}


.hero-gate {{
    position: relative;

    z-index: 3;

    font-size: 13px;

    padding:
        8px
        16px;
}}


/* --------------------------------------------------------------
   Interpretation
-------------------------------------------------------------- */

.interpretation {{
    display: flex;

    gap: 14px;

    padding:
        16px
        18px;

    border-radius: 16px;

    margin-bottom: 20px;

    border:
        1px
        solid
        var(--border);
}}


.interpretation-symbol {{
    width: 32px;
    height: 32px;

    flex:
        0 0 32px;

    border-radius:
        10px;

    display:
        flex;

    align-items:
        center;

    justify-content:
        center;

    font-weight:
        900;
}}


.interpretation strong {{
    font-size: 14px;
}}


.interpretation p {{
    margin:
        5px
        0
        0;

    color:
        var(--muted);

    font-size: 12px;

    line-height: 1.6;
}}


.interpretation-warn {{
    background:
        rgba(
            255,
            182,
            72,
            0.055
        );
}}


.interpretation-warn
.interpretation-symbol {{
    color:
        var(--warn);

    background:
        rgba(
            255,
            182,
            72,
            0.13
        );
}}


.interpretation-pass {{
    background:
        rgba(
            39,
            213,
            155,
            0.05
        );
}}


.interpretation-pass
.interpretation-symbol {{
    color:
        var(--pass);

    background:
        rgba(
            39,
            213,
            155,
            0.12
        );
}}


.interpretation-fail {{
    background:
        rgba(
            255,
            105,
            120,
            0.05
        );
}}


.interpretation-fail
.interpretation-symbol {{
    color:
        var(--fail);

    background:
        rgba(
            255,
            105,
            120,
            0.12
        );
}}


/* --------------------------------------------------------------
   SVG charts
-------------------------------------------------------------- */

.chart-header {{
    display: flex;

    align-items:
        flex-start;

    justify-content:
        space-between;

    gap: 18px;

    margin-bottom: 6px;
}}


.chart-header h3 {{
    margin:
        0
        0
        4px;

    font-size: 14px;
}}


.legend {{
    display: flex;

    flex-wrap: wrap;

    gap: 14px;

    color:
        var(--muted);

    font-size: 11px;
}}


.legend-item {{
    display:
        inline-flex;

    align-items:
        center;

    gap: 6px;
}}


.legend-line {{
    display:
        inline-block;

    width: 18px;
    height: 3px;

    border-radius:
        999px;
}}


.chart-scroll {{
    overflow-x: auto;
}}


.chart-svg {{
    display: block;

    width: 100%;

    min-width: 590px;
}}


.grid-line {{
    stroke:
        rgba(
            148,
            163,
            184,
            0.10
        );

    stroke-width: 1;
}}


.axis-line {{
    stroke:
        rgba(
            148,
            163,
            184,
            0.26
        );

    stroke-width: 1;
}}


.axis-text {{
    fill:
        var(--muted);

    font-size: 10px;
}}


.axis-label {{
    fill:
        #a9b7cb;

    font-size: 11px;

    font-weight: 700;
}}


.chart-line {{
    fill: none;

    stroke-width: 2.5;

    vector-effect:
        non-scaling-stroke;

    stroke-linecap: round;

    stroke-linejoin: round;
}}


.line-r1 {{
    stroke:
        var(--r1);

    background:
        var(--r1);
}}


.line-r2 {{
    stroke:
        var(--r2);

    background:
        var(--r2);
}}


.hover-point {{
    opacity: 0;

    cursor:
        crosshair;
}}


.hover-point:hover {{
    opacity: 0.88;
}}


.line-r1-point {{
    fill:
        var(--r1);
}}


.line-r2-point {{
    fill:
        var(--r2);
}}


.empty-chart {{
    color:
        var(--muted);

    padding:
        48px
        20px;

    text-align: center;

    border:
        1px
        dashed
        var(--border);

    border-radius: 13px;
}}


/* --------------------------------------------------------------
   Donut
-------------------------------------------------------------- */

.donut-layout {{
    display: flex;

    align-items: center;

    justify-content:
        center;

    gap: 26px;

    min-height: 240px;
}}


.donut {{
    --p1:
        calc(
            var(--pass)
            * 1%
        );

    --p2:
        calc(
            (
                var(--pass)
                + var(--warn)
            )
            * 1%
        );

    width: 148px;
    height: 148px;

    border-radius: 50%;

    position: relative;

    background:
        conic-gradient(
            var(--pass-color)
            0
            var(--p1),

            var(--warn-color)
            var(--p1)
            var(--p2),

            var(--fail-color)
            var(--p2)
            100%
        );

    --pass-color:
        #27d59b;

    --warn-color:
        #ffb648;

    --fail-color:
        #ff6978;
}}


.donut:after {{
    content: "";

    position: absolute;

    inset: 18px;

    border-radius: 50%;

    background:
        var(--panel2);
}}


.donut-center {{
    position: absolute;

    z-index: 2;

    inset: 0;

    display: flex;

    flex-direction: column;

    align-items: center;

    justify-content:
        center;
}}


.donut-center strong {{
    font-size: 29px;
}}


.donut-center span {{
    color:
        var(--muted);

    font-size: 10px;

    text-transform:
        uppercase;

    letter-spacing:
        0.08em;
}}


.donut-legend {{
    display: grid;

    gap: 11px;

    min-width: 115px;

    font-size: 12px;
}}


.donut-legend div {{
    display: grid;

    grid-template-columns:
        10px
        1fr
        auto;

    gap: 8px;

    align-items: center;

    color:
        var(--muted);
}}


.donut-legend strong {{
    color:
        var(--text);
}}


.dot {{
    width: 8px;
    height: 8px;

    border-radius: 50%;
}}


.dot-pass {{
    background:
        var(--pass);
}}


.dot-warn {{
    background:
        var(--warn);
}}


.dot-fail {{
    background:
        var(--fail);
}}


/* --------------------------------------------------------------
   Event cards
-------------------------------------------------------------- */

.events {{
    display: grid;

    gap: 10px;
}}


.event-card {{
    display: flex;

    align-items:
        flex-start;

    gap: 13px;

    padding:
        13px
        14px;

    border-radius: 14px;

    background:
        rgba(
            255,
            255,
            255,
            0.018
        );

    border:
        1px
        solid
        var(--border);
}}


.event-icon {{
    width: 29px;
    height: 29px;

    flex:
        0 0 29px;

    display: flex;

    align-items: center;

    justify-content:
        center;

    border-radius: 8px;

    font-weight: 900;
}}


.event-warn
.event-icon {{
    color:
        var(--warn);

    background:
        rgba(
            255,
            182,
            72,
            0.12
        );
}}


.event-fail
.event-icon {{
    color:
        var(--fail);

    background:
        rgba(
            255,
            105,
            120,
            0.12
        );
}}


.event-pass
.event-icon {{
    color:
        var(--pass);

    background:
        rgba(
            39,
            213,
            155,
            0.12
        );
}}


.event-top {{
    display: flex;

    gap: 9px;

    align-items: center;

    flex-wrap: wrap;
}}


.event-message {{
    color:
        var(--muted);

    font-size: 12px;

    margin-top: 6px;
}}


code {{
    color:
        #b9cced;

    font-size: 11px;
}}


/* --------------------------------------------------------------
   Tables
-------------------------------------------------------------- */

.table-wrap {{
    overflow-x: auto;
}}


table {{
    width: 100%;

    border-collapse:
        collapse;

    font-size: 12px;
}}


th {{
    text-align: left;

    color:
        var(--muted);

    font-size: 10px;

    text-transform:
        uppercase;

    letter-spacing:
        0.06em;

    font-weight: 800;
}}


th,
td {{
    padding:
        11px
        10px;

    border-bottom:
        1px
        solid
        var(--border);
}}


tbody tr {{
    transition:
        background
        0.15s
        ease;
}}


tbody tr:hover {{
    background:
        rgba(
            255,
            255,
            255,
            0.025
        );
}}


.filename {{
    color:
        #c4d0e2;

    font-family:
        ui-monospace,
        SFMono-Regular,
        Menlo,
        monospace;

    font-size: 10px;
}}


.read-pill {{
    display:
        inline-flex;

    padding:
        4px
        8px;

    border-radius:
        7px;

    color:
        #dbe9ff;

    background:
        rgba(
            104,
            167,
            255,
            0.14
        );

    font-weight: 900;
}}


/* --------------------------------------------------------------
   Collapsible details
-------------------------------------------------------------- */

details {{
    border:
        1px
        solid
        var(--border);

    border-radius:
        18px;

    margin-bottom: 20px;

    background:
        var(--panel);
}}


summary {{
    cursor: pointer;

    list-style: none;

    padding:
        18px
        21px;

    font-weight: 800;

    font-size: 14px;

    user-select: none;
}}


summary::-webkit-details-marker {{
    display: none;
}}


summary:after {{
    content: "+";

    float: right;

    color:
        var(--accent);

    font-size: 20px;
}}


details[open]
summary:after {{
    content: "–";
}}


.detail-content {{
    padding:
        0
        21px
        21px;
}}


/* --------------------------------------------------------------
   Footer
-------------------------------------------------------------- */

.footer {{
    margin-top: 33px;

    display: flex;

    align-items:
        flex-end;

    justify-content:
        space-between;

    gap: 20px;

    padding-top: 22px;

    border-top:
        1px
        solid
        var(--border);

    color:
        var(--muted);

    font-size: 11px;

    line-height: 1.65;
}}


.author {{
    color:
        #d8e4f6;

    font-weight: 800;
}}


.disclaimer {{
    max-width: 680px;
}}


/* --------------------------------------------------------------
   Responsive
-------------------------------------------------------------- */

@media (
    max-width: 1050px
) {{

    .metric-grid {{
        grid-template-columns:
            repeat(
                2,
                minmax(
                    0,
                    1fr
                )
            );
    }}

    .chart-grid {{
        grid-template-columns:
            1fr;
    }}

}}


@media (
    max-width: 620px
) {{

    .page {{
        padding:
            20px
            13px
            45px;
    }}

    .hero {{
        padding:
            23px
            20px;
    }}

    .hero-top {{
        align-items:
            flex-start;

        flex-direction:
            column;
    }}

    .metric-grid {{
        grid-template-columns:
            1fr;
    }}

    .footer {{
        flex-direction:
            column;

        align-items:
            flex-start;
    }}

    .donut-layout {{
        flex-direction:
            column;
    }}

}}

</style>

</head>


<body>

<div class="page">


<section class="hero">

    <div class="hero-top">

        <div>

            <div class="product">
                Modular Human Genomics Pipeline
            </div>

            <h1>
                Raw Sequencing QC
            </h1>

            <p>
                Sample
                <span class="sample-name">
                    {html.escape(sample_id)}
                </span>
                · Illumina paired-end WGS
            </p>

        </div>

        <span
            class="
                badge
                hero-gate
                {gate.lower()}
            "
        >
            QC {html.escape(gate)}
        </span>

    </div>

</section>


<section class="metric-grid">

    <div class="metric">

        <div class="metric-label">
            QC gate
        </div>

        <div class="metric-value">
            {html.escape(gate)}
        </div>

        <div class="metric-sub">
            Raw sequencing gate
        </div>

    </div>


    <div class="metric">

        <div class="metric-label">
            Read pairs
        </div>

        <div class="metric-value">
            {pair_count}
        </div>

        <div class="metric-sub">
            Matched R1/R2 pairs
        </div>

    </div>


    <div class="metric">

        <div class="metric-label">
            Reads / mate
        </div>

        <div class="metric-value">
            {read_count}
        </div>

        <div class="metric-sub">
            Sequencing depth input
        </div>

    </div>


    <div class="metric">

        <div class="metric-label">
            Read length
        </div>

        <div class="metric-value">
            {html.escape(read_length)}
            <span style="
                font-size:13px;
                color:var(--muted);
            ">
                bp
            </span>
        </div>

        <div class="metric-sub">
            FastQC basic statistics
        </div>

    </div>


    <div class="metric">

        <div class="metric-label">
            GC content
        </div>

        <div class="metric-value">
            {html.escape(gc_percent)}
        </div>

        <div class="metric-sub">
            Overall sequence GC
        </div>

    </div>

</section>


{warning_text}


<section class="chart-grid">

    <div class="panel">

        <div class="panel-title">

            <div>
                <h2>
                    Base quality profile
                </h2>

                <div class="panel-sub">
                    Actual FastQC per-base
                    mean Phred values
                </div>
            </div>

        </div>

        {quality_chart}

    </div>


    <div class="panel">

        <div class="panel-title">

            <div>
                <h2>
                    FastQC status
                </h2>

                <div class="panel-sub">
                    Across R1 and R2 modules
                </div>
            </div>

        </div>

        {status_donut(module_counts)}

    </div>

</section>


<section class="panel">

    <div class="panel-title">

        <div>

            <h2>
                GC distribution
            </h2>

            <div class="panel-sub">
                Actual per-sequence GC
                distribution from FastQC.
                Hover over the curves
                for individual values.
            </div>

        </div>

    </div>

    {gc_chart}

</section>


<section class="panel">

    <div class="panel-title">

        <div>

            <h2>
                Adapter signal
            </h2>

            <div class="panel-sub">
                Maximum detected adapter
                content at each read position
            </div>

        </div>

    </div>

    {adapter_chart}

</section>


<section class="panel">

    <div class="panel-title">

        <div>

            <h2>
                QC observations
            </h2>

            <div class="panel-sub">
                Rule-based findings generated
                by QC Engine v{QC_ENGINE_VERSION}
            </div>

        </div>

    </div>

    <div class="events">
        {''.join(event_cards)}
    </div>

</section>


<details>

    <summary>
        Read-level technical statistics
    </summary>

    <div class="detail-content">

        <div class="table-wrap">

            <table>

                <thead>
                    <tr>
                        <th>Read</th>
                        <th>Filename</th>
                        <th>Sequences</th>
                        <th>Length</th>
                        <th>GC</th>
                        <th>Encoding</th>
                    </tr>
                </thead>

                <tbody>
                    {''.join(read_rows)}
                </tbody>

            </table>

        </div>

    </div>

</details>


<details>

    <summary>
        Complete FastQC module results
    </summary>

    <div class="detail-content">

        <div class="table-wrap">

            <table>

                <thead>

                    <tr>

                        <th>
                            Module
                        </th>

                        {
                            ''.join(
                                f"<th>"
                                f"{html.escape(str(read['read']))}"
                                f"</th>"
                                for read in reads
                            )
                        }

                    </tr>

                </thead>

                <tbody>
                    {''.join(module_rows)}
                </tbody>

            </table>

        </div>

    </div>

</details>


<footer class="footer">

    <div>

        <div class="author">
            Pipeline & QC dashboard:
            {html.escape(AUTHOR)}
        </div>

        <div>
            Modular Human Genomics Pipeline
        </div>

        <div>
            QC Engine v{QC_ENGINE_VERSION}
            · Schema v{REPORT_SCHEMA_VERSION}
        </div>

        <div>
            Generated UTC:
            {html.escape(generated)}
        </div>

    </div>


    <div class="disclaimer">
        Research, benchmarking and pipeline
        validation use. Raw sequencing QC
        does not constitute clinical
        interpretation or diagnostic
        validation. Downstream alignment,
        coverage, sample integrity,
        contamination and variant-level QC
        are required.
    </div>

</footer>


</div>

</body>

</html>
"""

    output_path.write_text(
        html_text,
        encoding="utf-8",
    )


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Parse FastQC output and "
            "generate structured raw "
            "sequencing QC."
        )
    )

    parser.add_argument(
        "--sample-id",
        required=True,
    )

    parser.add_argument(
        "--fastqc-zip",
        nargs="+",
        required=True,
    )

    parser.add_argument(
        "--thresholds",
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

    setup_logging(
        Path(args.log)
    )

    try:
        LOGGER.info(
            "Starting raw QC for sample %s",
            args.sample_id,
        )

        threshold_path = Path(
            args.thresholds
        )

        if not threshold_path.exists():
            raise FileNotFoundError(
                f"Threshold file not found: "
                f"{threshold_path}"
            )

        with threshold_path.open(
            encoding="utf-8"
        ) as handle:
            thresholds = yaml.safe_load(
                handle
            )

        if not thresholds:
            raise ValueError(
                "Threshold configuration is empty"
            )

        reads = [
            parse_archive(
                Path(archive)
            )
            for archive
            in args.fastqc_zip
        ]

        reads.sort(
            key=lambda item: item[
                "read"
            ]
        )

        gate, events = evaluate_qc(
            reads,
            thresholds,
        )

        result = {
            "schema_version":
                REPORT_SCHEMA_VERSION,

            "qc_engine_version":
                QC_ENGINE_VERSION,

            "author":
                AUTHOR,

            "sample_id":
                args.sample_id,

            "generated_at":
                datetime.now(
                    timezone.utc
                ).isoformat(),

            "gate":
                gate,

            "reads":
                reads,

            "events":
                events,

            "module_summary":
                module_status_counts(
                    reads
                ),

            "event_summary":
                event_status_counts(
                    events
                ),
        }

        output_json = Path(
            args.output_json
        )

        output_json.write_text(
            json.dumps(
                result,
                indent=2,
                ensure_ascii=False,
            ),
            encoding="utf-8",
        )

        build_html(
            args.sample_id,
            gate,
            reads,
            events,
            Path(
                args.output_html
            ),
        )

        for event in events:

            if event["level"] == "FAIL":
                LOGGER.error(
                    "%s | %s",
                    event["code"],
                    event["message"],
                )

            elif event["level"] == "WARN":
                LOGGER.warning(
                    "%s | %s",
                    event["code"],
                    event["message"],
                )

            else:
                LOGGER.info(
                    "%s | %s",
                    event["code"],
                    event["message"],
                )

        LOGGER.info(
            "Raw QC completed. "
            "Gate=%s | Engine=%s | "
            "Author=%s",
            gate,
            QC_ENGINE_VERSION,
            AUTHOR,
        )

    except Exception as exc:
        LOGGER.exception(
            "RAW_QC_FATAL_ERROR: %s",
            exc,
        )

        sys.exit(1)


if __name__ == "__main__":
    main()
