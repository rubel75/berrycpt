"""Read and validate BerryCPT ``bcurv_ij.dat`` files.

Two BerryCPT output formats are supported.

Current format
--------------

A k-point header has the form

    # KP: <index> NBANDS_OUT: <reported bands> NBANDS_TRANS: <transition bands>

and each band-resolved row contains

    band_index  degenerate_block  Omega_yz  Omega_zx  Omega_xy

The occupation-weighted total is written as

    0  0  Omega_yz  Omega_zx  Omega_xy

Legacy format
-------------

A k-point header has the form

    # KP: <index> NVB: <reported bands> NEMAX: <transition bands>

and each band-resolved row contains

    band_index  Omega_yz  Omega_zx  Omega_xy

The legacy format has no degenerate-block column. Its ``NVB`` and ``NEMAX``
metadata are mapped to the current ``NBANDS_OUT`` and ``NBANDS_TRANS``
concepts, respectively.

For both formats, the row with band index 0 is the occupation-weighted total
written by BerryCPT. Since Berry curvature is an axial vector, the three
tensor columns are also the Cartesian vector components

    (Omega_x, Omega_y, Omega_z) = (Omega_yz, Omega_zx, Omega_xy).

This module only parses and validates the irreducible-point data. It does not
perform symmetry reconstruction or Brillouin-zone integration.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from pathlib import Path
import re
from typing import Iterator, Literal, Sequence

import numpy as np
from numpy.typing import NDArray

from vasp_ibzkpt import IBZ, read_ibzkpt


FloatVector = NDArray[np.float64]
BerryCPTFormat = Literal["legacy", "current"]


class BerryCurvatureFileError(ValueError):
    """Raised when a BerryCPT Berry-curvature file is malformed."""


_LEGACY_HEADER_RE = re.compile(
    r"^#\s*KP\s*:\s*(\d+)\s+NVB\s*:\s*(\d+)\s+NEMAX\s*:\s*(\d+)\s*$",
    flags=re.IGNORECASE,
)

_CURRENT_HEADER_RE = re.compile(
    r"^#\s*KP\s*:\s*(\d+)\s+NBANDS_OUT\s*:\s*(\d+)\s+"
    r"NBANDS_TRANS\s*:\s*(\d+)\s*$",
    flags=re.IGNORECASE,
)


def _readonly_vector(values: Sequence[float]) -> FloatVector:
    vector = np.asarray(values, dtype=np.float64)
    if vector.shape != (3,):
        raise BerryCurvatureFileError(
            f"Berry-curvature vector must have shape (3,), got {vector.shape}"
        )
    if not np.all(np.isfinite(vector)):
        raise BerryCurvatureFileError(
            "Berry-curvature vector contains non-finite values"
        )
    vector.setflags(write=False)
    return vector


@dataclass(frozen=True)
class BandBerryCurvature:
    """One band-resolved row from ``bcurv_ij.dat``.

    ``degenerate_block`` is ``None`` for legacy files, which did not contain
    this column. In the current format, positive bands have positive local
    block labels and the occupation-weighted total has block label 0.
    """

    band_index: int
    omega_xyz: FloatVector
    source_line: int
    degenerate_block: int | None = None

    @property
    def omega_yz_zx_xy(self) -> FloatVector:
        """Return the three tensor columns in their original order."""

        return self.omega_xyz


@dataclass(frozen=True)
class KPointBerryCurvature:
    """One complete irreducible-k-point block.

    The stored names ``nvb`` and ``nemax`` are retained for compatibility with
    the first utility version. Their format-independent aliases are
    ``nbands_out`` and ``nbands_trans``.
    """

    kp_index: int
    nvb: int
    nemax: int
    rows: tuple[BandBerryCurvature, ...]
    header_line: int
    format_version: BerryCPTFormat = "legacy"

    @property
    def nbands_out(self) -> int:
        """Number of positive-band rows written for this k-point."""

        return self.nvb

    @property
    def nbands_trans(self) -> int:
        """Number of transition-space bands used by BerryCPT."""

        return self.nemax

    @property
    def band_rows(self) -> tuple[BandBerryCurvature, ...]:
        return tuple(row for row in self.rows if row.band_index > 0)

    @property
    def total_row(self) -> BandBerryCurvature:
        totals = tuple(row for row in self.rows if row.band_index == 0)
        if len(totals) != 1:
            raise BerryCurvatureFileError(
                f"k-point {self.kp_index}: expected one band-0 total, "
                f"found {len(totals)}"
            )
        return totals[0]

    @property
    def total_omega_xyz(self) -> FloatVector:
        """Occupation-weighted total Berry curvature written by BerryCPT."""

        return self.total_row.omega_xyz

    def band(self, band_index: int) -> BandBerryCurvature:
        """Return one positive-band record by its 1-based band index."""

        if band_index <= 0:
            raise KeyError("band() accepts only positive band indices")
        for row in self.rows:
            if row.band_index == band_index:
                return row
        raise KeyError(
            f"k-point {self.kp_index}: band {band_index} is not present"
        )


@dataclass(frozen=True)
class BerryCurvatureData:
    """Validated contents of one BerryCPT Berry-curvature file."""

    records: tuple[KPointBerryCurvature, ...]
    comments: tuple[str, ...]
    source: Path

    def __len__(self) -> int:
        return len(self.records)

    def __iter__(self) -> Iterator[KPointBerryCurvature]:
        return iter(self.records)

    def __getitem__(self, index: int) -> KPointBerryCurvature:
        return self.records[index]

    @property
    def nbands_out_distribution(self) -> dict[int, int]:
        return dict(
            sorted(Counter(record.nbands_out for record in self.records).items())
        )

    @property
    def nbands_trans_distribution(self) -> dict[int, int]:
        return dict(
            sorted(Counter(record.nbands_trans for record in self.records).items())
        )

    @property
    def nvb_distribution(self) -> dict[int, int]:
        """Legacy alias for ``nbands_out_distribution``."""

        return self.nbands_out_distribution

    @property
    def nemax_distribution(self) -> dict[int, int]:
        """Legacy alias for ``nbands_trans_distribution``."""

        return self.nbands_trans_distribution

    @property
    def format_distribution(self) -> dict[str, int]:
        return dict(
            sorted(Counter(record.format_version for record in self.records).items())
        )

    def totals_xyz(self) -> NDArray[np.float64]:
        """Return an ``(nk, 3)`` array of occupation-weighted totals."""

        if not self.records:
            return np.empty((0, 3), dtype=np.float64)
        return np.stack(
            [record.total_omega_xyz for record in self.records], axis=0
        )


@dataclass(frozen=True)
class IBZBerryCurvature:
    """An IBZKPT object paired index-by-index with BerryCPT records."""

    ibz: IBZ
    bcurv: BerryCurvatureData

    def __len__(self) -> int:
        return len(self.bcurv)


def _parse_nonnegative_integer(token: str, context: str) -> int:
    try:
        value = int(token)
    except ValueError as exc:
        raise BerryCurvatureFileError(
            f"{context}: invalid integer {token!r}"
        ) from exc
    if value < 0:
        raise BerryCurvatureFileError(
            f"{context}: expected a nonnegative integer"
        )
    return value


def _parse_vector(fields: Sequence[str], line_number: int) -> FloatVector:
    if len(fields) != 3:
        raise BerryCurvatureFileError(
            f"line {line_number}: expected three Berry-curvature components"
        )
    try:
        values = [
            float(field.replace("D", "E").replace("d", "e"))
            for field in fields
        ]
    except ValueError as exc:
        raise BerryCurvatureFileError(
            f"line {line_number}: invalid Berry-curvature component"
        ) from exc
    return _readonly_vector(values)


def _validate_degenerate_blocks(
    kp_index: int,
    rows: list[BandBerryCurvature],
    format_version: BerryCPTFormat,
) -> None:
    if format_version == "legacy":
        if any(row.degenerate_block is not None for row in rows):
            raise BerryCurvatureFileError(
                f"k-point {kp_index}: legacy rows cannot contain block labels"
            )
        return

    if any(row.degenerate_block is None for row in rows):
        raise BerryCurvatureFileError(
            f"k-point {kp_index}: current-format rows require block labels"
        )

    total_rows = [row for row in rows if row.band_index == 0]
    if total_rows and total_rows[0].degenerate_block != 0:
        raise BerryCurvatureFileError(
            f"k-point {kp_index}: band-0 total must have degenerate block 0"
        )

    positive_rows = sorted(
        (row for row in rows if row.band_index > 0),
        key=lambda row: row.band_index,
    )
    labels = [int(row.degenerate_block) for row in positive_rows]
    if any(label < 1 for label in labels):
        raise BerryCurvatureFileError(
            f"k-point {kp_index}: positive bands require positive block labels"
        )

    if labels:
        if labels[0] != 1:
            raise BerryCurvatureFileError(
                f"k-point {kp_index}: first degenerate block must be 1"
            )
        for previous, current in zip(labels, labels[1:]):
            if current < previous or current > previous + 1:
                raise BerryCurvatureFileError(
                    f"k-point {kp_index}: degenerate-block labels must be "
                    "contiguous and nondecreasing"
                )


def _validate_block(
    kp_index: int,
    nvb: int,
    nemax: int,
    rows: list[BandBerryCurvature],
    header_line: int,
    require_total_last: bool,
    format_version: BerryCPTFormat,
) -> KPointBerryCurvature:
    if kp_index < 1:
        raise BerryCurvatureFileError(
            f"line {header_line}: KP index must be positive"
        )
    if nvb < 1:
        raise BerryCurvatureFileError(
            f"k-point {kp_index}: NBANDS_OUT must be positive"
        )
    if nemax < nvb:
        raise BerryCurvatureFileError(
            f"k-point {kp_index}: NBANDS_TRANS={nemax} is smaller than "
            f"NBANDS_OUT={nvb}"
        )
    if not rows:
        raise BerryCurvatureFileError(
            f"k-point {kp_index}: block contains no data rows"
        )

    band_indices = [row.band_index for row in rows]
    negative = [band for band in band_indices if band < 0]
    if negative:
        raise BerryCurvatureFileError(
            f"k-point {kp_index}: negative band indices are not allowed: "
            f"{negative}"
        )

    duplicate_counts = Counter(band_indices)
    duplicates = sorted(
        band for band, count in duplicate_counts.items() if count > 1
    )
    if duplicates:
        raise BerryCurvatureFileError(
            f"k-point {kp_index}: duplicate band indices {duplicates}"
        )

    positive = sorted(band for band in band_indices if band > 0)
    expected = list(range(1, nvb + 1))
    if positive != expected:
        missing = sorted(set(expected) - set(positive))
        extra = sorted(set(positive) - set(expected))
        raise BerryCurvatureFileError(
            f"k-point {kp_index}: positive bands do not match "
            f"1..NBANDS_OUT. Missing={missing}, extra={extra}"
        )

    if band_indices.count(0) != 1:
        raise BerryCurvatureFileError(
            f"k-point {kp_index}: expected exactly one band-0 total row"
        )
    if require_total_last and band_indices[-1] != 0:
        raise BerryCurvatureFileError(
            f"k-point {kp_index}: band-0 total must be the final row"
        )

    _validate_degenerate_blocks(kp_index, rows, format_version)

    return KPointBerryCurvature(
        kp_index=kp_index,
        nvb=nvb,
        nemax=nemax,
        rows=tuple(rows),
        header_line=header_line,
        format_version=format_version,
    )


def _parse_header(
    text: str,
) -> tuple[int, int, int, BerryCPTFormat] | None:
    current_match = _CURRENT_HEADER_RE.match(text)
    if current_match:
        return (
            int(current_match.group(1)),
            int(current_match.group(2)),
            int(current_match.group(3)),
            "current",
        )

    legacy_match = _LEGACY_HEADER_RE.match(text)
    if legacy_match:
        return (
            int(legacy_match.group(1)),
            int(legacy_match.group(2)),
            int(legacy_match.group(3)),
            "legacy",
        )

    return None


def read_bcurv_ij(
    filename: str | Path,
    *,
    require_total_last: bool = True,
    require_sequential_kpoints: bool = True,
) -> BerryCurvatureData:
    """Read and validate a legacy or current BerryCPT output file.

    Parameters
    ----------
    filename
        Path to the BerryCPT output file.
    require_total_last
        Require the unique band-0 total row to be the final row of each block.
    require_sequential_kpoints
        Require KP indices to be exactly 1, 2, ..., ``nk`` in file order.
    """

    path = Path(filename)
    if not path.is_file():
        raise FileNotFoundError(f"Berry-curvature file not found: {path}")

    comments: list[str] = []
    records: list[KPointBerryCurvature] = []

    current_header: tuple[int, int, int, int, BerryCPTFormat] | None = None
    current_rows: list[BandBerryCurvature] = []

    def finalize_current() -> None:
        nonlocal current_header, current_rows
        if current_header is None:
            return
        kp_index, nvb, nemax, header_line, format_version = current_header
        records.append(
            _validate_block(
                kp_index=kp_index,
                nvb=nvb,
                nemax=nemax,
                rows=current_rows,
                header_line=header_line,
                require_total_last=require_total_last,
                format_version=format_version,
            )
        )
        current_header = None
        current_rows = []

    with path.open("r", encoding="utf-8", errors="replace") as stream:
        for line_number, raw_line in enumerate(stream, start=1):
            text = raw_line.strip()
            if not text:
                continue

            parsed_header = _parse_header(text)
            if parsed_header is not None:
                finalize_current()
                kp_index, nvb, nemax, format_version = parsed_header
                current_header = (
                    kp_index,
                    nvb,
                    nemax,
                    line_number,
                    format_version,
                )
                continue

            if text.startswith("#"):
                comments.append(text)
                continue

            if current_header is None:
                raise BerryCurvatureFileError(
                    f"line {line_number}: data row appears before the first "
                    "KP header"
                )

            fields = text.split()
            format_version = current_header[4]

            if format_version == "legacy":
                if len(fields) != 4:
                    raise BerryCurvatureFileError(
                        f"line {line_number}: legacy row requires a band "
                        "index and three components, "
                        f"found {len(fields)} fields"
                    )
                band_index = _parse_nonnegative_integer(
                    fields[0], f"line {line_number} band index"
                )
                degenerate_block = None
                vector_fields = fields[1:]
            else:
                if len(fields) != 5:
                    raise BerryCurvatureFileError(
                        f"line {line_number}: current row requires band and "
                        "degenerate-block indices plus three components, "
                        f"found {len(fields)} fields"
                    )
                band_index = _parse_nonnegative_integer(
                    fields[0], f"line {line_number} band index"
                )
                degenerate_block = _parse_nonnegative_integer(
                    fields[1],
                    f"line {line_number} degenerate-block index",
                )
                vector_fields = fields[2:]

            current_rows.append(
                BandBerryCurvature(
                    band_index=band_index,
                    omega_xyz=_parse_vector(vector_fields, line_number),
                    source_line=line_number,
                    degenerate_block=degenerate_block,
                )
            )

    finalize_current()

    if not records:
        raise BerryCurvatureFileError("file contains no KP blocks")

    kp_indices = [record.kp_index for record in records]
    duplicates = sorted(
        index for index, count in Counter(kp_indices).items() if count > 1
    )
    if duplicates:
        raise BerryCurvatureFileError(
            f"duplicate KP block indices {duplicates}"
        )

    if require_sequential_kpoints:
        expected = list(range(1, len(records) + 1))
        if kp_indices != expected:
            raise BerryCurvatureFileError(
                "KP blocks must appear sequentially as 1..nk. "
                f"Found first indices {kp_indices[:10]}"
            )

    return BerryCurvatureData(
        records=tuple(records),
        comments=tuple(comments),
        source=path,
    )


def pair_with_ibzkpt(
    bcurv: BerryCurvatureData,
    ibz: IBZ,
) -> IBZBerryCurvature:
    """Validate and pair BerryCPT blocks with IBZ points by KP index.

    ``bcurv_ij.dat`` contains no k-point coordinates, so the correspondence can
    only be checked by the sequential KP index and total record count. The
    files must therefore come from the same VASP/BerryCPT calculation and
    retain the original IBZ ordering.
    """

    if ibz.nktot != len(ibz.kpoints):
        raise BerryCurvatureFileError(
            f"IBZ object reports {ibz.nktot} points but contains "
            f"{len(ibz.kpoints)}"
        )
    if len(bcurv.records) != ibz.nktot:
        raise BerryCurvatureFileError(
            f"BerryCPT file has {len(bcurv.records)} KP blocks, but IBZKPT "
            f"has {ibz.nktot} points"
        )

    for expected_index, record in enumerate(bcurv.records, start=1):
        if record.kp_index != expected_index:
            raise BerryCurvatureFileError(
                f"BerryCPT record {expected_index} has KP index "
                f"{record.kp_index}"
            )

    return IBZBerryCurvature(ibz=ibz, bcurv=bcurv)


def load_ibz_bcurv(
    ibzkpt_filename: str | Path,
    bcurv_filename: str | Path,
) -> IBZBerryCurvature:
    """Read and pair matching ``IBZKPT`` and ``bcurv_ij.dat`` files."""

    ibz = read_ibzkpt(ibzkpt_filename)
    bcurv = read_bcurv_ij(bcurv_filename)
    return pair_with_ibzkpt(bcurv, ibz)


def _format_distribution(distribution: dict[int, int]) -> str:
    return ", ".join(
        f"{value}: {count}" for value, count in distribution.items()
    )


def _main() -> int:
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "Read and validate legacy or current BerryCPT band-resolved and "
            "band-0 total Berry curvature, optionally checking the KP "
            "ordering against IBZKPT."
        )
    )
    parser.add_argument("bcurv", help="path to bcurv_ij.dat")
    parser.add_argument(
        "--ibzkpt",
        help="matching IBZKPT file used to validate k-point count and ordering",
    )
    parser.add_argument(
        "--allow-total-not-last",
        action="store_true",
        help="accept a unique band-0 total row when it is not the final row",
    )
    args = parser.parse_args()

    data = read_bcurv_ij(
        args.bcurv,
        require_total_last=not args.allow_total_not_last,
    )

    print(f"BerryCPT k-point blocks: {len(data)}")
    print(f"Format distribution: {data.format_distribution}")
    print(
        "NBANDS_OUT distribution (bands: points): "
        f"{_format_distribution(data.nbands_out_distribution)}"
    )
    print(
        "NBANDS_TRANS distribution (bands: points): "
        f"{_format_distribution(data.nbands_trans_distribution)}"
    )
    print(f"Band-0 totals: {data.totals_xyz().shape[0]}")
    print(
        "Component convention: (Omega_x, Omega_y, Omega_z) = "
        "(Omega_yz, Omega_zx, Omega_xy)"
    )

    if args.ibzkpt:
        ibz = read_ibzkpt(args.ibzkpt)
        pair_with_ibzkpt(data, ibz)
        print(f"IBZKPT points: {ibz.nktot}")
        print("BerryCPT KP ordering matches IBZKPT index ordering")
        print(
            "Warning: bcurv_ij.dat has no k coordinates, so matching is "
            "index-based and assumes both files come from the same "
            "calculation."
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
