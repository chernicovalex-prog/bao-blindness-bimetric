#!/usr/bin/env python3
"""Build research-preview PDF/source/ZIP packages from the current text tree."""
from __future__ import annotations
import hashlib
import json
import os
from pathlib import Path
import subprocess
import zipfile

ROOT = Path(__file__).resolve().parent
OUT = ROOT.parent / "release"
OUT.mkdir(exist_ok=True)
ARTICLE = ROOT / "01_Statya"
PDF = ARTICLE / "main.pdf"


def sha_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha_file(path: Path) -> str:
    return sha_bytes(path.read_bytes())


def write_zip(path: Path, files: dict[str, bytes]) -> None:
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as archive:
        for name, data in sorted(files.items()):
            archive.writestr(name, data)


def source_files() -> dict[str, bytes]:
    excluded_names = {"main.aux", "main.blg", "main.log", "main.out", "main.txt"}
    result = {}
    for path in sorted(ROOT.rglob("*")):
        if not path.is_file():
            continue
        relative = path.relative_to(ROOT)
        if "__pycache__" in relative.parts or path.suffix == ".pyc":
            continue
        if path.name in excluded_names or path.name.endswith(".tmp"):
            continue
        if path.name == "MANIFEST.sha256":
            continue
        result[relative.as_posix()] = path.read_bytes()
    manifest = "".join(f"{sha_bytes(data)}  {name}\n" for name, data in sorted(result.items()))
    result["MANIFEST.sha256"] = manifest.encode()
    return result


def main() -> None:
    if not PDF.is_file():
        raise RuntimeError("main.pdf is missing; build the article first")
    required = ["results/summary.json", "results/anchor_scan.json",
                "results/cmb_lambda_theta_sweep.json", "results/scalar_k_map.json",
                "results/scalar_mfp_scan.json"]
    for name in required:
        if not (ROOT / name).is_file():
            raise RuntimeError(f"missing required result: {name}")
    pages_line = subprocess.check_output(["pdfinfo", str(PDF)], text=True)
    pages = int(next(line.split()[1] for line in pages_line.splitlines()
                     if line.startswith("Pages:")))
    figure_text = (ROOT / "results/figure_checks.txt").read_text()
    if "ALL FIGURE CHECKS PASSED (117)" not in figure_text:
        raise RuntimeError("figure assertion completion marker is absent")

    arxiv = {}
    for name in ["main.tex", "refs.bib", "main.bbl"]:
        arxiv[name] = (ARTICLE / name).read_bytes()
    for folder in ["sections", "figures"]:
        for path in sorted((ARTICLE / folder).glob("*")):
            if path.is_file():
                arxiv[path.relative_to(ARTICLE).as_posix()] = path.read_bytes()
    arxiv_path = OUT / "H-M_Anchor_Decomposed_r28_arxiv.zip"
    write_zip(arxiv_path, arxiv)

    source_path = OUT / "bimetric_r28_reviewed_source_only.zip"
    source = source_files()
    write_zip(source_path, source)

    pdf_out = OUT / "H-M_Anchor_Decomposed_r28_reviewed.pdf"
    pdf_out.write_bytes(PDF.read_bytes())
    qa = {
        "status": "RESEARCH_PREVIEW_BUILD_COMPLETE",
        "scientific_pass": False,
        "article_pages": pages,
        "figure_assertions": "117/117",
        "source_files_including_manifest": len(source),
        "arxiv_files": len(arxiv),
        "external_publication_is_separate": True,
    }
    qa_path = OUT / "QA_SUMMARY.json"
    qa_path.write_text(json.dumps(qa, indent=2) + "\n")
    readme_out = OUT / "PRERELEASE_README.md"
    readme_out.write_bytes((ROOT / "PRERELEASE_README.md").read_bytes())

    deliver = [pdf_out, arxiv_path, source_path, qa_path, readme_out]
    sums = "".join(f"{sha_file(path)}  {path.name}\n" for path in sorted(deliver))
    sums_path = OUT / "SHA256SUMS.txt"
    sums_path.write_text(sums)
    deliver.append(sums_path)
    full_path = OUT / "H-M_Anchor_Decomposed_r28_research_preview.zip"
    with zipfile.ZipFile(full_path, "w", zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(deliver):
            archive.write(path, path.name)
    print(json.dumps({"status": qa["status"], "pages": pages,
                      "source_files": len(source), "arxiv_files": len(arxiv),
                      "release_files": [p.name for p in deliver],
                      "full_zip_sha256": sha_file(full_path)}, indent=2))


if __name__ == "__main__":
    main()
