#!/usr/bin/env python3
"""Verify that environment.yml still covers everything the code uses.

environment.yml is the single source of truth for the genDiv runtime
environment. This script proves that claim rather than trusting it: it
scans the pipeline for the Python modules, R packages, and CLI tools that
are actually invoked, maps each to its conda package name, and reports any
that environment.yml does not declare.

It exists because the failure it catches is silent. A package can be
missing from the spec and the pipeline still runs, because some other
dependency happened to pull it in — that is exactly how r-dplyr, r-tidyr,
and tqdm went undeclared while the pipeline worked fine on the machine
where they had been installed transitively. Nothing breaks until someone
builds the environment from scratch.

Runs on a bare python3 (stdlib only), so it works without the genDiv
environment activated and in CI.

Exit status:
    0  every used dependency is declared
    1  at least one is missing (or --strict and something is unused)
"""

import argparse
import re
import sys
from pathlib import Path

# Import name / command name -> conda package name, where they differ.
PACKAGE_ALIASES = {
    "allel": "scikit-allel",
    "matplotlib": "matplotlib-base",
    "sklearn": "scikit-learn",
    "yaml": "pyyaml",
    "PIL": "pillow",
    "cv2": "opencv",
    # CLI tools shipped by a differently-named package
    "tabix": "htslib",
    "bgzip": "htslib",
    "Rscript": "r-base",
}

# CLI tools the pipeline may invoke that we expect conda to provide.
# Shell builtins and coreutils (awk, sort, paste, …) are deliberately absent:
# they come from the OS, not the environment spec.
KNOWN_CLI_TOOLS = {
    "plink", "plink2", "bcftools", "bedtools", "beagle",
    "tabix", "bgzip", "Rscript",
}

# Supplied outside conda; see the header note in environment.yml.
EXTERNALLY_PROVIDED = {"rclone"}

R_PACKAGE_PREFIX = "r-"


def parse_environment_yml(path):
    """Return declared package names, lowercased, without version pins.

    Deliberately regex-based rather than pyyaml: this script must run on a
    bare interpreter, and environment.yml is a flat list we control.
    """
    names = set()
    in_deps = False
    for raw in path.read_text().splitlines():
        line = raw.rstrip()
        if re.match(r"^dependencies:\s*$", line):
            in_deps = True
            continue
        if in_deps and re.match(r"^[a-zA-Z]", line):
            break  # a new top-level key ends the dependency block
        if not in_deps:
            continue
        m = re.match(r"^\s*-\s*([A-Za-z0-9_.-]+)", line)
        if m:
            names.add(m.group(1).lower())
    return names


def python_imports(py_files):
    """Third-party top-level imports, excluding stdlib and local modules."""
    stdlib = set(getattr(sys, "stdlib_module_names", ()))
    local = {p.stem for p in py_files}
    found = {}
    pattern = re.compile(r"^\s*(?:import|from)\s+([A-Za-z0-9_]+)", re.M)
    for path in py_files:
        for mod in pattern.findall(path.read_text()):
            if mod in stdlib or mod in local or mod.startswith("_"):
                continue
            found.setdefault(mod, set()).add(path)
    return found


def r_packages(r_files):
    found = {}
    pattern = re.compile(r"(?:library|require)\(\s*([A-Za-z0-9._]+)")
    for path in r_files:
        for pkg in pattern.findall(path.read_text()):
            found.setdefault(pkg, set()).add(path)
    return found


def cli_tools(sh_files):
    """CLI tools invoked at the start of a line or after a pipe."""
    found = {}
    pattern = re.compile(
        r"(?:^|\||&&|\$\()\s*([A-Za-z][A-Za-z0-9_]*)\b", re.M)
    for path in sh_files:
        text = "\n".join(
            l for l in path.read_text().splitlines()
            if not l.lstrip().startswith("#"))
        for tool in pattern.findall(text):
            if tool in KNOWN_CLI_TOOLS:
                found.setdefault(tool, set()).add(path)
    return found


def to_conda_name(name, kind):
    if name in PACKAGE_ALIASES:
        return PACKAGE_ALIASES[name]
    if kind == "R":
        return R_PACKAGE_PREFIX + name.lower()
    return name.lower()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=str(Path(__file__).resolve().parent.parent),
                    help="Repo root (default: parent of this script's directory).")
    ap.add_argument("--env-file", default=None,
                    help="Path to environment.yml (default: <root>/environment.yml).")
    ap.add_argument("--strict", action="store_true",
                    help="Also fail when environment.yml declares a package no "
                         "code appears to use (off by default: a package can be "
                         "needed at runtime without any literal reference).")
    args = ap.parse_args()

    root = Path(args.root).resolve()
    env_file = Path(args.env_file) if args.env_file else root / "environment.yml"
    if not env_file.is_file():
        sys.exit(f"[check_env] no environment.yml at {env_file}")

    declared = parse_environment_yml(env_file)

    py_files = sorted(root.glob("scripts/**/*.py")) + sorted(root.glob("explore/*.py"))
    py_files = [p for p in py_files if p.name != Path(__file__).name]
    r_files = sorted(root.glob("scripts/**/*.R"))
    sh_files = (sorted(root.glob("genDiversity*.sh"))
                + sorted(root.glob("*.sh"))
                + sorted(root.glob("scripts/**/*.sh")))
    sh_files = sorted(set(sh_files))

    used = {}   # conda name -> (kind, original name, files)
    for mod, files in python_imports(py_files).items():
        used[to_conda_name(mod, "py")] = ("python", mod, files)
    for pkg, files in r_packages(r_files).items():
        used[to_conda_name(pkg, "R")] = ("R", pkg, files)
    for tool, files in cli_tools(sh_files).items():
        if tool in EXTERNALLY_PROVIDED:
            continue
        used[to_conda_name(tool, "cli")] = ("CLI", tool, files)

    missing = {k: v for k, v in used.items() if k not in declared}
    unused = sorted(declared - set(used) - {"python", "r-base"})

    print(f"[check_env] environment.yml declares {len(declared)} packages")
    print(f"[check_env] code references {len(used)} packages across "
          f"{len(py_files)} .py, {len(r_files)} .R, {len(sh_files)} .sh files")

    if missing:
        print("\nMISSING from environment.yml — used but not declared:")
        for conda_name in sorted(missing):
            kind, orig, files = missing[conda_name]
            where = ", ".join(sorted(str(f.relative_to(root)) for f in files)[:3])
            label = f"{orig} -> {conda_name}" if orig != conda_name else conda_name
            print(f"  [{kind:6}] {label:28} {where}")

    if unused:
        print("\nDeclared but no literal use found (may still be needed "
              "transitively or at runtime):")
        for name in unused:
            print(f"  {name}")

    if missing:
        noun = "dependency" if len(missing) == 1 else "dependencies"
        print(f"\n[check_env] FAIL: {len(missing)} undeclared {noun}")
        return 1
    if unused and args.strict:
        print(f"\n[check_env] FAIL (--strict): {len(unused)} unused declaration(s)")
        return 1
    print("\n[check_env] OK: environment.yml covers every dependency the code uses")
    return 0


if __name__ == "__main__":
    sys.exit(main())
