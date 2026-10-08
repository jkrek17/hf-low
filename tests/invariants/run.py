#!/usr/bin/env python3
"""Project invariants: corrections that were turned into checks.

    python3 tests/invariants/run.py            # run every check
    python3 tests/invariants/run.py --list     # print every violation, known or new

Each check below exists because of an entry in failures/. A check reports
violations as short keys (usually a file path). A violation is one of:

  NEW     not in known_violations.json           -> the run fails
  KNOWN   listed in known_violations.json        -> reported, does not fail
  FIXED   listed, the file still exists, and it  -> the run fails until the line
          no longer violates                        is removed from the list

known_violations.json is a ratchet, not an excuse: it records what was already
wrong when a check was introduced so the check can start catching new cases
straight away. The list should only get shorter. A listed file that does not
exist on the branch being checked is ignored, because main and the research
branch do not hold the same files.

Standard library only. Runs from any directory; reads the working tree and asks
git which files are tracked.
"""

import ast
import json
import os
import re
import subprocess
import sys

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
HERE = os.path.dirname(os.path.abspath(__file__))


def tracked():
    out = subprocess.check_output(["git", "-C", REPO, "ls-files", "-z"])
    return [p for p in out.decode("utf-8", "replace").split("\0") if p]


def read(path):
    with open(os.path.join(REPO, path), encoding="utf-8", errors="replace") as handle:
        return handle.read()


def exists(path):
    return os.path.exists(os.path.join(REPO, path))


# --- stale-claims ------------------------------------------------------------
# failures/2026-10-07-stale-readme-no-gust-index: a README went on saying
# something did not exist after it had been committed. claims.json pairs a
# sentence with the path whose existence makes it false.
def check_stale_claims(files):
    claims = json.load(open(os.path.join(HERE, "claims.json"), encoding="utf-8"))
    violations = []
    for claim in claims:
        if not exists(claim["false_when_exists"]):
            continue
        for path in claim["files"]:
            if exists(path) and claim["sentence"].lower() in " ".join(read(path).lower().split()):
                violations.append("%s :: %s" % (path, claim["sentence"]))
    return violations


# --- seasons-complete ----------------------------------------------------------
# failures/2026-06-24-missing-atlantic-seasons: an analysis ran with four
# Atlantic seasons silently absent. Every season from the record start to the
# latest season in the payload must hold at least one low in each basin.
def check_seasons_complete(files):
    path = "docs/data/hf-lows.json"
    if not exists(path):
        return []
    payload = json.loads(read(path))
    fields = payload["lowFields"]
    basin_at, season_at = fields.index("basin"), fields.index("season")
    seen = {}
    for low in payload["lows"]:
        seen.setdefault(low[basin_at], set()).add(low[season_at])
    start = payload["recordStart"]
    last = max(max(seasons) for seasons in seen.values())
    violations = []
    for basin, seasons in sorted(seen.items()):
        for season in range(start, last + 1):
            if season not in seasons:
                violations.append("%s season %d-%02d has no lows" % (basin, season, (season + 1) % 100))
    return violations


# --- record-start --------------------------------------------------------------
# failures/2026-10-08-validation-floor-stated-two-ways: the first season the
# proxy can be checked against was written as 2001 in some files and 2004 in
# others. The build constant and the payload must agree on 2004, because every
# "fit and test from 2004-05" rule leans on it.
def check_record_start(files):
    violations = []
    build = "tools/build_hf_lows.py"
    if exists(build):
        value = None
        for node in ast.parse(read(build)).body:
            if (isinstance(node, ast.Assign) and len(node.targets) == 1
                    and getattr(node.targets[0], "id", None) == "RECORD_START"):
                value = ast.literal_eval(node.value)
        if value != 2004:
            violations.append("%s RECORD_START is %r, expected 2004" % (build, value))
    payload = "docs/data/hf-lows.json"
    if exists(payload) and json.loads(read(payload)).get("recordStart") != 2004:
        violations.append("%s recordStart is not 2004" % payload)
    return violations


# --- no-session-paths ----------------------------------------------------------
# failures/2026-10-08-hardcoded-home-paths and
# failures/2026-10-08-verifier-scripts-hardcode-clone-path: a script that opens
# a path under one session's home directory, or under /tmp/, runs nowhere else.
# This catches the path. It cannot tell whether a script reads committed files
# only; that is a rule in hf-result-closeout and a question for hf-reviewer.
HOME_PATH = re.compile(r"""["'](/home/[A-Za-z0-9_.-]+/|/Users/[A-Za-z0-9_.-]+/|/root/|/tmp/)""")
CODE = (".py", ".js", ".sh", ".gs")


def check_no_session_paths(files):
    violations = []
    for path in files:
        if not path.endswith(CODE) or path.startswith("tests/invariants/"):
            continue
        if HOME_PATH.search(read(path)):
            violations.append(path)
    return violations


# --- proxy-label ---------------------------------------------------------------
# The ERA5 record is a proxy and must be called one wherever it is described.
# Every README under research/era5/ that mentions ERA5 has to use the word.
def check_proxy_label(files):
    violations = []
    for path in files:
        if not (path.startswith("research/era5/") and path.endswith("README.md")):
            continue
        text = read(path)
        if "ERA5" in text and "proxy" not in text.lower():
            violations.append(path)
    return violations


# --- forbidden-files -----------------------------------------------------------
# failures/2026-10-07-large-files-nearly-committed: the pre-commit hook guards
# a clone that has it enabled; this guards the history whatever the clone did.
LIMIT_BYTES = 10 * 1024 * 1024
FORBIDDEN = (
    re.compile(r"\.(whl|egg|tar\.gz|zip)$"),
    re.compile(r"(^|/)hsf_lows\.csv$"),
    re.compile(r"^data/hsf_cache/"),
    re.compile(r"(^|/)publish\.local\.json$"),
    re.compile(r"(^|/)\.clasp\.json$"),
    re.compile(r"(^|/)\.clasprc\.json$"),
    re.compile(r"^research/era5/event_fields_[0-9]+\.csv$"),
)


def check_forbidden_files(files):
    violations = []
    for path in files:
        if any(pattern.search(path) for pattern in FORBIDDEN):
            violations.append("%s (must never be committed)" % path)
            continue
        full = os.path.join(REPO, path)
        if os.path.isfile(full) and os.path.getsize(full) > LIMIT_BYTES:
            violations.append("%s (%.1f MB, limit 10 MB)" % (path, os.path.getsize(full) / 1048576))
    return violations


# --- failure-library -----------------------------------------------------------
# The library has rules of its own (failures/README.md); run its checker here so
# one command covers everything.
def check_failure_library(files):
    script = os.path.join(REPO, "tools", "failures.py")
    if not os.path.exists(script):
        return []
    result = subprocess.run([sys.executable, script, "check"], capture_output=True, text=True)
    return [line[6:] for line in result.stdout.splitlines() if line.startswith("FAIL  ")]


CHECKS = [
    ("stale-claims", check_stale_claims),
    ("seasons-complete", check_seasons_complete),
    ("record-start", check_record_start),
    ("no-session-paths", check_no_session_paths),
    ("proxy-label", check_proxy_label),
    ("forbidden-files", check_forbidden_files),
    ("failure-library", check_failure_library),
]


def key_path(key):
    """The file a violation key refers to: the text before the first space."""
    return key.split(" ", 1)[0]


def main(argv):
    show_all = "--list" in argv
    known = json.load(open(os.path.join(HERE, "known_violations.json"), encoding="utf-8"))
    files = tracked()
    failed = False
    totals = {"new": 0, "known": 0, "fixed": 0}
    for name, check in CHECKS:
        found = check(files)
        listed = {item["key"]: item for item in known.get(name, [])}
        new = [v for v in found if v not in listed]
        still = [v for v in found if v in listed]
        fixed = [k for k in listed if k not in found and exists(key_path(k))]
        totals["new"] += len(new)
        totals["known"] += len(still)
        totals["fixed"] += len(fixed)
        status = "FAIL" if (new or fixed) else "ok  "
        print("%s %-17s %d new, %d known, %d fixed-but-still-listed"
              % (status, name, len(new), len(still), len(fixed)))
        for violation in new:
            print("       NEW    " + violation)
        for key in fixed:
            print("       FIXED  %s  (remove it from known_violations.json)" % key)
        if show_all:
            for violation in still:
                print("       KNOWN  %s  [%s]" % (violation, listed[violation].get("failure", "")))
        failed = failed or bool(new or fixed)
    print()
    print("invariants: %(new)d new, %(known)d known, %(fixed)d fixed-but-still-listed" % totals)
    if failed:
        print("A NEW violation is a correction waiting to happen: fix it, or if it is a "
              "false alarm, fix the check. Do not add it to known_violations.json.")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
