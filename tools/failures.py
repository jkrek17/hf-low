#!/usr/bin/env python3
"""Check and summarise the failure library in failures/.

    python3 tools/failures.py check    # format, and that every repeat changed the workflow
    python3 tools/failures.py tally    # counts by week, class, who caught it, and what is open

The library's rules are in failures/README.md. This script is standard library
only, like the rest of tools/, so it runs in CI and on the production box alike.

Why this exists: a list of lessons that nobody counts cannot show whether the
project is repeating itself. `tally` is the number to watch: entries per week,
how many were caught by Jason rather than by a check or a reviewer, and how many
repeat a class that was already in the library.
"""

import datetime
import os
import sys
from collections import Counter, defaultdict

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FAILURES_DIR = os.path.join(REPO, "failures")

REQUIRED = ("id", "date", "class", "caught_by", "became", "became_ref")
OPTIONAL = ("repeat_of", "workflow_change")
CAUGHT_BY = {"jason", "verifier", "reviewer", "ci", "thread"}
BECAME = {"test", "hook", "skill", "rule", "gate", "none"}


def parse(path):
    """Return (fields, body) for one entry. The header is `key: value` lines
    between two `---` markers; values are taken literally to the end of line."""
    with open(path, encoding="utf-8") as handle:
        lines = handle.read().split("\n")
    if not lines or lines[0].strip() != "---":
        raise ValueError("does not start with a --- header")
    try:
        end = lines.index("---", 1)
    except ValueError:
        raise ValueError("header is not closed with ---")
    fields = {}
    for line in lines[1:end]:
        if not line.strip():
            continue
        if ":" not in line:
            raise ValueError("header line without a colon: %r" % line)
        key, value = line.split(":", 1)
        fields[key.strip()] = value.strip()
    return fields, "\n".join(lines[end + 1:]).strip()


def load():
    """All entries, oldest first, plus a list of problems found while reading."""
    entries, problems = [], []
    if not os.path.isdir(FAILURES_DIR):
        return entries, ["failures/ does not exist"]
    for name in sorted(os.listdir(FAILURES_DIR)):
        if not name.endswith(".md") or name == "README.md":
            continue
        path = os.path.join(FAILURES_DIR, name)
        try:
            fields, body = parse(path)
        except ValueError as error:
            problems.append("%s: %s" % (name, error))
            continue
        fields["_file"] = name
        fields["_body"] = body
        entries.append(fields)
    entries.sort(key=lambda e: (e.get("date", ""), e.get("id", "")))
    return entries, problems


def known_classes():
    """Classes listed in failures/README.md: the first, backticked column of
    its table. A new class is added there in the same pull request that first
    uses it, so a class missing from the table is a typo or a skipped step."""
    path = os.path.join(FAILURES_DIR, "README.md")
    if not os.path.exists(path):
        return set()
    classes = set()
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("| `") and "` |" in line:
                classes.add(line.split("`")[1])
    return classes


def check(entries, problems):
    """Format problems, and repeats that did not change the workflow."""
    classes = known_classes()
    for e in entries:
        if classes and e.get("class") and e["class"] not in classes:
            problems.append("%s: class %r is not in the table in failures/README.md"
                            % (e["_file"], e["class"]))
    by_id = {e.get("id"): e for e in entries}
    for e in entries:
        name = e["_file"]
        for key in REQUIRED:
            if not e.get(key):
                problems.append("%s: missing %s" % (name, key))
        for key in e:
            if not key.startswith("_") and key not in REQUIRED + OPTIONAL:
                problems.append("%s: unknown field %s" % (name, key))
        if e.get("id") and e["id"] + ".md" != name:
            problems.append("%s: id %r does not match the file name" % (name, e["id"]))
        try:
            datetime.date.fromisoformat(e.get("date", ""))
        except ValueError:
            problems.append("%s: date %r is not YYYY-MM-DD" % (name, e.get("date")))
        if e.get("caught_by") and e["caught_by"] not in CAUGHT_BY:
            problems.append("%s: caught_by %r is not one of %s"
                            % (name, e["caught_by"], ", ".join(sorted(CAUGHT_BY))))
        if e.get("became") and e["became"] not in BECAME:
            problems.append("%s: became %r is not one of %s"
                            % (name, e["became"], ", ".join(sorted(BECAME))))
        for heading in ("**What happened.**", "**How it was caught.**",
                        "**What it was turned into.**"):
            if heading not in e["_body"]:
                problems.append("%s: body has no paragraph starting %s" % (name, heading))
    # Rule 2 in the README. Which entry came first in a class is stated by the
    # entries themselves (the one with no repeat_of), not inferred from dates:
    # two failures can share a day, and file order would then decide by accident.
    by_class = defaultdict(list)
    for e in entries:
        if e.get("class"):
            by_class[e["class"]].append(e)
    for cls, group in sorted(by_class.items()):
        firsts = [e for e in group if not e.get("repeat_of")]
        if len(firsts) != 1:
            problems.append("class %s has %d entries without repeat_of (%s); exactly one "
                            "is the first, the others set repeat_of to it"
                            % (cls, len(firsts), ", ".join(e["_file"] for e in firsts) or "none"))
        for e in group:
            if not e.get("repeat_of"):
                continue
            target = by_id.get(e["repeat_of"])
            if target is None:
                problems.append("%s: repeat_of %r is not an entry" % (e["_file"], e["repeat_of"]))
            elif target.get("class") != cls:
                problems.append("%s: repeat_of points at a %s entry, not %s"
                                % (e["_file"], target.get("class"), cls))
            if not e.get("workflow_change"):
                problems.append("%s: second or later %s failure with no workflow_change; "
                                "fixing the instance is not enough" % (e["_file"], cls))
    return problems


def week_of(date_text):
    year, week, _ = datetime.date.fromisoformat(date_text).isocalendar()
    return "%d-W%02d" % (year, week)


def tally(entries):
    by_week = defaultdict(list)
    for e in entries:
        by_week[week_of(e["date"])].append(e)
    print("Failure library: %d entries, %d classes" %
          (len(entries), len({e["class"] for e in entries})))
    print()
    print("%-9s %7s %9s %8s   %s" % ("week", "entries", "by Jason", "repeats", "classes"))
    for week in sorted(by_week):
        rows = by_week[week]
        repeats = sum(1 for e in rows if e.get("repeat_of"))
        print("%-9s %7d %9d %8d   %s" % (
            week, len(rows), sum(1 for e in rows if e["caught_by"] == "jason"), repeats,
            ", ".join("%s x%d" % (c, n) if n > 1 else c
                      for c, n in sorted(Counter(e["class"] for e in rows).items()))))
    print()
    print("By class (two or more means the workflow should already have changed):")
    for cls, n in sorted(Counter(e["class"] for e in entries).items(),
                         key=lambda item: (-item[1], item[0])):
        print("  %2d  %s" % (n, cls))
    print()
    print("What entries were turned into:")
    for form, n in sorted(Counter(e["became"] for e in entries).items(),
                          key=lambda item: (-item[1], item[0])):
        print("  %2d  %s" % (n, form))
    open_entries = [e for e in entries if e["became"] == "none"]
    if open_entries:
        print()
        print("Not yet converted:")
        for e in open_entries:
            print("  %s" % e["id"])
    print()
    print("Read it as: entries per week and 'by Jason' should fall; 'repeats' should stay "
          "near zero. A week with repeats means a workflow change did not hold.")


def main(argv):
    command = argv[1] if len(argv) > 1 else "tally"
    entries, problems = load()
    if command == "check":
        problems = check(entries, problems)
        for problem in problems:
            print("FAIL  " + problem)
        print("failure library: %d entries, %d problem%s"
              % (len(entries), len(problems), "" if len(problems) == 1 else "s"))
        return 1 if problems else 0
    if command == "tally":
        if problems:
            for problem in problems:
                print("unreadable: " + problem)
        tally(entries)
        return 0
    print(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv))
