import os
import pandas as pd

path = r"Inputs\rendered_scenarios\PlanetMaker.csv"



root, ext = os.path.splitext(path)
out_path = root + "_fix" + (ext if ext else ".csv")

# line_no -> {"actual": int, "missing": int, "extra": int}
mismatch_info = {}

chunk_rows = 100_000
buffer_rows = []
written_any = False

def flush_buffer():
    global written_any, buffer_rows
    if not buffer_rows:
        return
    df = pd.DataFrame(buffer_rows)  # store as df before saving
    df.to_csv(out_path, mode="a", header=False, index=False)
    written_any = True
    buffer_rows = []

with open(path, "rb") as fin:
    # Read + write first line exactly as-is
    first_b = fin.readline()
    first_s = first_b.decode("utf-8", errors="replace").rstrip("\r\n")
    first_stripped = first_b.strip()
    expected = (first_stripped.count(b",") + 1) if first_stripped else 0

    # overwrite output and write first line verbatim
    with open(out_path, "w", encoding="utf-8", newline="") as fout:
        fout.write(first_s + "\n")

    total_lines = 1

    # Process remaining lines
    for line_no, bline in enumerate(fin, start=2):
        total_lines = line_no
        stripped = bline.strip()

        if not stripped:
            actual = 0
            values = []
        else:
            actual = stripped.count(b",") + 1
            line_s = bline.decode("utf-8", errors="replace").rstrip("\r\n")
            values = line_s.split(",")

        # record mismatch info (missing relative to expected)
        if actual != expected:
            missing = max(0, expected - actual)
            extra   = max(0, actual - expected)
            mismatch_info[line_no] = {"actual": actual, "missing": missing, "extra": extra}

        # apply fixes + buffer rows
        if expected == 0:
            # nothing sensible to do; just write empty rows as-is
            buffer_rows.append(values)
        elif actual == expected:
            buffer_rows.append(values)
        elif actual < expected:
            buffer_rows.append(values + [""] * (expected - actual))
        else:
            # actual > expected
            if actual == 2 * expected:
                buffer_rows.append(values[:expected])
                buffer_rows.append(values[expected:2 * expected])
            else:
                # not double: keep first expected fields (still consistent width)
                print(f"not double! {expected} expected, but found {actual}")
                buffer_rows.append(values[:expected])


        if len(buffer_rows) >= chunk_rows:
            flush_buffer()

flush_buffer()

print(f"EXPECTED FIELD COUNT (line 1): {expected}")
print(f"TOTAL LINES (input): {total_lines}")
print(f"MISMATCHING LINES: {len(mismatch_info)}")
print(f"WRITTEN FIXED CSV: {out_path}")

# If you really want the full list (can be huge):
print(list(mismatch_info.keys()))

print("stop")