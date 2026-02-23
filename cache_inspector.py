import pandas as pd
from itertools import islice

path = r"Inputs\rendered_scenarios\PlanetMaker_broken.csv"
target_lines = set(range(1, 17000)) | {9862203, 9862204, 9862205, 9862206}
# ---- read only the needed lines (plus line 1) ----
wanted = sorted(target_lines)
raw_lines = {}  # line_no -> decoded line string

with open(path, "rb") as f:
    want_idx = 0
    current_target = wanted[want_idx] if wanted else None

    for i, bline in enumerate(f, start=1):
        if current_target is None:
            break
        if i == current_target:
            raw_lines[i] = bline.decode("utf-8", errors="replace").rstrip("\r\n")
            want_idx += 1
            current_target = wanted[want_idx] if want_idx < len(wanted) else None

# ---- split into fields and pad to max width so different lengths work ----
rows = []
max_fields = 0
for ln in wanted:
    s = raw_lines.get(ln, "")
    fields = s.split(",") if s != "" else []
    max_fields = max(max_fields, len(fields))
    rows.append({"line_no": ln, "raw": s, "fields": fields})

# pad fields to same length
for r in rows:
    r["fields"] = r["fields"] + [""] * (max_fields - len(r["fields"]))
    r["n_fields"] = len(r["raw"].split(",")) if r["raw"] != "" else 0

# build dataframe: one row per selected line, columns col_1..col_N
df = pd.DataFrame({
    "line_no": [r["line_no"] for r in rows],
    "n_fields": [r["n_fields"] for r in rows],
    "raw": [r["raw"] for r in rows],
})

field_df = pd.DataFrame([r["fields"] for r in rows],
                        columns=[f"col_{j+1}" for j in range(max_fields)])

df = pd.concat([df, field_df], axis=1)
df_sub = df[df["col_2"] == 12]
df_sub2 = df[df["col_2"].astype(str) == "3211"]
print(df)
print("stop")