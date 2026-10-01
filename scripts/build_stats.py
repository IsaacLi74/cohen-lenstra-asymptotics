import csv
import math
from pathlib import Path

from src.dataset import iter_all_fields


OUTPUT = Path("data/processed/p3_counts.csv")

X_MAX = 2**28
N_CHECKPOINTS = 200

TARGETS = {
    "trivial": (),
    "C3": (3,),
    "C9": (9,),
    "C3xC3": (3, 3),
}


# Cohen-Lenstra normalization constant for p = 3
eta3 = 1.0
for j in range(1, 100):
    eta3 *= 1 - 3 ** (-j)

CL = {
    "trivial": eta3,
    "C3": eta3 / 2,
    "C9": eta3 / 6,
    "C3xC3": eta3 / 48,
}


# Logarithmically spaced X-values
log_min = math.log10(1_000)
log_max = math.log10(X_MAX)

checkpoints = sorted(set(
    int(10 ** (
        log_min
        + i * (log_max - log_min) / N_CHECKPOINTS
    ))
    for i in range(N_CHECKPOINTS + 1)
))

if checkpoints[-1] != X_MAX:
    checkpoints.append(X_MAX)


counts = {name: 0 for name in TARGETS}
total = 0
checkpoint_index = 0
rows = []


def save_row(X):
    row = {
        "X": X,
        "N_fields": total,
    }

    for name in TARGETS:
        n = counts[name]
        p_cl = CL[name]

        row[f"N_{name}"] = n
        row[f"P_{name}"] = n / total if total else float("nan")

        # Xu-Zhu style:
        # N_G(X) - P_CL(G) * (3/pi^2) * X
        row[f"E_xuzhu_{name}"] = (
            n - p_cl * (3 / math.pi**2) * X
        )

        # Pure class-group-distribution error:
        # N_G(X) - P_CL(G) * actual number of fields
        row[f"E_conditional_{name}"] = (
            n - p_cl * total
        )

    rows.append(row)


for record in iter_all_fields():
    x = abs(record.D_K)

    while (
        checkpoint_index < len(checkpoints)
        and checkpoints[checkpoint_index] < x
    ):
        save_row(checkpoints[checkpoint_index])
        checkpoint_index += 1

    total += 1

    for name, group in TARGETS.items():
        if record.Cl3_K_invariants == group:
            counts[name] += 1

    while (
        checkpoint_index < len(checkpoints)
        and checkpoints[checkpoint_index] == x
    ):
        save_row(checkpoints[checkpoint_index])
        checkpoint_index += 1

    if total % 5_000_000 == 0:
        print(f"processed {total:,} fields...")


while checkpoint_index < len(checkpoints):
    save_row(checkpoints[checkpoint_index])
    checkpoint_index += 1


OUTPUT.parent.mkdir(parents=True, exist_ok=True)

with OUTPUT.open("w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=rows[0].keys(),
    )
    writer.writeheader()
    writer.writerows(rows)


print()
print(f"fields processed: {total:,}")
print(f"checkpoints:      {len(rows)}")
print(f"saved:            {OUTPUT}")

print()
print("Final checkpoint:")
print(rows[-1])
