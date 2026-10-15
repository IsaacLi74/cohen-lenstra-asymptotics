import argparse

from src.dataset import group_label, iter_all_fields
from src.statistics import DATA_BOUND


def main():
    parser = argparse.ArgumentParser(
        description="List class numbers and class-group decompositions from the committed data."
    )
    bounds = parser.add_mutually_exclusive_group(required=True)
    bounds.add_argument(
        "--max-D",
        type=int,
        help="List all fields with |D_K| <= this bound.",
    )
    bounds.add_argument(
        "--max-d",
        type=int,
        help="List all fields Q(sqrt(d)) with squarefree d < 0 and |d| <= this bound.",
    )
    args = parser.parse_args()

    radicand = args.max_d is not None
    bound = args.max_d if radicand else args.max_D
    maximum = DATA_BOUND//4 if radicand else DATA_BOUND
    if not 1 <= bound <= maximum:
        parser.error(f"Bound must be between 1 and {maximum} for complete coverage by the committed data.")
    discriminant_bound = 4*bound if radicand else bound

    print(
        f"{'d':>6} "
        f"{'D':>6} "
        f"{'h_D':>5} "
        f"{'Cl_D':>15} "
        f"{'Sylow3':>15} "
        f"{'|Cl[3]|':>8} "
        f"{'3|h':>4}"
    )
    print("-" * 68)

    count = 0

    for record in iter_all_fields():
        if abs(record.D_K) > discriminant_bound:
            break
        if radicand and abs(record.d) > bound:
            continue

        print(
            f"{record.d:6d} "
            f"{record.D_K:6d} "
            f"{record.h_K:5d} "
            f"{record.Cl_K:>15} "
            f"{group_label(record.sylow3_invariants):>15} "
            f"{record.torsion3_size:8d} "
            f"{int(record.h_K % 3 == 0):4d}"
        )
        count += 1

    print("-" * 68)
    print(f"total fields: {count}")


if __name__ == "__main__":
    main()
