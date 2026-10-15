import argparse

from src.dataset import group_label, iter_all_fields


def main():
    parser = argparse.ArgumentParser(
        description="List imaginary quadratic fields by discriminant bound."
    )
    parser.add_argument(
        "--max-D",
        type=int,
        required=True,
        help="List all fields with |D_K| <= this bound.",
    )
    args = parser.parse_args()

    bound = args.max_D

    print(
        f"{'d':>6} "
        f"{'D_K':>6} "
        f"{'h_K':>5} "
        f"{'Cl_K':>15} "
        f"{'Sylow3':>15} "
        f"{'|Cl[3]|':>8} "
        f"{'3|h':>4}"
    )
    print("-" * 68)

    count = 0

    for record in iter_all_fields():
        if abs(record.D_K) > bound:
            break

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
