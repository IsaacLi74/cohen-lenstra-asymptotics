import argparse
from itertools import islice

from src.dataset import iter_all_fields


parser = argparse.ArgumentParser()
parser.add_argument("--n", type=int, default=20)
args = parser.parse_args()


print(
    f"{'d':>10} "
    f"{'D_K':>10} "
    f"{'h_K':>6} "
    f"{'Cl_K':>20} "
    f"{'Cl_K[3∞]':>20}"
)

print("-" * 72)

for record in islice(iter_all_fields(), args.n):
    print(
        f"{record.d:10d} "
        f"{record.D_K:10d} "
        f"{record.h_K:6d} "
        f"{record.Cl_K:>20} "
        f"{record.Cl3_K:>20}"
    )
