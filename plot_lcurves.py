import sys

import matplotlib.pyplot as plt
import pandas as pd


df = pd.read_fwf(sys.argv[1])

plt.figure(figsize=(8, 4))
plt.plot(df["batch"], df["l2_e_trn"], marker=".", label="train loss")
plt.plot(df["batch"], df["l2_e_tst"], marker=".", label="test loss")
plt.yscale("log")
plt.xlabel("Batch")
plt.ylabel("Energy loss")
plt.legend()
plt.title("Energy Loss Curves")
plt.tight_layout()
plt.savefig("lcurve.svg", bbox_inches="tight")
plt.show()
