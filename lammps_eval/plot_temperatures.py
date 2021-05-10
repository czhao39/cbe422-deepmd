import matplotlib.pyplot as plt
import pandas as pd


plt.figure(figsize=(8, 4))

lj_df = pd.read_table("thermo.ar.metal", delim_whitespace=True)
plt.plot(lj_df["Step"], lj_df["Temp"], marker=".", label="exact")

deepmd_df = pd.read_table("thermo.deepmd.ar.metal", delim_whitespace=True)
plt.plot(deepmd_df["Step"], deepmd_df["Temp"], marker=".", label="deepmd")

plt.legend()
plt.xlabel("Step")
plt.ylabel("Temperature (K)")
plt.title("LAMMPS Molecular Dynamics Simulations")
plt.tight_layout()
plt.savefig("lammps.svg", bbox_inches="tight")
plt.show()
