#!/usr/bin/env python3
"""Plot wMel genomic titer vs day, colored by treatment group."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from io import StringIO

rng = np.random.default_rng(0)

data = """color,day,wmel_titer
#818181,100,46.59034653
#f58020,0,0
#f58020,0,0
#f58020,100,15.11688835
#f58020,100,15.79281647
#f58020,0,0
#f58020,0,0
#f58020,1,0.2272727273
#f58020,1,0.2
#f58020,7,0.3574191552
#f58020,7,0.4682739436
#818181,0,0
#818181,100,13.45035219
#818181,13,2.968533669
#818181,13,0.5367078301
#818181,13,1.413956239
#818181,13,4.105902778
#818181,0,0
#818181,100,25.69849875
#818181,20,10.14374226
#818181,20,4.418145957
#818181,20,11.80342651
#818181,20,16.77629488
#f58020,28,37.44752624
#f58020,28,60.00745712
#818181,0,0
#818181,100,30.67985612
#818181,34,109.4312169
#818181,34,109.1885442
#818181,34,95.32941176
#818181,34,63.72996301
#f58020,0,0
#f58020,100,37.06967742
#818181,42,47.84722222
#818181,42,118.3116883
#818181,42,83.13953488
#818181,42,23.79966887
#f58020,56,7.727070677
#f58020,56,6.642842962
#1779b8,0,0
#1779b8,0,0
#1779b8,100,50.94326726
#1779b8,100,52.01349432
#1779b8,1,6.762354651
#1779b8,1,5.259340659
#1779b8,0,0
#1779b8,0,0
"""

df = pd.read_csv(StringIO(data))

fig, ax = plt.subplots(figsize=(4, 2))

# Plot each color group separately so points of the same group share a legend handle
for color, group in df.groupby("color", sort=False):
    jitter = rng.uniform(-0.6, 0.6, size=len(group))
    ax.scatter(
        group["day"] + jitter,
        group["wmel_titer"],
        c=color,
        edgecolors="black",
        linewidths=0.4,
        s=20,
        alpha=0.5,
    )

ax.set_xlabel("Day")
ax.set_ylabel("wMel genomic titer")
ax.set_yscale("symlog", linthresh=0.1)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

fig.tight_layout()
fig.savefig("/mnt/user-data/outputs/wmel_titer.svg", format="svg")
print("Wrote /mnt/user-data/outputs/wmel_titer.svg")