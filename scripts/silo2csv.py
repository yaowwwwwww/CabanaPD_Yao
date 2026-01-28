# silo2csv_paraview.py
from paraview.simple import *
import glob, os

# 你要导出的点变量
arrays = [
    "rank_0/displacements",
    "rank_0/velocities",
    "rank_0/forces",
    "rank_0/damage",
    "rank_0/strain_energy",
    "rank_0/type",
]

files = sorted(glob.glob("particles_*.silo"))

for f in files:
    print("Exporting", f)
    base = os.path.splitext(os.path.basename(f))[0]

    reader = OpenDataFile(f)

    # 勾选所有需要的数组
    reader.PointArrayStatus = arrays
    reader.UpdatePipeline()

    # 合并成单一 block
    merged = MergeBlocks(Input=reader)
    merged.UpdatePipeline()

    # 保存为 CSV
    outcsv = base + "_all.csv"
    SaveData(outcsv, proxy=merged, PointDataArrays=arrays)
    print(" -> wrote", outcsv)

