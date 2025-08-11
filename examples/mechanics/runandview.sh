#!/bin/bash



# 运行模拟程序
cd /home/wuwen/program/CabanaPD_Yao/output/
echo "Running simulation..."
/home/wuwen/program/CabanaPD_Yao/build/examples/mechanics/ColdSprayImpact \
  /home/wuwen/program/CabanaPD_Yao/examples/mechanics/inputs/simple_impact.json

# 等待关键结果文件
if [ $? -eq 0 ]; then
  echo "Simulation completed."
else
  echo "Waring: Simulation failed with error code $?"
  # 你可以选择在这里 exit 1 来阻止 ParaView 启动，如果模拟失败了。
  # exit 1
fi
echo ""

# 直接调用你自己写好的 Python 脚本可视化
/home/wuwen/Downloads/ParaView-6.0.0-RC1-MPI-Linux-Python3.12-x86_64/bin/paraview \
  --script=/home/wuwen/program/CabanaPD_Yao/examples/mechanics/rload_and_play.py
