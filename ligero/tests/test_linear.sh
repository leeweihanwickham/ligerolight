# 第一个命令行参数
cmd=$1

#!/bin/bash
# 使用方法
## 线性约束检查单元-总体功能测试
# ./test_linear.sh 1
## 线性约束检查单元-通信量优化对比测试
## 優化前
# ./test_linear.sh 2
## 优化后
# ./test_linear.sh 3

case $cmd in
  1)
    echo 
    cd ../../build/ligero
    ./test_linear-1
    ;;
  2)
    echo "Running test_linear-2-1"
    cd ../../build/ligero
    ./test_linear-2-1
    ;;
  3)
    echo "Running test_linear-2-2"
    cd ../../build/ligero
    ./test_linear-2-2
    ;;
esac
