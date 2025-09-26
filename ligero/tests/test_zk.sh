# 第一个命令行参数
cmd=$1
input_size=$2

#!/bin/bash
# 使用方法
## 基于安全多方计算的零知识证明模块-总体功能测试
# ./test_zk.sh 1
## 基于安全多方计算的零知识证明模块-输入规模扩展测试
# ./test_zk.sh 2 [1,2,4,8] 
## 基于安全多方计算的零知识证明模块-可插拔及重复执行测试
# ./test_zk.sh 3 ""  
## 基于安全多方计算的零知识证明模块-主要性能指标测试
# ./test_zk.sh 4 [1,2,4,8] 
 
case $cmd in
  1)
    echo "Running test_zk-1"
    cd ../../build/ligero
    ./test_zk-1
    ;;
  2)
    echo "Running test_zk-2 with input_size = $input_size"
    cd ../../build/ligero
    ./test_zk-2 $input_size
    ;;
  3)
    echo "Running test_zk-3"
    cd ../../build/ligero
    ./test_zk-3
    ;;
  4)
    echo "Running test_zk-4 with input_size = $input_size"
    cd ../../build/ligero
    ./test_zk-4 $input_size
    ;;
esac
