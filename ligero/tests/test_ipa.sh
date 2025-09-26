# 第一个命令行参数
cmd=$1
var_num=$2
instance_num=$3

#!/bin/bash
# 使用方法
## 批量化零知识内积论证模块-总体功能测试
# ./test_ipa.sh 1
## 批量化零知识内积论证模块-输入规模扩展测试
# ./test_ipa.sh 2 [8,9,10,11,12]
## 批量化零知识内积论证模块-批量处理测试
# ./test_ipa.sh 3 "" [8,9,10,11,12]
## 批量化零知识内积论证模块-主要性能指标测试
# ./test_ipa.sh 4

case $cmd in
  1)
    echo 
    cd ../../build/ligero
    ./test_ipa-1
    ;;
  2)
    echo "Running test_ipa-2 with var_num = $var_num"
    cd ../../build/ligero
    ./test_ipa-2 $var_num
    ;;
  3)
    echo "Running test_ipa-3 with instance_num = $instance_num"
    cd ../../build/ligero
    ./test_ipa-3 $instance_num
    ;;
  4)
    echo "Running test_ipa-4"
    cd ../../build/ligero
    ./test_ipa-4
    ;;
esac
