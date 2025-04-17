# 编译

运行

```bash
sudo apt-get install build-essential cmake git libgmp3-dev libprocps4-dev libboost-all-dev libssl-dev libsodium-dev --fix-missing
git submodule init && git submodule update
```

然后编译
```bash
mkdir build
cd build
cmake ..
make
```

# 批量化零知识内积论证模块

## （1）	批量化零知识内积论证模块-总体功能测试

运行
```bash
cd build
cmake ..
make
./ligero/test_PCS
```

可以看到控制台打印

```bash
这是批量化内积论证总体功能测试 !
批量化内积论证运行成功!
```
字样

这代表批量化零知识内积论证模块-总体功能测试通过

## （2）	批量化零知识内积论证模块-输入规模扩展测试

修改 [ligero/test_PCS.cpp](ligero/tests/test_PCS.cpp) 中 variable_num 的值以调整输入向量规模。
可以调整为8,9,10,11,12.

对于每一个 variable_num，运行
```bash
cd build
cmake ..
make
./ligero/test_PCS
```

可以看到控制台打印

```bash
这是批量化内积论证总体功能测试 !
批量化内积论证运行成功!
向量长度: xx
```
字样

这代表批量化零知识内积论证模块-输入规模扩展测试通过

## （3）	批量化零知识内积论证模块-异常处理测试

修改 [ligero/test_PCS.cpp](ligero/tests/test_PCS.cpp) 中 variable_num = 10。

将参数 test_type 设置为 0。

运行
```bash
cd build
cmake ..
make
./ligero/test_PCS
```

可以看到控制台打印

```bash
这是批量化内积论证异常处理测试!
批量化内积论证运行失败!
批量化内积论证异常处理测试运行成功!
```
字样

这代表批量化零知识内积论证模块-异常处理测试通过

## （4）	批量化零知识内积论证模块-批量处理测试

修改 [ligero/test_PCS.cpp](ligero/tests/test_PCS.cpp) 中 test_type = 1。

将 instance 设置为 8, 9 ,10, 11, 12.

对于每一个 instance 值，运行
```bash
cd build
cmake ..
make
./ligero/test_PCS
```

可以看到控制台打印

```bash
这是批量化内积论证总体功能测试!
批量化内积论证运行成功!
批处理规模: xx
```
字样

这代表批量化零知识内积论证模块-批量处理测试通过

## （5）	线性约束检查单元-总体功能测试

修改 [ligero/test_PCS.cpp](ligero/tests/test_PCS.cpp) 中 instance = 1。

运行
```bash
cd build
cmake ..
make
./ligero/test_PCS
```

可以看到控制台打印

```bash
线性约束检查模块通过！
线性约束检查通信量：xxx
```
字样

这代表线性约束检查单元-总体功能测试通过

# 基于安全多方计算的零知识证明模块

## （1）	基于安全多方计算的零知识证明模块-总体功能测试


运行
```bash
cd build
cmake ..
make
./ligero/test_debugg
```

可以看到控制台打印

```bash
这是基于安全多方计算的零知识证明在正常情况下的测试！
基于安全多方计算的零知识证明运行成功！
基于安全多方计算的零知识证明总体功能测试通过!
```
字样

这代表基于安全多方计算的零知识证明模块-总体功能测试通过

## （2）	基于安全多方计算的零知识证明模块-输入规模扩展测试

修改 [ligero/test_debugg.cpp](ligero/tests/test_debugg.cpp) 中 input_size = 1, 2, 4或 8。

对于每一个input_size，运行
```bash
cd build
cmake ..
make
./ligero/test_debugg
```

可以看到控制台打印

```bash
这是基于安全多方计算的零知识证明在正常情况下的测试！
输入规模为深度为 xx 的默克尔树!
基于安全多方计算的零知识证明运行成功！
```
字样

这代表基于安全多方计算的零知识证明模块-输入规模扩展测试通过

## （3）	基于安全多方计算的零知识证明模块-异常处理测试

修改 [ligero/test_debugg.cpp](ligero/tests/test_debugg.cpp) 中 input_size = 1。

修改 test_type = 0。

运行
```bash
cd build
cmake ..
make
./ligero/test_debugg
```

可以看到控制台打印

```bash
这是基于安全多方计算的零知识证明异常处理测试！
基于安全多方计算的零知识证明运行失败！
基于安全多方计算的零知识证明异常处理测试通过！
```
字样

这代表基于安全多方计算的零知识证明模块-异常处理测试通过

## （4）	基于安全多方计算的零知识证明模块-可插拔及重复执行测试

修改 [ligero/test_debugg.cpp](ligero/tests/test_debugg.cpp) 中  test_type = 1。

运行
```bash
cd build
cmake ..
make
./ligero/test_debugg
```

可以看到控制台打印

```bash
基于安全多方计算的零知识证明Ligerolight模块（模块1）执行通过！
```
字样

这代表基于安全多方计算的零知识证明模块-模块1可插拔测试通过

多次运行
```bash
cd build
cmake ..
make
./ligero/test_debugg
```

可以看到控制台打印

```bash
基于安全多方计算的零知识证明Ligerolight模块（模块1）执行通过！
```
字样

这代表基于安全多方计算的零知识证明模块-模块1重复执行测试通过
