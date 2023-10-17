> 学习并测试多晶-晶体塑性模拟的笔记及备忘
# TODO
1. 完善并理解 `CrystalPlasticityKalidindi` + `ComputeMultipleCrystalPlasticityStress` 的计算逻辑框架
2. 调研在晶体塑性中加入 `back stree` 的目的以及背后的物理机制是什么？
   1. 是否back stress 是一个特征参量，与GND相关，需要考虑各种微结构的硬化效应，参考北京大学张寅的文章；
   2. 系统化测试晶体塑性部分的代码

# 构建初始模型 neper
> Here's a step-by-step guide for using Neper-generated polycrystalline structures as initial microstructures for crystal plasticity simulations in Moose:

**Step 1:** Generate a polycrystalline structure with Neper by running the following command:
```bash
neper -T -n 20 -dim 2 -morpho graingrowth -domain "square(10,10)"
```
This command generates a polycrystalline structure with a size of 10x10 containing 20 grains.

**Step 2:** Generate a mesh file named `n20-id1.msh` using Neper by running the following command:
```bash
neper -M n20-id1.tess -elttype tri -rcl "body>5?0.35:1"
```
Afterward, make a copy of `n20-id1.msh` and name it `n20-id2.msh`.

**Step 3:** Open the `n20-id2.msh` file using software like Gmsh on Windows 11. Make any necessary modifications if needed and save it, overwriting the existing `n20-id2.msh` file.

**Step 4:** Copy the `PhysicalNames` section from `n20-id1.msh` and paste it into the `n20-id2.msh` file. This step ensures that the material information is consistent in the mesh file.

**Step 5:** In Moose, run the `polycrystal_neper_1.i` input file to perform the crystal plasticity simulation for the uniaxial tension case. Note that this example uses both mesh adaptation and time adaptation during the calculation.

By following these steps, you will be able to generate and utilize Neper-generated polycrystalline structures as initial microstructures for crystal plasticity simulations in Moose.

# 建立晶体塑性模型
1. 2023.09.07 - 基于 `c_pfor_am` 迁移对象 `ComputeElasticityTensorCPGrain` & `GrainPropertyReadFile`:
   1. ComputeElasticityTensorCPGrain - 
   2. GrainPropertyReadFile - 

# 脚本

## 脚本1-获取内存使用情况
```bash 
while true; do
    # 获取当前时间
    current_time=$(date +"%Y-%m-%d %H:%M:%S")
    
    # 获取内存使用情况并提取已使用内存（以G为单位）
    memory_info=$(free -h | awk '/Mem:/ {print $3}')
    
    # 获取总内存大小
    total_memory=$(free -h | awk '/Mem:/ {print $2}')
    
    # 计算已使用内存占总内存的比重
    memory_percentage=$(awk "BEGIN {printf \"%.2f\", ${memory_info}/${total_memory} * 100}")
    
    # 将数据写入memory.txt文件
    echo "${current_time} ${memory_info}G ${memory_percentage}%" >> memory.txt
    
    # 每隔一段时间采集数据，可以根据需要调整时间间隔
    sleep 10  # 10秒钟采集一次数据
done
```

# 后处理
1. `E:\PhD\prm3_GNS_mechinicalStability_coupledModeling_2023_new\p5_threeType_coupledStype_2023\scripts\m1_draw_tensile_curve.m`
   1. 用于绘制不同背应力系数下的晶体塑性模型参数的测试结果输出
2. 
# 问题
1. 内置多晶初始条件，如何加速？
2. 内置多晶初始条件，是否可以使用网格自适应技术？
3. 采用neper创建的多晶结构以及网格文件，计算效率是否会有所提高？


# 理论笔记：CPFEM-MOOSE

