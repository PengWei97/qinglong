qinglong
=====

"Fork qinglong" to create a new MOOSE-based application.

For more information see: [http://mooseframework.org/create-an-app/](http://mooseframework.org/create-an-app/)

# Project Plan
> moose::next-new 是完全和idaholab/moose::next分支代码完全相同的分支
> moose::next-dev 是配合qinglong::next需要同步合并的分支
> qinglong::next 是目前演进版本的分支，CuAGG2022/TiGNS2022/CPPF2023完成代码构建之后创建RP的分支


## CuAGG2022
> Refactor the code in the article [acta_copperThinFilm_2023](https://www.sciencedirect.com/science/article/pii/S1359645423005669?via%3Dihub)

### 已完成
1. 考虑GB anisotropy ~ `GBAnisotropyMisori` and `MisorienationAngleCalculator`;
2. 添加输出相邻晶粒数目及特征ID ~ `FeatureDataVectorPostprocessor` and  `FeatureFloodCount`;
3. 测试算例 `GBAnisotropyMisori`

### TODO
1. 文章相关部分的代码及算例代码的重构
2. `GBAnisotropyMisoriInit` 材料类的重构

## TiGNS2022
> 用于研究梯度结构热稳定性的分支

### 已完成
1. 为了完成当相邻晶粒取向差低于某个阈值时自动合并的操作：
   1. 修改了 `GrainTracker`，创建了两个钩子 ~ moose::TiGNS2022；
   2. 创建了派生类 `GrainTrackerMerge`来具体化判定自动合并的标准；
2. 合并
   1. 将 `GrainTracker` & `GrainTrackerMerge` 从 moose::TiGNS2022 合并到 moose::next-dev中


### TODO
2. 添加考虑存储能的晶粒长大模块：
   1. 材料类: DeformedGrainEBSDMaterial, 
   2. kernel类: ACSEDGPolyEBSD
   3. action类：PolycrystalStoredEnergyEBSD
3. 相关算例
   1. 550/700du

## CPFP2023
1. 考虑背应力的晶体塑性模型：
   1. materials：

