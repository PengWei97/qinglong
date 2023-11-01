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
> moose::next-dev
> E:\PhD\prm3_GNS_mechinicalStability_coupledModeling_2023_new\p5_threeType_coupledStype_2023\moose_建模.pptx
>
> 
> [notes_coupled](./graingrowth/p4_CPPF_threeCoupled2023/coupled_tests/notes_coupled.md)

### 学习资料
1. [Compute Finite Strain Elastic Stress](https://mooseframework.inl.gov/source/materials/crystal_plasticity/ComputeMultipleCrystalPlasticityStress.html)
2. Input file:
   1. /home/pw-moose/projects/moose/modules/combined/test/tests/ACGrGrElasticDrivingForce/bicrystal.i
   2. 

### 建模流程
2. 考虑背应力的晶体塑性模型 (后续加入)：
   1. materials：ComputeElasticityTensorCPPF
      1. 借鉴 - ComputeElasticityTensorCPGrain & ComputePolycrystalElasticityTensor
4. VectorPostprocessor: FeatureMatePropVectorPostprocessor <-- FeatureVolumeVectorPostprocessorCopy
   1. TIMESTEP_BEGIN 时统计每个晶粒的平均材料参数
5. UserObjects: GrainTrackerMatProp <-- GrainTrackerData + GrainTracker
6. Materials: ComputePolyMultCPStressCpl <-- ComputePolyMultCPStressCopy + ComputeMultipleCrystalPlasticityStress + computePolyMultCPElasticStressCpl
7. Actions: PolyElasticEnergyDrivingCpl <-- PolycrystalElasticDrivingForce
8. Kernels: ACGGElasticEnergyCpl <-- ACGrGrElasticDrivingForce
9. input_file: coupled_elastic_energy_cp.i ~ 只考虑弹性能3
10. 结果测试 ~ 失败
11. 查看每个对象调用的顺序
    1.  初始时刻-初始化（构造函数），先调用 GrainTrackerMatProp ，之后执行材料类： GBEvolution ，ComputePolyElasticTensorCpl ， ComputePolyMultCPStressCplV2 个三次，之后执行 FeatureMatePropVectorPostprocessorV2 ，最后kernels： ACGrGrPoly  +  ACGrGrElasticDrivingForce 各5次
    2.  初始化之后， GrainTrackerMatProp 在设定 time_begin 执行一次，包括函数 `initialize execute finalize`.
    3.  计算过程中，每个线性求解和非线性求解中都存在 computeQpElasticityTensor() 等
    4.  判定该[_qp]存在新激活的晶粒，修改 GrainTracker, 在 GrainTracker中添加 _entity_var_to_features_old
    5.  
12. TODO-
    1.  在 CrystalPlasticityKalidindiUpdate 对 _slip_resistance 进行插值处理
    2.  结果测试： 去掉弹性能驱动力试试
    3.  
    4.  step 1: 对CPFEM计算所得的 Sigma， E^e 进行插值处理之后，再计算残差
    5.  step 2：需要修改GrainTracker来确定每个[_qp]处新激活的晶粒，并设定CPFEM的内部状态变量
~~~~


# 脚本
cp /home/pw-moose/projects/moose/modules/phase_field/src/kernels/ACGrGrElasticDrivingForce.C /home/pw-moose/projects/qinglong/src/kernels/ACGGElasticEnergyCpl.C

cp /home/pw-moose/projects/moose/modules/phase_field/include/kernels/ACGrGrElasticDrivingForce.h /home/pw-moose/projects/qinglong/include/kernels/ACGGElasticEnergyCpl.h

mkdir /home/pw-moose/projects/qinglong/src/kernels/
mkdir /home/pw-moose/projects/qinglong/include/kernels/

tar -jcvf 

1. ll ex_t1_gg_elastic_cp/*.e-s*.20
2. tar -cvf - ex_t1_gg_elastic_cp/*.e.* ex_case4_recovery_v42/*.e-s0095.* ex_case4_recovery_v42/*.e-s0045.* | pigz -9 -p 20 > ex_case4_recovery_v4_v42.tgz
3. tar -cvf - ex_t2_gg_elastic/* | pigz -9 -p 20 > ex_t2_gg_elastic.tgz

# 其他

1. ComputePolycrystalElasticityTensor ~ ComputeElasticityTensorBase ~ DerivativeMaterialInterface<Material>
2. ComputeElasticityTensorCP ~ ComputeElasticityTensor ~ ComputeRotatedElasticityTensorBaseTempl<is_ad> ~ ComputeElasticityTensorBaseTempl<is_ad> ~ DerivativeMaterialInterface<Material>



GrainPropertyReadFile

ComputeElasticityTensorCPGrain