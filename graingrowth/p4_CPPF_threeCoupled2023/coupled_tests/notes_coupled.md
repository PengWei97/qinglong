> 耦合模型构建的笔记
>
# TODO
1. 解码考虑弹性能的GG
2. 解码multi-APP的逻辑思路.
3. 创建自定义类 GrainTrackerMatProp 

# 建模计划
1. 基于 ComputeMultipleCrystalPlasticityStress, CrystalPlasticityKalidindiUpdate 代码创建 ComputePolycrystalCrystalPlasticityStress & PloycrystalCrystalPlasticityKalidindiUpdate ；
   1. https://mooseframework.inl.gov/source/materials/crystal_plasticity/ComputeMultipleCrystalPlasticityStress.html
2. 简化 Kalidindi cp model，不考虑塑性部分 - PolycrystalNoPlasticityUpdate；；
   1. : public CrystalPlasticityStressUpdateBase
3. 创建测试算例 test_PolycrystalNoPlasticityUpdate.i


# 建模思路
1. 简化晶体塑性模型，去掉塑性部分的模拟
2. 只考虑弹性部分的插值，不考虑塑性部分


# refs
## moose-web
1. https://mooseframework.inl.gov/modules/tensor_mechanics/

## c++ objects
1. ComputeMultipleCrystalPlasticityStress, CrystalPlasticityKalidindiUpdate, ComputeElasticityTensorCPGrain
2. GrainTrackerElasticity, ComputePolycrystalElasticityTensor, 

# 脚本
cp /home/pw-moose/projects/moose/modules/phase_field/src/postprocessors/GrainTrackerElasticity.C ~/projects/qinglong/src/postprocessors/GrainTrackerMatProp.C
cp /home/pw-moose/projects/moose/modules/phase_field/include/postprocessors/GrainTrackerElasticity.h ~/projects/qinglong/include/postprocessors/GrainTrackerMatProp.h

code  /home/pw-moose/projects/moose/modules/phase_field/include/vectorpostprocessors/FeatureVolumeVectorPostprocessor.h

code  ~/projects/qinglong/include/vectorpostprocessors/FeatureMatePropVectorPostprocessor.h
code  ~/projects/qinglong/src/vectorpostprocessors/FeatureMatePropVectorPostprocessor.C