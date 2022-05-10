# 计算逻辑

## 数据准备

分染色体处理,进行质量过滤,查找bam中的成对reads,加上insert size,组成片段fragments包含start end qname barcode,
对fragments数据,去除PCR重复(根据start end barcode唯一性),并统计barcode的次数

## 确定HQ barcodes

根据次数,使用KDE(核密度估计)方法计算骤降点,从而过滤出来HQ barcodes

## 计算reads overlap

对frags数据做过滤,只要barcode符合HQ的,去除高NC值(6),以及去除PCR重复,通过判断start/end相同,将认为是同一序列的不同barcode进行两两组合,统计次数,去掉个数太少(4)的数据

## barcode融合

对两两的barcode组合,计算雅卡尔相似系数,计算数据来源是barcode出现次数以及两者相交的次数,并使用雅卡尔相似系数根据KDE计算的阈值过滤数据

根据两两的barcode组合,统计每个barcode对应的相似barcode,开始融合,将相似的barcode融合为一个drop barcode,形式为run_BCxxxx_Nxxxx,并统计原始barcode与drop barcode的对应关系,多个原始barcode对应到一个新的barcode,表示这些原始barcode在droplet-level barcode是重复的

## 更新bam并输出fragments

将frag的qname删掉,把旧barcode替换成新barcode,并去重,就是要输出的frags

bam文件添加tags,增加一列DB数据

根据用户输入参数计算饱和度

## 统计与绘图

针对细胞核和线粒体,找出两者相交的cell barcode,计算各自的cell barcode的uniq(chr start end cellbarcode)和total指标,
insert size的均值,中值,与tss区域有overlap的比率,frip,dup率,library size等,输出统计文件

绘图包括barcode/jaccard的阈值等

## workflow

![流程图](workflow.png)
