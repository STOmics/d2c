# D2C

Drop to Cell

目前版本 **V1.0.0**

## 环境

### 运行环境

* centos 7.0+
* gcc-9.1 library
* R library

### 脚本设置

需引入gcc动态库和R的可执行文件路径,形式如下:

```sh
export LD_LIBRARY_PATH="/hwfssz5/ST_BIGDATA/USER/zhaofuxiang/lib/gcc-9.1.0/lib:/hwfssz5/ST_BIGDATA/USER/zhaofuxiang/lib/gcc-9.1.0/lib64:$LD_LIBRARY_PATH"
export PATH="/hwfssz1/ST_MCHRI/STEMCELL/USER/yuanyue/bin:$PATH"
```

## 输入/输出

### 输入
linux命令行运行软件 **./bin/d2c -h** 查看参数解释

```
$./bin/d2c -h
D2C: Drop to Cell.
Usage: ./install/bin/d2c [OPTIONS]

Options:
  -h,--help                             Print this help message and exit
  -i TEXT:FILE REQUIRED                 Input bam filename
  -o TEXT REQUIRED                      Output result path
  -b TEXT:FILE                          Barcode list file
  --bt TEXT                             Barcode tag in bam file, default 'XB'
  --mapq INT                            Filter thrshold of mapping quality, default 30
  -c INT                                CPU core number, default detect
  -n TEXT                               Name for the all output files, default prefix of input bam file
  --tn5                                 Process data knowing that the barcodes were generated with a barcoded Tn5
  --bf FLOAT                            Minimum number of fragments to be thresholded for doublet merging
  --ji FLOAT                            Minimum jaccard index for collapsing bead barcodes to cell barcodes
  --log TEXT                            Set logging path, default is './logs'
  -r TEXT                               Specify supported reference genome, default hg19
  --mc TEXT                             Name of the mitochondrial chromosome
  --bg TEXT:FILE                        Bedtools genome file
  --bl TEXT:FILE                        Blacklist bed file
  --ts TEXT:FILE                        Path bed file of transcription start sites
  --bp FLOAT:FLOAT in [0 - 1]           Percentage of minimum number of fragments to be thresholded for doublet merging
  --jp FLOAT:FLOAT in [0 - 1]           Percentage of minimum jaccard index for collapsing bead barcodes to cell barcodes
  --sat                                 Output sequencing saturation file, default False
  --br TEXT:FILE                        Barcode runname list file, default detect

D2C version: 1.1.0
```

三个必需参数:

* -i 输入bam文件,如无index,程序自动创建index文件
* -o 输出目录
  
可选参数:

* -b barcode列表文件,即存储了1536种barcode类型的文本文件
* --bt bam文件中barcode的标签,默认 *XB*
* --mapq 质量过滤阈值,默认 30
* -c 设置CPU核数,暂不支持该参数
* -n 设置run name,默认使用bam的文件名
* --tn5 设置高通量模式
* --bf 设置过滤阈值,用于过滤fragments
* --ji 设置过滤阈值,用于过滤jaccard index
* --bp 设置过滤百分比,用于过滤fragments
* --jp 设置过滤百分比,用于过滤jaccard index
* -r 设置参考基因组,默认 *hg19*
* --mc 设置线粒体染色体
* --bg --bl --ts 这三个参数用于非模式生物
* --sat 计算测序饱和度
* --br 存储barcode后缀的run name的列表文件(如果未指定,程序需花费一定时间遍历bam文件来查找所有run name)
* --logs 设置日志存储路径,默认为当前路径下的 *logs*

### 输出

假设 run name 为 *ABC*, 正常情况在设置的 *-o* 路径输出以下文件:

数据文件:
* ABC.bam
* ABC.barcodeQuantSimple.csv 
* ABC.barcodeTranslate.tsv
* ABC.fragments.tsv.gz
* ABC.HQbeads.tsv
* ABC.implicatedBarcodes.csv.gz

统计文件:
* ABC.NCsumstats.tsv
* ABC.QCstats.csv
* ABC.basicQC.tsv
* ABC.d2cParam.csv
* ABC.sequenceSaturation.tsv 测序饱和度输出文件,只有给定 *--sat* 选项才生成,共四列,分别是采样比率,每个细胞的平均fragment个数,对应的测序饱和度,对应的每个cell下的唯一barcode的中值

图表文件:
* ABC.beadBarcodeKneeCurve.pdf
* ABC.beadBarcodeKneeDensity.pdf
* ABC.beadBarcodeKnee.pdf
* ABC.beadBarcodeKnee.png
* ABC.jaccardOverlapKnee.pdf
* ABC.jaccardOverlapKnee.png
* ABC.kneesPlotted.txt

日志文件:
* bin/logs/D2C_20200813_140243.log 日志在程序目录下的logs文件夹,按程序启动时间建立文件名


## 示例

模式生物

除了三个必须参数之外,指定了barcode 标签为 *CB*, run name为 *ABC*, 过滤质量为 *20*, 参考基因组为 *mm10*, barcode后缀的run name列表文件 *runname.list* 

```
./bin/d2c \
    -i input.bam \
    -o output \
    -b barcode.list \
    --bt CB \
    -n ABC \
    --mapq 20 \
    -r mm10 \
    --br runname.list
```

非模式生物

与模式生物不同之处就是设置了 *--bg --bl --ts* 三个参数

```
./bin/d2c \
    -i input.bam \
    -o output \
    --mc chrMT \
    --bt CB \
    -n ABC \
    --bg ABC.sizes \
    --bl ABC.bed \
    --ts ABC.bed \
    --mapq 30 \
    -b barcode.list \
    --br runname.list
```

## TODO

功能:

* 支持 *-c* 参数指定CPU核数,程序运行更加灵活
* 使用python画图替代现有的R脚本,图表形式变更
* 增加测序饱和度的图片输出


性能:

* 降低大数据量的内存消耗
* 提高程序运行速度