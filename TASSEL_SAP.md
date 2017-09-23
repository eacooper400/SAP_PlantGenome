# Command line history from TASSEL 5 on SAP

## Tag Counts
```bash
lizcooper$ mkdir TagCounts
lizcooper$ /Users/lizcooper/tassel5.0_standalone/run_pipeline.pl -Xmx50G -fork1 -FastqToTagCountPlugin -i /Volumes/Kresovich/DataArchives/DNA/GBS/CornellRaw/BREAD/ -k ../SAP_13_14_key.txt -e ApeKI -s 700000000 -c 1 -o TagCounts/ -endPlugin -runfork1 >TagCounts/FastqToTagCount.log 2>TagCounts/FastqToTagCount.err
```

## Merge Multiple Tag Counts
```bash
lizcooper$ mkdir mergedTagCounts
lizcooper$ /Users/lizcooper/tassel5.0_standalone/run_pipeline.pl -fork1 -MergeMultipleTagCountPlugin -Xmx32g -i TagCounts/ -o mergedTagCounts/MasterSAPtags.cnt -c 10 -endPlugin -runfork1 >mergedTagCounts/MergeMultipleTags.log 2>mergedTagCounts/MergeMultipleTags.err
lizcooper$ /Users/lizcooper/tassel5.0_standalone/run_pipeline.pl -fork1 -TagCountToFastqPlugin -Xmx32g -i mergedTagCounts/MasterSAPtags.cnt -o mergedTagCounts/MasterSAPtags.fq -c 10 -endPlugin -runfork1 >mergedTagCounts/TagsToFastq.log 2>mergedTagCounts/TagsToFastq.err
```

## Tag Alignment

```bash
lizcooper$ mkdir bwa_alignment
lizcooper$ /Users/lizcooper/bwa-0.7.8/bwa aln -t 4 /Users/lizcooper/Sorghum_Genome/Sbicolor_v2.1_255.fa mergedTagCounts/MasterSAPtags.fq >bwa_alignment/mergedSAPtags.sai
lizcooper$ /Users/lizcooper/bwa-0.7.8/bwa samse /Users/lizcooper/Sorghum_Genome/Sbicolor_v2.1_255.fa bwa_alignment/mergedSAPtags.sai mergedTagCounts/MasterSAPtags.fq >bwa_alignment/mergedSAPtags.sam
lizcooper$ sed 's/Chr0//g' mergedSAPtags.sam | sed 's/Chr//g' | sed 's/super_/1/g' >mergedSAPtags_renamed.sam
```

## TOPM
```bash
lizcooper$ mkdir topm
lizcooper$ /Users/lizcooper/tassel5.0_standalone/run_pipeline.pl -fork1 -SAMConverterPlugin -i bwa_alignment/mergedSAPtags_renamed.sam -o topm/MasterSAPtags.topm -endPlugin -runfork1 >topm/SAMConverter.log 2>topm/SAMConverter.err
```

## TBT
```bash
lizcooper$ mkdir tbt
lizcooper$ /Users/lizcooper/tassel5.0_standalone/run_pipeline.pl -fork1 -FastqToTBTPlugin -i /Volumes/Kresovich/DataArchives/DNA/GBS/CornellRaw/BREAD/ -k ../SAP_13_14_key.txt -e ApeKI -o tbt/ -y -t mergedTagCounts/MasterSAPtags.cnt -endPlugin -runfork1 >tbt/FastqToTBT.log 2>tbt/FastqToTBT.err
lizcooper$ /Users/lizcooper/tassel5.0_standalone/run_pipeline.pl -fork1 -MergeTagsByTaxaFilesPlugin -Xmx32g -i tbt/ -o tbt/mergedSAP.tbt.byte -endPlugin -runfork1 >tbt/MergeTagsBYTaxa.log 2>tbt/MergeTagsByTaxa.err
```

##  SNP Calling
```bash
lizcooper$ mkdir hapmap
lizcooper$ /Users/lizcooper/tassel5.0_standalone/run_pipeline.pl -Xmx10g -fork1 -DiscoverySNPCallerPlugin -i tbt/mergedSAP.tbt.byte -m topm/MasterSAPtags.topm  -o hapmap/SAP_chr+.hmp -sC 1 -eC 10 -mnMAF 0.01 -endPlugin -runfork1 >vcf/SNPCaller.log 2>vcf/SNPCaller.err
```
