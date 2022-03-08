## ONT R9.4.1 Guppy 5.0.7 "Sup" HG003 whole genome performance evaluation
We hold out `HG003` sample while training `PEPPER-Margin-DeepVariant` so we use `HG003` to demonstrate our whole genome performance. We report both runtime and accuracy for this evaluation between PEPPER-Margin-DeepVariant r0.7 and r0.8.

### Setup
We used the following dataset:
```
Sample:                   HG003 (Whole genome)
Coverage:                 ~85x
Chemistry:                R9.4.1
Basecaller:               Guppy 5.0.7 "Sup"
```
#### Downsampling the alignment file:
We downsampled the `~85x` variant calling data using the following command:
```bash
samtools view -s 0.71 -b -@${THREADS} HG003_guppy_507_2_GRCh38_pass.bam > HG003_guppy_507_2_GRCh38_pass.60x.bam
samtools view -s 0.36 -b -@${THREADS} HG003_guppy_507_2_GRCh38_pass.bam > HG003_guppy_507_2_GRCh38_pass.30x.bam
```
#### Calling variants with PEPPER-Margin-DeepVariant (PEPPER r0.8)

```bash
time docker run -it -v /data:/data \
-u `id -u`:`id -g` \
kishwars/pepper_deepvariant:r0.8 \
run_pepper_margin_deepvariant call_variant \
-b $BAM \
-f $REF \
-o $OUTPUT_DIR \
-t $THREADS \
-s HG003 \
--ont_r9_guppy5_sup 2>&1 | tee $LOG_FILE
```

### Results

In all stratified coverages `(30x, 60x, 85x)`, PEPPER-Margin-DeepVariant r0.8 shows increased accuracy:
<p align="center">
<img src="img/pepper_r8_ont_HG003_wgs.png" alt="PEPPER performance whole genome">
</p>

##### HG003 30x performance:
<p align="center">
<table><thead><tr><th>Sample</th><th>Version</th><th>Type</th><th>Truth<br>total</th><th>True<br>positives</th><th>False<br>negatives</th><th>False<br>positives</th><th>Recall</th><th>Precision</th><th>F1-Score</th></tr></thead><tbody><tr><td rowspan="4">HG003 30x</td><td rowspan="2">r0.7</td><td>INDEL</td><td>504501</td><td>317621</td><td>186880</td><td>35084</td><td>0.629575</td><td>0.902714</td><td>0.7418</td></tr><tr><td>SNP</td><td>3327495</td><td>3310002</td><td>17493</td><td>11986</td><td>0.994743</td><td>0.996393</td><td>0.995567</td></tr><tr><td rowspan="2">r0.8</td><td>INDEL</td><td>504501</td><td>337206</td><td>167295</td><td>53674</td><td>0.668395</td><td>0.865863</td><td>0.754422</td></tr><tr><td>SNP</td><td>3327495</td><td>3313043</td><td>14452</td><td>12451</td><td>0.995657</td><td>0.996257</td><td>0.995957</td></tr></tbody></table>
</p>

##### HG003 60x performance:
<p align="center">
<table><thead><tr><th>Sample</th><th>Version</th><th>Type</th><th>Truth<br>total</th><th>True<br>positives</th><th>False<br>negatives</th><th>False<br>positives</th><th>Recall</th><th>Precision</th><th>F1-Score</th></tr></thead><tbody><tr><td rowspan="4">HG003 60x</td><td rowspan="2">r0.7</td><td>INDEL</td><td>504501</td><td>366144</td><td>138357</td><td>33484</td><td>0.725755</td><td>0.91827</td><td>0.810741</td></tr><tr><td>SNP</td><td>3327495</td><td>3317492</td><td>10003</td><td>8548</td><td>0.996994</td><td>0.99743</td><td>0.997212</td></tr><tr><td rowspan="2">r0.8</td><td>INDEL</td><td>504501</td><td>390595</td><td>113906</td><td>47118</td><td>0.77422</td><td>0.895066</td><td>0.830269</td></tr><tr><td>SNP</td><td>3327495</td><td>3318785</td><td>8710</td><td>9212</td><td>0.997382</td><td>0.997233</td><td>0.997308</td></tr></tbody></table>
</p>

##### HG003 85x performance:
<p align="center">
<table><thead><tr><th>Sample</th><th>Version</th><th>Type</th><th>Truth<br>total</th><th>True<br>positives</th><th>False<br>negatives</th><th>False<br>positives</th><th>Recall</th><th>Precision</th><th>F1-Score</th></tr></thead><tbody><tr><td rowspan="4">HG003 85x</td><td rowspan="2">r0.7</td><td>INDEL</td><td>504501</td><td>383384</td><td>121117</td><td>30595</td><td>0.759927</td><td>0.927982</td><td>0.835588</td></tr><tr><td>SNP</td><td>3327495</td><td>3318437</td><td>9058</td><td>8032</td><td>0.997278</td><td>0.997586</td><td>0.997432</td></tr><tr><td rowspan="2">r0.8</td><td>INDEL</td><td>504501</td><td>409096</td><td>95405</td><td>40539</td><td>0.810892</td><td>0.912201</td><td>0.858568</td></tr><tr><td>SNP</td><td>3327495</td><td>3319475</td><td>8020</td><td>8449</td><td>0.99759</td><td>0.997462</td><td>0.997526</td></tr></tbody></table>
</p>
