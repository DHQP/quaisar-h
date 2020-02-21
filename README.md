README for quaisar

How to run: (updated 3/21/19)  
Project:
```bash
qsub ./parallel_quaisar.sh -p run_ID
```
List
```bash
qsub ./parallel_quaisar.sh -i /path/to/folder/containing/zipped/fastqs [1,2,3,4] -o /path/to/output/directory/ outputFolderName

[1,2,3,4] - options that correlate to postfix files types (1=_L001_SX_RX_00X_fastq.gz. 2=_RX.fastq.gz. 3=X.fastq.gz 4=_SX_RX_00X.fastq.gz)"
*don't include brackets in command
```
