# bw4c — Cython-optimized Bisulfite/TAPS Alignment Pipeline

Cython version for high-performance bisulfite sequencing alignment (C2T conversion, alignment, and back-conversion).

## Changelog

- Fixed an issue where bw4 would mark alignments with low match percentage as failed (MAPQ=20). Since BWA itself does not fail these alignments, bw4 no longer marks them as failed. (2022-04-01)

## Installation

1. Configure Python3 with the required packages.
2. Update the Python path in `make.sh`.
3. Run:

```bash
bash make.sh
```

## Usage

### Step 1: Align FQ reads to Watson and Crick references

Use `simpleScripter.pl` to configure program paths and input file paths, which generates shell scripts for job submission. Each FQ pair is aligned separately to the Watson and Crick reference genomes.

**PE (paired-end) example:**

```bash
# Watson strand
python3 bw4c.py c2t fq1.gz fq2.gz | \
  bwa mem -p -t 12 -T 40 -B 2 -L 10 -CM -Y $REF_DIR/bw4/ref.c2t.watson.fa /dev/stdin | \
  python3 bw4c.py seqBack /dev/stdin | \
  samtools view -T $REF_DIR/bw4/ref.c2t.watson.fa -b -o sample.watson.pe.bam

# Crick strand
python3 bw4c.py c2t fq1.gz fq2.gz | \
  bwa mem -p -t 12 -T 40 -B 2 -L 10 -CM -Y $REF_DIR/bw4/ref.c2t.crick.fa /dev/stdin | \
  python3 bw4c.py seqBack /dev/stdin | \
  samtools view -T $REF_DIR/bw4/ref.c2t.crick.fa -b -o sample.crick.pe.bam
```

**SE (single-end) example:**

```bash
# Watson strand
python3 bw4c.py c2t fq.se.gz NA | \
  bwa mem -p -t 12 -T 40 -B 2 -L 10 -CM -Y $REF_DIR/bw4/ref.c2t.watson.fa /dev/stdin | \
  python3 bw4c.py seqBack /dev/stdin | \
  samtools view -T $REF_DIR/bw4/ref.c2t.watson.fa -b -o sample.watson.se.bam

# Crick strand
python3 bw4c.py c2t fq.se.gz NA | \
  bwa mem -p -t 12 -T 40 -B 2 -L 10 -CM -Y $REF_DIR/bw4/ref.c2t.crick.fa /dev/stdin | \
  python3 bw4c.py seqBack /dev/stdin | \
  samtools view -T $REF_DIR/bw4/ref.c2t.crick.fa -b -o sample.crick.se.bam
```

### Step 2: Merge BAMs and post-processing

```bash
python3 bw4c.py bamMerge watson.bam crick.bam out.bam pe
# or: se (depending on the alignment mode above)
```

Full pipeline example (with markdup):

```bash
python3 bw4c.py bamMerge watson.bam crick.bam /dev/stdout pe | \
  samtools fixmate -m -@ 4 /dev/stdin /dev/stdout | \
  samtools sort -l 0 -O BAM -@ 4 -m 3G /dev/stdin -T sample.fixmate.sort.markdup.bam | \
  samtools markdup -@ 4 /dev/stdin sample.fixmate.sort.markdup.bam
```

> **Note on NM/MD tags:** If you need correct NM and MD tags relative to the *unconverted* reference, run `samtools calmd` as a final step. Without it, NM/MD are computed against the C2T-converted reference. Most downstream tools recompute these tags, but some may read them directly.

```bash
samtools calmd -@ 4 -b -Q input.bam $REF_DIR/ref.fa > calmd.bam
```

### Resources

Memory usage is primarily determined by BWA thread count, samtools thread count, and sort memory.

If `samtools markdup` runs out of memory, use Picard MarkDuplicates instead:

```bash
java -Djava.io.tmpdir=tmp -jar -Xmx6g -Xms6g -Xss256k \
  -XX:-MaxFDLimit -XX:+UseParallelOldGC -XX:ParallelGCThreads=4 \
  picard.jar MarkDuplicates I=sort.bam O=markdup.bam \
  VALIDATION_STRINGENCY=SILENT M=markdup.matrix ASO=coordinate
```

## Known Issues

- Slightly different clip behavior at both ends of reads (fewer S clips) compared to standard BWA. This is related to bw4's penalty parameter settings and is kept as-is for consistency with the original bw4 parameters.
- If `samtools sort` fails to allocate memory, reduce or disable the sort thread parameter.

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE).
