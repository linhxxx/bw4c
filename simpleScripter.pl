my $fq1 = shift;
my $fq2 = shift;
my $outputPrefix = shift;
my $bwa = shift;
my $samtools = shift;
my $python = shift;
my $bw4c = shift;
my $refwat = shift;
my $refcrk = shift;
my $pe = shift;

if ($pe){
   print "$python $bw4c c2t $fq1 $fq2 | $bwa mem -p -t 12 -T 40 -B 2 -L 10 -CM -Y $refwat /dev/stdin | $python $bw4c seqBack /dev/stdin | $samtools view -T $refwat -b -o  $outputPrefix.watson.pe.bam\n";
   print "$python $bw4c c2t $fq1 $fq2 | $bwa mem -p -t 12 -T 40 -B 2 -L 10 -CM -Y $refcrk /dev/stdin | $python $bw4c seqBack /dev/stdin | $samtools view -T $refcrk -b -o  $outputPrefix.crick.pe.bam\n";
}
else{
    print "$python $bw4c c2t $fq1 NA | $bwa mem  -t 12 -T 40 -B 2 -L 10 -CM -Y $refwat /dev/stdin | $python $bw4c seqBack /dev/stdin | $samtools view -T $refwat -b -o  $outputPrefix.watson.se.bam\n";
    print "$python $bw4c c2t $fq1 NA | $bwa mem  -t 12 -T 40 -B 2 -L 10 -CM -Y $refcrk /dev/stdin | $python $bw4c seqBack /dev/stdin | $samtools view -T $refcrk -b -o  $outputPrefix.crick.se.bam\n";
}

