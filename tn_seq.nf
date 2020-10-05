#!/usr/bin/env nextflow

/****************
**
** NEXTFLOW TN-Seq pipeline.
**
** This pipeline is developed for use at MIMS.
*****************/

// FASTA genome sequence
params.refseq = ["NZ_AP014524.1", "NZ_AP014525.1"]

// FASTA sequence of the TA plasmid. Used to check for contamination of plasmid
// sequence
params.plasmid = ["AY115560.1"]

params.filter_tasites = "gene"

// Minimum memory to be allocated for tasks. Lower will be faster, but at a risk
// of jobs crashing. Should probably not ever be lower than 256m
params.min_memory = "1g"

params.max_memory = "4g"

params.map_cpus = 6
params.trim_cpus = 4
params.duk_cpus = 4

/*************
**
** END OF USER OPTIONS
** DO NOT CHANGE ANYTHING BELOW THIS LINE UNLESS YOU KNOW WHAT YOU ARE DOING
**
**************/


process verify_sampleinfo {
  input:
  path sample_file from "$baseDir/datafiles/sampledata.csv"
  path rawdata from "$baseDir/rawdata/"
  path verify_script from "$baseDir/scripts/verify_samples.py"

  output:
  val true into verify_ch

  """
  #!/bin/bash

  python3 ${verify_script}
  """ 
}

@Grab('com.xlson.groovycsv:groovycsv:1.1')
import static com.xlson.groovycsv.CsvParser.parseCsv

file('./rawdata').mkdir()
def sample_data = []
def sample_file = file('./datafiles/sampledata.csv')
def data = parseCsv(sample_file.text, autoDetect:true)
def record = 0
for (line in data) {
  def sample = ""
  def name = ""
  def id = ""
  if ("$line.Sample" == "" && "$line.BasespaceId" == "") {
    println "You must provide either 'Sample' or 'BasespaceId' for each row"
  }

  name = "$line.Name"
  id = "$line.BasespaceId"
  if ("$line.Sample" == "") {
    sample = "$line.Name" + ".fastq.gz"
  } else {
    sample = "$line.Sample"
  }
  sample_data.add([record, sample, name, id])
  record += 1
}
sampleinfo_ch = Channel.fromList(sample_data)

process get_samples {
  input:
  val sampledata from sampleinfo_ch
  val check_samples from verify_ch
  path rawdata from "$baseDir/rawdata/"

  output:
  path '*.fastq.gz' into reads_ch_1
  path '*.fastq.gz' into reads_ch_2
  path '*.json' optional true

  """
  #!/bin/bash

  ID=\"${sampledata[3]}\"
  SN=\"${sampledata[1]}\"

  if [ ! -z \$ID ]; then
    bs download dataset -i \"\$ID\" -o .
    mv -f *.fastq.gz \"\$SN\"
  else
    x=\$(find -L ./rawdata -name \"\$SN\" -type f)
    ln -sf \"\$x\" .
  fi
  """
}

process genome_download {
  cpus 1
  memory params.min_memory

  input:
  val glist from params.refseq.size > 1 ? params.refseq.join(',') : params.refseq[0]

  output:
  path "genome.fa" into genome_ch

  publishDir 'genome'

  """
  #!/bin/sh
  curl 'https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=fasta&id=${glist}' > genome.fa
  """
}

process annotation_download {
  cpus 1
  memory params.min_memory

  input:
  val glist from params.refseq.size > 1 ? params.refseq.join(',') : params.refseq[0]

  output:
  path "genome.gff3" into gff_ch

  publishDir 'genome'

  """
  #!/bin/sh
  curl 'https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=${glist}' > genome.gff3
  """
}

process plasmid_download {
  cpus 1
  memory params.min_memory

  input:
  val glist from params.plasmid.size > 1 ? '%5B' + params.plasmid.join('%2C') + '%5D' : params.plasmid[0]

  output:
  path "plasmid.fa" into plasmid_ch

  publishDir 'plasmid'

  """
  #!/bin/sh
  curl 'https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=fasta&id=${glist}' > plasmid.fa
  """
}


// Module to create a genome index. Uses BBMap
process genome_index {
  cpus params.map_cpus
  memory params.max_memory

  input:
  path genome_fasta from genome_ch
  val cpus from params.map_cpus
  val mem from params.max_memory

  output:
  tuple path(genome_fasta),path('ref') into genome_idx_ch

  """
  #!/bin/sh
  bbmap.sh -Xmx${mem} trd ref=$genome_fasta vslow=t threads=${cpus}
  """
}

// First QC using FastQC on raw data
process fastqc1 {
  cpus 1
  memory params.min_memory

  input:
  path reads from reads_ch_1

  output:
  path '*.zip'

  publishDir 'analysis/00-qc-raw'

  """
  #!/bin/sh

  /usr/local/bin/FastQC/fastqc -o . --noextract $reads
  """
}

// Trim common adapter sequences and very bad quality bases using trim_galore
process trim_galore {
  cpus params.trim_cpus
  memory params.min_memory

  input:
  path reads from reads_ch_2
  val cpus from params.trim_cpus

  output:
  path '*trimming_report.txt'
  path '*trimmed.fq.gz' into trimmed_ch_1
  path '*trimmed.fq.gz' into trimmed_ch_2

  publishDir 'analysis/01-qc-trimmed'

  """
  #!/bin/sh

  trim_galore -j ${cpus} -q 5 --gzip $reads 
  """
}

// Second QC (post trimming) using FastQC on raw data
process fastqc2 {
  cpus 1
  memory params.min_memory

  input:
  path reads from trimmed_ch_1

  output:
  path '*.zip'

  publishDir 'analysis/01-qc-trimmed'

  """
  #!/bin/sh

  /usr/local/bin/FastQC/fastqc -o . --noextract $reads
  """
}

// Detect and remove plasmid sequences (failed selection) using BBDuk 
process bbduk {
  cpus params.duk_cpus
  memory params.min_memory

  input:
  path reads from trimmed_ch_2
  path plasmid from plasmid_ch
  val mem from params.min_memory
  val cpus from params.duk_cpus

  publishDir 'analysis/02-qc-decontaminate'

  output:
  path '*_dec.fq.gz' into duk_ch_1
  path '*_dec.fq.gz' into duk_ch_2
  path '*bbduk*'

  """
  I=$reads
  bbduk.sh -Xmx${mem} in=$reads out=\${I/.fq.gz/_dec.fq.gz} ref=$plasmid \
    stats=${reads}.bbduk.stats refstats=${reads}.bbduk.refstats k=23 ktrim=r \
    mink=11 threads=${cpus}
  """

}

// Third QC (post decontamination) using FastQC on raw data
process fastqc3 {
  cpus 1
  memory params.min_memory

  input:
  path reads from duk_ch_1

  output:
  path '*.zip'

  publishDir 'analysis/02-qc-decontaminate'

  """
  #!/bin/sh

  /usr/local/bin/FastQC/fastqc -o . --noextract $reads
  """
}

// Align using BBMap, very slow & accurate settings
process bbmap {
  cpus params.map_cpus
  memory params.min_memory

  input:
  path reads from duk_ch_2
  path genome from genome_idx_ch
  val mem from params.min_memory
  val cpus from params.map_cpus

  output:
  tuple path(genome),path('*.sam') into map_ch
  path '*.bbmap.*'

  publishDir 'analysis/03-align'

  """
  #!/bin/bash

  bbmap.sh -Xmx${mem} ref=${genome[0]} in=$reads out=${reads}.sam vslow=t \
    covstats=${reads}.bbmap.covstats statsfile=${reads}.bbmap.stats \
    threads=${cpus} trd
  """
}

// Convert SAM to BAM, get more alignment quality stats from samtools
process samtools {
  cpus 1
  memory params.min_memory

  input:
  tuple path(genome),path(aln) from map_ch

  output:
  tuple path('*.bam'),path('*.bai') into samtools_ch_1
  tuple path('*.bam'),path('*.bai') into samtools_ch_2
  path '*.flagstat'
  path '*.idxstats'
  path '*.stats'
  

  publishDir 'analysis/03-align'

  """
  #!/bin/bash

  samtools view -bT ${genome[0]} $aln | samtools sort -o ${aln}.bam -O BAM
  samtools index ${aln}.bam
  samtools flagstat ${aln}.bam > ${aln}.flagstat
  samtools idxstats ${aln}.bam > ${aln}.idxstats
  samtools stats ${aln}.bam > ${aln}.stats
  """
}

// Final all-encompassing preprocessing QC report from MultiQC
process multiqc {
  cache false
  cpus 1
  memory params.min_memory

  input:
    path dummy from samtools_ch_1.collect()

  output:
    path 'multiqc_report.html'

  publishDir 'analysis/'

  """
  multiqc -d $baseDir/analysis
  """
}

process tasites {
  cpus 1
  memory params.min_memory

  input:
  path genome from genome_ch

  output:
  path 'genome_tasites.bed' into tasites_ch

  publishDir 'genome/'

  """
  #!/bin/sh
  TASites --one --ref ${genome} > genome_tasites.bed
  """
}

process filter_ta {
  cpus 1
  memory params.min_memory

  input:
  path gff from gff_ch
  val feat from params.filter_tasites

  output:
  path 'genome_filtered.gff3' into gff_filtered_ch

  publishDir 'genome/'

  """
  #!/bin/sh
  awk '{ if (\$3 == \"${feat}\") { print } }' ${gff} > genome_filtered.gff3
  """
}

process overlap_ta {
  cpus 1
  memory params.min_memory

  input:
  path gff from gff_filtered_ch
  path tas from tasites_ch

  output:
  path 'overlap_bed.tsv' into overlap_ch

  """
  #!/bin/sh

  bedtools intersect -loj -a ${tas} -b ${gff} > overlap_bed.tsv
  """
}

process generate_saf {
  cpus 1
  memory params.min_memory

  input:
  path tsv from overlap_ch

  output:
  path 'genes_tasites.saf'

  publishDir 'genome/'

  """
  #!/bin/sh
  cat ${tsv} | \
  tr ';' '\t' | \
  awk '{ gsub(/ID=/,\"\"); \
         print \$12 \"\\t\" \$1 \"\\t\" \$2 \"\\t\" \$3 \"\\t\" \$10}' | \
  grep -E -v '^\\.' > genes_tasites.saf
  """
}