task Picard_CreateSequenceDictionary {
  File ref_fasta
  String output_filename
  
  command {
    java -Xmx4g -jar /home/biodocker/bin/picard.jar \
    CreateSequenceDictionary \
    REFERENCE=${ref_fasta} \
    OUTPUT=${output_filename}
  }
  output {
    File ref_dict = output_filename
  }
  runtime {
    docker: "fzkhan/picard-1.136-gatk-2.8"
  }
}

task BWA_mem {
  File ref_fasta
  File ref_bwt
  File ref_sa
  File ref_amb
  File ref_ann
  File ref_pac
  Array[File] reads
  String output_sam_name
  String rg
  
  command {
    bwa mem -v 3 -R '${rg}' -t 1 ${ref_fasta} ${sep=' ' reads} > ${output_sam_name}
  }
  output {
    File output_sam = output_sam_name
  }
  
  runtime {
    docker: "scidap/bwa:v0.7.12"
  }
  
}

task samtools_view {
  File input_sam
  String output_filename
  
  command {
    samtools view -Sb ${input_sam} > ${output_filename}
  } 
  output {
    File output_bam = output_filename
  }
  runtime {
    docker: "isaacliao/samtools-0.1.19"
  }
}

task samtools_sort {
  File input_bam
  String output_filename
  
  command {
    samtools sort -f ${input_bam} ${output_filename}
  } 
  output {
    File output_bam = output_filename
  }
  runtime {
    docker: "isaacliao/samtools-0.1.19"
  }
}

task Picard_MarkDuplicates {
  File input_bam
  String output_filename
  String output_index_filename
  String metrics_filename
  
  command {
    java -Xmx4g -jar /home/biodocker/bin/picard.jar \
    MarkDuplicates \
    INPUT=${input_bam} \
    OUTPUT=${output_filename} \
    METRICS_FILE=${metrics_filename} \
    ASSUME_SORTED=true \
    REMOVE_DUPLICATES=true \
    CREATE_INDEX=true
  }
  output {
    File metrics = metrics_filename
    File output_bam = output_filename
    File output_bam_index = output_index_filename
  }
  runtime {
    docker: "fzkhan/picard-1.136-gatk-2.8"
  }
}

task GATK_RealignerTargetCreator {
  File input_bam
  File input_bam_index
  String output_filename
  File reference_fasta
  File reference_fasta_index
  File reference_dict
  Array[File] known_variants
  Array[File] known_variants_indexes
  String? interval
  
  command {
    ln -s ${reference_fasta} ref.fasta
    ln -s ${reference_fasta_index} ref.fasta.fai
    ln -s ${reference_dict} ref.dict
        
    java -Xmx4g -jar /home/biodocker/bin/GenomeAnalysisTK-2.8-1-g932cd3a/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -R ref.fasta \
    -known ${sep=" -known " known_variants} \
    --disable_auto_index_creation_and_locking_when_reading_rods \
    ${"-L " + interval} \
    -o ${output_filename}
  }
  output {
    File output_intervals = output_filename
  }
  runtime {
    docker: "fzkhan/picard-1.136-gatk-2.8"
  }
}

task GATK_IndelRealigner {
  File reference_fasta
  File reference_fasta_index
  File reference_dict
  File target_intervals
  File input_bam
  File input_bam_index
  String output_bam_filename
  String output_bam_filename_index
  Array[File] known_variants
  Array[File] known_variants_indexes
  
  
  command {
    ln -s ${reference_fasta} ref.fasta
    ln -s ${reference_fasta_index} ref.fasta.fai
    ln -s ${reference_dict} ref.dict
    
    java -Xmx4g -jar /home/biodocker/bin/GenomeAnalysisTK-2.8-1-g932cd3a/GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -R ref.fasta \
    -I ${input_bam} \
    -known ${sep=" -known " known_variants} \
    --disable_auto_index_creation_and_locking_when_reading_rods \
    -targetIntervals ${target_intervals} \
    -o ${output_bam_filename}
  }
  output {
    File output_bam = output_bam_filename
    File output_bam_index = output_bam_filename_index
  }
  runtime {
    docker: "fzkhan/picard-1.136-gatk-2.8"
  }  
}

task GATK_BaseRecalibrator {
  File input_bam
  File input_bam_index
  String output_filename
  File reference_fasta
  File reference_fasta_index
  File reference_dict
  Array[File] known_variants
  Array[File] known_variants_indexes
  
  
  command {
    ln -s ${reference_fasta} ref.fasta
    ln -s ${reference_fasta_index} ref.fasta.fai
    ln -s ${reference_dict} ref.dict
        
    java -Xmx4g -jar /home/biodocker/bin/GenomeAnalysisTK-2.8-1-g932cd3a/GenomeAnalysisTK.jar \
    -T BaseRecalibrator \
    -R ref.fasta \
    -I ${input_bam} \
    -knownSites ${sep=" -knownSites " known_variants} \
    --disable_auto_index_creation_and_locking_when_reading_rods \
    --covariate CycleCovariate --covariate ContextCovariate \
    -o ${output_filename}
  }
  output {
    File output_bqsr_table = output_filename
  }
  runtime {
    docker: "fzkhan/picard-1.136-gatk-2.8"
  }
}

task GATK_PrintReads {
  File reference_fasta
  File reference_fasta_index
  File reference_dict
  File input_bam
  File input_bam_index
  File input_bqsr_table
  String output_bam_filename
  String output_bam_filename_index
  
  command {
    ln -s ${reference_fasta} ref.fasta
    ln -s ${reference_fasta_index} ref.fasta.fai
    ln -s ${reference_dict} ref.dict
    
    java -Xmx4g -jar /home/biodocker/bin/GenomeAnalysisTK-2.8-1-g932cd3a/GenomeAnalysisTK.jar \
    -T PrintReads \
    -R ref.fasta \
    -I ${input_bam} \
    -BQSR ${input_bqsr_table} \
    -o ${output_bam_filename}
  }
  output {
    File output_bam = output_bam_filename
    File output_bam_index = output_bam_filename_index
  }
  runtime {
    docker: "fzkhan/picard-1.136-gatk-2.8"
  }  
}

task GATK_HaplotypeCaller {
  File reference_fasta
  File reference_fasta_index
  File reference_dict
  File input_bam
  File input_bam_index
  File dbsnp_vcf
  File dbsnp_vcf_index
  String? interval
  String output_gvcf_filename
    
  command {
    ln -s ${reference_fasta} ref.fasta
    ln -s ${reference_fasta_index} ref.fasta.fai
    ln -s ${reference_dict} ref.dict
    
    java -Xmx4g -jar /home/biodocker/bin/GenomeAnalysisTK-2.8-1-g932cd3a/GenomeAnalysisTK.jar \
    -T HaplotypeCaller \
    -R ref.fasta \
    -I ${input_bam} \
    --dbsnp ${dbsnp_vcf} \
    ${"-L " + interval} \
    -o ${output_gvcf_filename}
  }
  output {
    File output_gvcf = output_gvcf_filename
  }
  runtime {
    docker: "fzkhan/picard-1.136-gatk-2.8"
  }  
}
workflow GATK_complete_Workflow {
  String basename = "example"
  File ref_fasta
  File ref_fasta_index
  Array[File] reads
  Array[File] known_variants
  Array[File] known_variants_indexes
  File dbsnp_vcf
  File dbsnp_vcf_index
  String interval
  
  call Picard_CreateSequenceDictionary {
    input:
      ref_fasta = ref_fasta
  }
  
  call BWA_mem {
    input:
      reads = reads,
      ref_fasta = ref_fasta,
      output_sam_name = "bwa-mem-2016-08-04.sam"
  }
  
  call samtools_view {
    input:
      input_sam = BWA_mem.output_sam,
      output_filename = "${basename}.unsorted.bam"
  }

  call samtools_sort {
    input:
      input_bam = samtools_view.output_bam,
      output_filename = "${basename}.sorted.bam"
  }
  
  call Picard_MarkDuplicates {
    input:
      input_bam = samtools_sort.output_bam,
      output_filename = "${basename}.sorted.md.bam",
      output_index_filename = "${basename}.sorted.md.bai",
      metrics_filename = "${basename}.markdups.metrics"
  }
  
  call GATK_RealignerTargetCreator {
    input:
      reference_fasta = ref_fasta,
      reference_fasta_index = ref_fasta_index,
      reference_dict = Picard_CreateSequenceDictionary.ref_dict,
      input_bam = Picard_MarkDuplicates.output_bam,
      input_bam_index = Picard_MarkDuplicates.output_bam_index,
      known_variants = known_variants,
      known_variants_indexes = known_variants_indexes,
      interval = interval,
      output_filename = "${basename}.realigner_targets.intervals"
  }

  call GATK_IndelRealigner {
    input:
      reference_fasta = ref_fasta,
      reference_fasta_index = ref_fasta_index,
      reference_dict = Picard_CreateSequenceDictionary.ref_dict,
      input_bam = Picard_MarkDuplicates.output_bam,
      input_bam_index = Picard_MarkDuplicates.output_bam_index,
      target_intervals = GATK_RealignerTargetCreator.output_intervals,
      known_variants = known_variants,      
      known_variants_indexes = known_variants_indexes,
      output_bam_filename = "${basename}.realigned.bam",
      output_bam_filename_index = "${basename}.realigned.bai"
  }
  
  call GATK_BaseRecalibrator {
    input:
      reference_fasta = ref_fasta,
      reference_fasta_index = ref_fasta_index,
      reference_dict = Picard_CreateSequenceDictionary.ref_dict,
      input_bam = GATK_IndelRealigner.output_bam,
      input_bam_index = GATK_IndelRealigner.output_bam_index,
      known_variants = known_variants,
      known_variants_indexes = known_variants_indexes,
      output_filename = "${basename}.bqsr.table"
  }
  
  call GATK_PrintReads {
    input:
      reference_fasta = ref_fasta,
      reference_fasta_index = ref_fasta_index,
      reference_dict = Picard_CreateSequenceDictionary.ref_dict,
      input_bam = GATK_IndelRealigner.output_bam,
      input_bam_index = GATK_IndelRealigner.output_bam_index,
      input_bqsr_table = GATK_BaseRecalibrator.output_bqsr_table,
      output_bam_filename = "${basename}.bam",
      output_bam_filename_index = "${basename}.bai"    
  }
  
  call GATK_HaplotypeCaller {
    input:
      reference_fasta = ref_fasta,
      reference_fasta_index = ref_fasta_index,
      reference_dict = Picard_CreateSequenceDictionary.ref_dict,
      input_bam = GATK_PrintReads.output_bam,
      input_bam_index = GATK_PrintReads.output_bam_index,
      interval = interval,
      dbsnp_vcf = dbsnp_vcf,
      dbsnp_vcf_index = dbsnp_vcf_index,      
      output_gvcf_filename = "${basename}.g.vcf",
  }
  
  output {
    File bam = GATK_PrintReads.output_bam
    File bam_index = GATK_PrintReads.output_bam_index
    File gvcf = GATK_HaplotypeCaller.output_gvcf
  }
  
}