process DOWNLOAD_INDEX_REFERENCE_DATA {

  cpus 8
  memory {60.GB * task.attempt}
  time {3.hour * task.attempt}
  errorStrategy 'retry' 
  maxRetries 1
  
  container 'oras://community.wave.seqera.io/library/bwakit_aria2:834e73b8be9d404a'

  publishDir "${projectDir}/reference_data/", pattern:"GRCm39.fa.*", mode:'copy', overwrite: false
  publishDir "${projectDir}/reference_data/", pattern:"*.vcf.gz", mode:'copy', overwrite: false
  
  output:
  tuple file("GRCm39.fa"), file("*.vcf.gz"), emit: ref_data
  tuple file("GRCm39.fa.*"), file("*.vcf.gz"), emit: ref_genome_indices
  
  script:
  
  """
    if [ -f ${projectDir}/reference_data/GRCm39.fa.* ]; then
        echo "Files already exist. Using the existing file."
        cp ${projectDir}/reference_data/GRCm39* .
        cp ${projectDir}/reference_data/mgp_REL2021_snps.vcf.gz .
    else
        echo "Downloading reference genome."
        aria2c -c -x 10 -s 10 ${params.ref_genome_url}
        
        echo "Decompressing reference genome."
        gzip -d -c -v Mus_musculus.GRCm39.dna.toplevel.fa.gz > GRCm39.fa

        echo "Downloading Sanger SNPs."
        aria2c -c -x 10 -s 10 ${params.ref_vcf_url}

        echo "Indexing reference genome."
        bwa index -a is GRCm39.fa GRCm39.fa

    fi
  """

  stub:

  """
    touch GRCm39.fa
    touch test.vcf.gz
  """
}