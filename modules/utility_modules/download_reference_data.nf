process DOWNLOAD_REFERENCE_DATA {

  memory 20.GB
  time '1:00:00'
  
  container 'oras://community.wave.seqera.io/library/aria2:1.34.0--a666a8e80d294f96'

  publishDir "${projectDir}/reference_data/", pattern:"GRCm39.fa", mode:'copy', overwrite: false
  publishDir "${projectDir}/reference_data/", pattern:"*.vcf.gz", mode:'copy', overwrite: false
  
  output:
  tuple file("GRCm39.fa"), file("*.vcf.gz"), emit: ref_data
  
  script:
  
  """
    if [ -f ${projectDir}/reference_data/GRCm39.fa ]; then
        echo "Files already exist. Using the existing file."
        cp ${projectDir}/reference_data/GRCm39.fa .
        cp ${projectDir}/reference_data/mgp_REL2021_snps.vcf.gz .
    else
        echo "Downloading reference genome."
        aria2c -c -x 10 -s 10 ${params.ref_genome_url}
        
        echo "Decompressing reference genome."
        gzip -d -c -v Mus_musculus.GRCm39.dna.toplevel.fa.gz > GRCm39.fa

        echo "Downloading Sanger SNPs."
        aria2c -c -x 10 -s 10 ${params.ref_vcf_url}
    fi
  """

  stub:

  """
    touch GRCm39.fa
    touch test.vcf.gz
  """
}