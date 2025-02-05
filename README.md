# quilt-nf: A Nextflow pipeline for haplotype reconstruction using low-coverage whole-genome sequencing data

These workflows support the service offered by JAX Genome Technologies for haplotype reconstruction with low-pass WGS. Please contact [Sam Widmayer](mailto:samuel.widmayer@jax.org)  or [Dan Gatti](mailto:dan.gatti@jax.org) for more information.

JAX users are required to have access to the Sumner cluster, and to have Nextflow installed in their home directory. Any setup for external users will require additional support, and those wishing to share these workflows are encouraged to contact the maintainers of this repository.

This pipeline is implemented using [Nextflow](https://www.nextflow.io/), a scalable, reproducible, and increasingly common language used in the development and maintenance of bioinformatics workflows. The modular nature of the workflow is enabled by software containers, such as [Docker](https://www.docker.com/) and [Singularity](https://sylabs.io/singularity), with all the software requirements for executing each step. Specific combinations and versions of software are specified in each container making analyses perfectly reproducible over time as long as the source data is unchanged.
