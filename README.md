Temp location to build the scripts for creating amptk dbs from ncbi refseq target loci db

requires
 - pytaxonkit
 - biopython
 - requests

Download and install NCBI taxonDB - see https://bioinf.shenwei.me/taxonkit/usage/ 


```
git clone https://github.com/hyphaltip/make_amptkdb_refseq
conda create -n amptkdb_refseq -y -c bioconda -c conda-forge biopython pytaxonkit taxonkit requests
conda activate amptkdb_refseq
cd make_amptkdb_refseq
./install_taxonkit.sh

# download and build NCBI refseq target loci for ITS, 18S, 28S for amptk installation
./fungirefseq_to_amptk.py

# install NCBI ITS refseq
amptk database -i fungi.ITS.fasta  --primer_required none -o FungiITSrefseq_SINTAX --create_db sintax  --install --source NCBI:20231107 --format off
amptk database -i fungi.ITS.fasta  --primer_required none -o FungiITSrefseq --create_db usearch  --install --source NCBI:20231107 --format off
# if you want to use utax and have usearch9
amptk database -i fungi.ITS.fasta  --primer_required none -o FungiITSrefseq_UTAX --create_db utax --install --source NCBI:20231107 --format off 

# install NCBI 18S refseq fungi
amptk database -i fungi.18SrRNA.fasta --primer_required none -o Fungi18Srefseq_SINTAX --create_db sintax --install --source NCBI:20231107 --format off  --fwd_primer None --rev_primer ''
amptk database -i fungi.18SrRNA.fasta --primer_required none -o Fungi18Srefseq --create_db usearch --install --source NCBI:20231107 --format off  --fwd_primer None --rev_primer ''
# if you want to use utax and have usearch9
amptk database -i fungi.18SrRNA.fasta --primer_required none -o Fungi18Srefseq_UTAX --create_db utax --install --source NCBI:20231107 --format off  --fwd_primer None --rev_primer ''


# install NCBI 28S refseq fungi
amptk database -i fungi.28SrRNA.fasta  --primer_required none -o Fungi28Srefseq --create_db usearch --install --source NCBI:20231107 --format off  --fwd_primer None --rev_primer None
amptk database -i fungi.28SrRNA.fasta  --primer_required none -o Fungi28Srefseq_SINTAX --create_db sintax --install --source NCBI:20231107 --format off  --fwd_primer None --rev_primer ''
# if you want to use utax and have usearch9
amptk database -i fungi.28SrRNA.fasta --primer_required none -o Fungi28Srefseq_UTAX --create_db utax --install --source NCBI:20231107 --format off  --fwd_primer None --rev_primer ''

# download and build EukRibo DB - but use NCBI taxonomy so this isn't totally taking on their new taxonomy
./eukribo_to_amptkdb.py
# install EukRibo 
amptk database -i eukribo.fasta  --primer_required none -o EukRibo --create_db usearch  --install --source EukRibo:v2 --format off  --fwd_primer None --rev_primer ''
amptk database -i eukribo.fasta  --primer_required none -o EukRibo_SINTAX --create_db sintax  --install --source EukRibo:v2 --format off  --fwd_primer None --rev_primer ''
# if you want to use utax and have usearch9
amptk database -i eukribo.fasta  --primer_required none -o EukRibo_UTAX --create_db utax  --install --source EukRibo:v2 --format off  --fwd_primer None --rev_primer ''
```
