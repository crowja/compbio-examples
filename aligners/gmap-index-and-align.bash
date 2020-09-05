#! /bin/bash

gmap_home="/mnt/common/gtsw/pkg/gmap/default"

mkdir -p ./my-gmap-dbs

# With fake-genome.fa, create a GMAP index under ./my-gmap-dbs, call it foo-ref.
"${gmap_home}/bin/gmap_build" \
--dir=./my-gmap-dbs \
--genomedb=foo-ref \
fake-genome.fa

# Now align fake-reads.fa to this index.
cat fake-reads.fa \
| "${gmap_home}/bin/gmap" \
--align \
--dir=./my-gmap-dbs \
--db=foo-ref \
> fake-reads-aligned-to-fake-genome-gmap.txt

# Same thing but with gene gff3 output
cat fake-reads.fa \
| "${gmap_home}/bin/gmap" \
--format=gff3_gene \
--dir=./my-gmap-dbs \
--db=foo-ref \
> fake-reads-aligned-to-fake-genome-gmap.gff3


