python3 transposon_analyzer.py -i GCF_000005005.2_B73_RefGen_v4_genomic.fna -o B73_transposons \
--min-tsd-pattern-size 5 \
--max-tsd-pattern-size 10 \
--gap-size 5000 \
--min-tir-size 10 \
--max-tir-size 20 \
--max-tir-mismatch 2 \
--search-dde-range 30 \
--subterminal-threshold 4.8