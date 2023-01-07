# usr/bin/env bash
# getDomains.sh


# cut the 13th column and sort the names
cut -f13 proteins.tsv | sort | uniq > domains.txt
