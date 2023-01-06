#!/usr/bin/env python3
# get_path

import subprocess

kegg = open('ko.txt')
kegg_ids = {}
# Open files for stdout and stderr
out = open('kegg_path.txt','w')
err = open('kegg_path.err','w')
# Iterate over lines of KEGG IDs
for kegg_line in kegg:
    # Remove line terminator
    kegg_line = kegg_line.rstrip()
    # Split line on whitespace
    fields = kegg_line.split()
    if len(fields) > 1:
       kegg = fields[1]
       kegg_ids[kegg] = 1
for kegg_id in kegg_ids:
    print(kegg_id)
    result = subprocess.call("curl http://rest.kegg.jp/link/pathway/" + kegg_id, stdout=out, stderr=err, shell=True)

