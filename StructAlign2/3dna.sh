#!/bin/bash
export X3DNA=/usr/local/src/x3dna-v2.1
export PATH=/usr/local/src/x3dna-v2.1/bin:$PATH
cd /var/www/tools/tmp/StructAlign/$2
find_pair /var/www/tools/cgi-bin/$1 out
